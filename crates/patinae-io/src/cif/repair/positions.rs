//! Explicit reference sequence view, never a rewrite of source atom labels.

use super::identity::{atom_entity, canonical_id, entry_id, polymers};
use super::*;

#[derive(Debug)]
pub struct ChainSequenceMap {
    pub target_asym_id: String,
    pub reference_asym_id: String,
    pub reference_entity_id: String,
    pub polymer_type: String,
    pub sequence: Vec<String>,
}

#[derive(Debug)]
pub struct ResiduePositionMap {
    pub target_asym_id: String,
    pub target_entity_id: String,
    pub target_label_seq_id: String,
    pub auth_asym_id: String,
    pub auth_seq_id: String,
    pub insertion_code: String,
    pub monomer: String,
    pub reference_asym_id: String,
    pub reference_entity_id: String,
    pub reference_seq_id: u32,
}

#[derive(Debug)]
pub struct PolymerPositionMap {
    pub chains: Vec<ChainSequenceMap>,
    pub residues: Vec<ResiduePositionMap>,
}

/// Build a separate, strict position map for explicitly selected label chains.
/// Each observed residue requires exactly one reference chain/entity/position,
/// consistent across its matching atom keys. Source labels (including zero or
/// missing label positions) are kept literally. No sequence alignment, chemical
/// equivalence, coordinate comparison, or automatic reader fallback is used.
pub fn map_polymer_positions(
    target_text: &str,
    reference_text: &str,
    expected_entry_id: &str,
    target_chains: &[String],
) -> IoResult<PolymerPositionMap> {
    let expected = canonical_id(expected_entry_id);
    if expected.len() != 4
        || !expected.as_bytes()[0].is_ascii_digit()
        || !expected.bytes().all(|b| b.is_ascii_alphanumeric())
    {
        return Err(invalid(
            "position mapping requires a four-character expected PDB ID",
        ));
    }
    let target = tables(target_text)?;
    let source = tables(reference_text)?;
    if canonical_id(entry_id(&target)?) != expected || canonical_id(entry_id(&source)?) != expected
    {
        return Err(invalid(
            "position mapping entry IDs differ from expected identity",
        ));
    }
    super::super::read_cif_str(target_text)?;
    super::super::read_cif_str(reference_text)?;
    let selected: BTreeSet<_> = target_chains.iter().cloned().collect();
    if selected.is_empty()
        || selected.len() != target_chains.len()
        || selected.iter().any(|c| missing(c))
    {
        return Err(invalid(
            "position mapping requires unique explicit nonempty target chains",
        ));
    }
    let atoms = table(&target, "_atom_site")?;
    let ref_atoms = table(&source, "_atom_site")?;
    let ref_polymers = polymers(&source)?;
    let sequences = sequence(table(&source, "_entity_poly_seq")?)?;
    type ReferencePosition = (String, String, u32);
    let mut index: BTreeMap<Vec<String>, BTreeSet<ReferencePosition>> = BTreeMap::new();
    for row in &ref_atoms.rows {
        let entity = atom_entity(&source, ref_atoms, row)?;
        if !ref_polymers.contains(entity) {
            continue;
        }
        let chain = ref_atoms.value(row, "_atom_site.label_asym_id")?;
        let position = ref_atoms
            .value(row, "_atom_site.label_seq_id")?
            .parse::<u32>()
            .ok()
            .filter(|n| *n > 0)
            .ok_or_else(|| invalid("reference polymer has invalid label sequence position"))?;
        if missing(chain) {
            return Err(invalid("reference polymer chain missing"));
        }
        let key = atom_key(ref_atoms, row)?;
        if key[..4].iter().any(|v| v.is_empty()) {
            return Err(invalid("reference atom lacks author/chemical anchor"));
        }
        index
            .entry(key)
            .or_default()
            .insert((chain.to_owned(), entity.to_owned(), position));
    }
    // source chain, raw entity, raw sequence label, author chain/number,
    // insertion code and exact CCD component. No numeric labels are normalized.
    type SourceResidue = (String, String, String, String, String, String, String);
    let mut grouped: BTreeMap<SourceResidue, BTreeSet<ReferencePosition>> = BTreeMap::new();
    for row in &atoms.rows {
        let chain = atoms.value(row, "_atom_site.label_asym_id")?;
        if !selected.contains(chain) {
            continue;
        }
        let key = atom_key(atoms, row)?;
        if key[..4].iter().any(|v| v.is_empty()) {
            return Err(invalid("target atom lacks author/chemical anchor"));
        }
        let residue = (
            chain.to_owned(),
            atoms.value(row, "_atom_site.label_entity_id")?.to_owned(),
            atoms.value(row, "_atom_site.label_seq_id")?.to_owned(),
            key[0].clone(),
            key[1].clone(),
            key[4].clone(),
            key[2].clone(),
        );
        let matches = grouped.entry(residue).or_default();
        if let Some(candidates) = index.get(&key) {
            matches.extend(candidates.iter().cloned());
        }
    }
    let mut chain_mapping: BTreeMap<String, (String, String)> = BTreeMap::new();
    let mut occupied = BTreeSet::new();
    let mut residues = Vec::new();
    for ((chain, entity, label, auth, number, insertion, component), candidates) in grouped {
        if candidates.len() != 1 {
            return Err(invalid(format!("position mapping refuses unanchored/ambiguous residue {chain}/{auth}:{number}{insertion}/{component} ({} positions)",candidates.len())));
        }
        let (ref_chain, ref_entity, ref_position) = candidates.into_iter().next().unwrap();
        if sequences.get(&(ref_entity.clone(), ref_position)) != Some(&component) {
            return Err(invalid(format!("reference observed/declared component disagreement at {ref_chain}:{ref_position}/{component}")));
        }
        let mapped = (ref_chain.clone(), ref_entity.clone());
        if let Some(previous) = chain_mapping.insert(chain.clone(), mapped.clone()) {
            if previous != mapped {
                return Err(invalid(
                    "source chain maps to multiple reference chains/entities",
                ));
            }
        }
        // Allow chain fragments to map to one reference chain, but never let
        // distinct source residues claim the same canonical chain position.
        if !occupied.insert((ref_chain.clone(), ref_entity.clone(), ref_position)) {
            return Err(invalid(
                "multiple source residues claim one reference chain position",
            ));
        }
        residues.push(ResiduePositionMap {
            target_asym_id: chain,
            target_entity_id: entity,
            target_label_seq_id: label,
            auth_asym_id: auth,
            auth_seq_id: number,
            insertion_code: insertion,
            monomer: component,
            reference_asym_id: ref_chain,
            reference_entity_id: ref_entity,
            reference_seq_id: ref_position,
        });
    }
    if chain_mapping.keys().cloned().collect::<BTreeSet<_>>() != selected {
        return Err(invalid(
            "selected target chain has no observed atom anchors",
        ));
    }
    let poly = table(&source, "_entity_poly")?;
    let mut chains = Vec::new();
    for (chain, (ref_chain, entity)) in chain_mapping {
        let seq: BTreeMap<_, _> = sequences
            .iter()
            .filter(|((e, _), _)| e == &entity)
            .map(|((_, p), m)| (*p, m.clone()))
            .collect();
        if seq.is_empty() || seq.keys().copied().ne(1..=seq.len() as u32) {
            return Err(invalid(
                "position mapping requires complete reference sequence metadata",
            ));
        }
        let rows = poly
            .rows
            .iter()
            .filter(|r| {
                poly.value(r, "_entity_poly.entity_id")
                    .is_ok_and(|id| id == entity)
            })
            .collect::<Vec<_>>();
        if rows.len() != 1 {
            return Err(invalid("ambiguous/missing reference polymer type row"));
        }
        let kind = poly.value(rows[0], "_entity_poly.type")?;
        if missing(kind) || kind == "other" {
            return Err(invalid("reference polymer type uninformative"));
        }
        chains.push(ChainSequenceMap {
            target_asym_id: chain,
            reference_asym_id: ref_chain,
            reference_entity_id: entity,
            polymer_type: kind.to_owned(),
            sequence: seq.into_values().collect(),
        });
    }
    Ok(PolymerPositionMap { chains, residues })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn doc(entity: &str, sequence: &str, atoms: &str) -> String {
        format!("data_1abc\n_entry.id 1abc\n_entity.id {entity}\n_entity.type polymer\n_entity_poly.entity_id {entity}\n_entity_poly.type 'polypeptide(L)'\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n{sequence}\n#\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n{atoms}\n")
    }
    #[test]
    fn keeps_raw_labels_and_maps_independent_positions_without_coordinates() {
        let target = doc(
            "5",
            "5 1 ALA\n5 2 GLY",
            "1 C CA GLY X 5 1 A 10 ? 1.234500 2.000 3.0",
        );
        let reference = doc(
            "3",
            "3 1 ALA\n3 2 GLY",
            "1 C CA GLY Q 3 2 A 10 ? 90.0 80.0 70.0",
        );
        let result = map_polymer_positions(&target, &reference, "1abc", &["X".into()]).unwrap();
        assert_eq!(result.residues[0].target_label_seq_id, "1");
        assert_eq!(result.residues[0].reference_seq_id, 2);
        assert_eq!(result.chains[0].reference_entity_id, "3");
        assert_eq!(result.chains[0].sequence, vec!["ALA", "GLY"]);
        assert!(target.contains("1.234500 2.000 3.0"));
    }
    #[test]
    fn rejects_wrong_identity_ambiguous_and_modified_component_anchors() {
        let target = doc("5", "5 1 GLY", "1 C CA GLY X 5 1 A 10 ? 1.0 2.0 3.0");
        let reference = doc("3", "3 1 GLY", "1 C CA GLY Q 3 1 A 10 ? 1.0 2.0 3.0");
        assert!(map_polymer_positions(&target, &reference, "2abc", &["X".into()]).is_err());
        assert!(map_polymer_positions(&target, &reference, "1abc", &["Y".into()]).is_err());
        let ambiguous = format!("{reference}2 C CA GLY R 3 1 A 10 ? 1.0 2.0 3.0\n");
        assert!(map_polymer_positions(&target, &ambiguous, "1abc", &["X".into()]).is_err());
        assert!(map_polymer_positions(
            &target.replace("GLY", "DAL"),
            &reference,
            "1abc",
            &["X".into()]
        )
        .is_err());
    }
    #[test]
    fn permits_explicit_fragment_mapping_but_rejects_position_collisions() {
        let target = doc(
            "5",
            "5 1 ALA\n5 2 GLY",
            "1 C CA ALA X 5 1 A 10 ? 1.0 2.0 3.0\n2 C CA GLY Y ? 0 A 11 ? 1.0 2.0 3.0",
        );
        let reference = doc(
            "3",
            "3 1 ALA\n3 2 GLY",
            "1 C CA ALA Q 3 1 A 10 ? 1.0 2.0 3.0\n2 C CA GLY Q 3 2 A 11 ? 1.0 2.0 3.0",
        );
        let result =
            map_polymer_positions(&target, &reference, "1abc", &["X".into(), "Y".into()]).unwrap();
        assert_eq!(result.chains.len(), 2);
        assert_eq!(
            result
                .residues
                .iter()
                .find(|r| r.target_asym_id == "Y")
                .unwrap()
                .target_label_seq_id,
            "0"
        );
        let collision = format!("{target}3 C CA GLY Z ? 0 A 11 ? 1.0 2.0 3.0\n");
        assert!(map_polymer_positions(
            &collision,
            &reference,
            "1abc",
            &["X".into(), "Y".into(), "Z".into()]
        )
        .is_err());
    }
}
