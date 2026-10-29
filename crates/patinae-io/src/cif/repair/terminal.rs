//! Restore only missing C-terminal sequence rows with identical label indexing.

use super::identity::{anchored_chains, atom_entity, canonical_id, entry_id, polymers};
use super::*;

#[derive(Debug)]
pub struct TerminalSequenceRepair {
    pub cif: String,
    pub sequence_rows_added: usize,
    pub chain_mapping: Vec<(String, String, String, String)>,
    pub changes: Vec<(usize, usize, String, String)>,
}

/// Explicit metadata-only repair against an independently pinned reference.
/// Target sequences must be exact prefixes of complete reference sequences.
/// Every selected entity needs observed residue anchors, and every matching
/// atom anchor must retain its original label sequence position. Entity IDs
/// need not match. Existing sequence rows and all atom tokens stay untouched.
/// Missing N-termini, internal holes, changed constructs, microheterogeneity,
/// are outside this API. Populated one-letter metadata requires a separate
/// explicit synchronization opt-in and is checked against both monomer tables.
pub fn repair_anchored_terminal_sequences(
    redo: &str,
    reference: &str,
    expected_entry_id: &str,
    entities: &[String],
    synchronize_one_letter: bool,
) -> IoResult<TerminalSequenceRepair> {
    let expected = canonical_id(expected_entry_id);
    if expected.len() != 4
        || !expected.as_bytes()[0].is_ascii_digit()
        || !expected.bytes().all(|b| b.is_ascii_alphanumeric())
    {
        return Err(invalid("expected entry ID must be a four-character PDB ID"));
    }
    let target = tables(redo)?;
    let source = tables(reference)?;
    if canonical_id(entry_id(&target)?) != expected || canonical_id(entry_id(&source)?) != expected
    {
        return Err(invalid(
            "terminal repair requires matching expected entry IDs",
        ));
    }
    super::super::read_cif_str(redo)?;
    super::super::read_cif_str(reference)?;
    let selected: BTreeSet<_> = entities.iter().cloned().collect();
    if selected.is_empty()
        || selected.len() != entities.len()
        || !selected.is_subset(&polymers(&target)?)
    {
        return Err(invalid(
            "select unique observed polymer entities explicitly",
        ));
    }
    let mapping = anchored_chains(&target, &source, &selected)?;
    let old_table = table(&target, "_entity_poly_seq")?;
    let ref_table = table(&source, "_entity_poly_seq")?;
    if !old_table.is_loop {
        return Err(invalid("terminal repair requires a sequence loop"));
    }
    let old = sequence(old_table)?;
    let full = sequence(ref_table)?;
    let atoms = table(&target, "_atom_site")?;
    let ref_atoms = table(&source, "_atom_site")?;
    let source_polymers = polymers(&source)?;
    // Unlike type anchoring, terminal repair also proves label positions.
    let mut anchors: BTreeMap<Vec<String>, BTreeSet<(String, String, u32)>> = BTreeMap::new();
    for row in &ref_atoms.rows {
        let entity = atom_entity(&source, ref_atoms, row)?;
        if !source_polymers.contains(entity) {
            continue;
        }
        let position = ref_atoms
            .value(row, "_atom_site.label_seq_id")?
            .parse::<u32>()
            .map_err(|_| invalid("invalid reference label sequence position"))?;
        anchors
            .entry(atom_key(ref_atoms, row)?)
            .or_default()
            .insert((
                ref_atoms.value(row, "_atom_site.label_asym_id")?.to_owned(),
                entity.to_owned(),
                position,
            ));
    }
    for row in &atoms.rows {
        let entity = atom_entity(&target, atoms, row)?;
        if !selected.contains(entity) {
            continue;
        }
        let chain = atoms.value(row, "_atom_site.label_asym_id")?;
        let position = atoms
            .value(row, "_atom_site.label_seq_id")?
            .parse::<u32>()
            .map_err(|_| invalid("invalid target label sequence position"))?;
        if old.get(&(entity.to_owned(), position)).map(String::as_str)
            != Some(atoms.value(row, "_atom_site.label_comp_id")?)
        {
            return Err(invalid(
                "observed component differs from target sequence metadata",
            ));
        }
        if let Some(matches) = anchors.get(&atom_key(atoms, row)?) {
            let mapped = mapping
                .iter()
                .find(|m| m.0 == chain && m.1 == entity)
                .ok_or_else(|| invalid("missing selected chain anchor"))?;
            if matches
                .iter()
                .any(|(c, e, p)| c != &mapped.2 || e != &mapped.3 || *p != position)
            {
                return Err(invalid(
                    "terminal repair refuses shifted/ambiguous label sequence anchors",
                ));
            }
        }
    }
    let mut added = String::new();
    let mut edits: Vec<(Range<usize>, String)> = Vec::new();
    let mut rows_added = 0;
    for entity in &selected {
        let refs: BTreeSet<_> = mapping
            .iter()
            .filter(|m| &m.1 == entity)
            .map(|m| &m.3)
            .collect();
        if refs.len() != 1 {
            return Err(invalid(
                "selected entity maps to multiple reference entities",
            ));
        }
        let ref_entity = refs.into_iter().next().unwrap();
        let declared: BTreeMap<_, _> = old
            .iter()
            .filter(|((e, _), _)| e == entity)
            .map(|((_, p), m)| (*p, m))
            .collect();
        let reference_seq: BTreeMap<_, _> = full
            .iter()
            .filter(|((e, _), _)| e == ref_entity)
            .map(|((_, p), m)| (*p, m))
            .collect();
        if declared.is_empty()
            || reference_seq.len() <= declared.len()
            || declared.keys().copied().ne(1..=declared.len() as u32)
            || reference_seq
                .keys()
                .copied()
                .ne(1..=reference_seq.len() as u32)
            || declared
                .iter()
                .any(|(p, m)| reference_seq.get(p) != Some(m))
        {
            return Err(invalid(
                "target must be an exact numbered prefix of a longer complete reference",
            ));
        }
        // Do not leave another, stale sequence representation behind.
        if let Some(poly) = target.get("_entity_poly") {
            if poly
                .rows
                .iter()
                .filter(|r| {
                    poly.value(r, "_entity_poly.entity_id")
                        .is_ok_and(|e| e == entity)
                })
                .count()
                > 1
            {
                return Err(invalid("duplicate target one-letter sequence entity"));
            }
            for row in &poly.rows {
                if poly.value(row, "_entity_poly.entity_id")? != entity {
                    continue;
                }
                for field in [
                    "_entity_poly.pdbx_seq_one_letter_code",
                    "_entity_poly.pdbx_seq_one_letter_code_can",
                ] {
                    if let Ok(i) = poly.col(field) {
                        if !missing(&row[i].value) {
                            if !synchronize_one_letter {
                                return Err(invalid("selected entity has one-letter sequence metadata; explicit synchronization required"));
                            }
                            let ref_poly = table(&source, "_entity_poly")?;
                            let ref_rows = ref_poly
                                .rows
                                .iter()
                                .filter(|r| {
                                    ref_poly
                                        .value(r, "_entity_poly.entity_id")
                                        .is_ok_and(|e| e == ref_entity)
                                })
                                .collect::<Vec<_>>();
                            if ref_rows.len() != 1 {
                                return Err(invalid(
                                    "ambiguous/missing reference one-letter entity",
                                ));
                            }
                            let ref_row = ref_rows[0];
                            let kind = patinae_mol::PolymerKind::from_mmcif(
                                ref_poly.value(ref_row, "_entity_poly.type")?,
                            );
                            let primary =
                                &ref_row[ref_poly.col("_entity_poly.pdbx_seq_one_letter_code")?];
                            let decoded =
                                crate::entity::one_letter_sequence(&primary.value, Some(&kind))?;
                            if decoded
                                .iter()
                                .ne(reference_seq.values().map(|m| m.as_str()))
                            {
                                return Err(invalid("reference one-letter sequence disagrees with full monomer metadata"));
                            }
                            let replacement = &ref_row[ref_poly.col(field)?];
                            if field.ends_with("_can") {
                                // Canonical strings are symbols, not CCD identities. Use
                                // the existing decoder as a symbol tokenizer: no parent
                                // chemistry or D/L equivalence is inferred here.
                                let symbols_kind = patinae_mol::PolymerKind::PeptideL;
                                let before = crate::entity::one_letter_sequence(
                                    &row[i].value,
                                    Some(&symbols_kind),
                                )?;
                                let after = crate::entity::one_letter_sequence(
                                    &replacement.value,
                                    Some(&symbols_kind),
                                )?;
                                let nucleotide = matches!(
                                    kind,
                                    patinae_mol::PolymerKind::Dna
                                        | patinae_mol::PolymerKind::Rna
                                        | patinae_mol::PolymerKind::DnaRnaHybrid
                                );
                                if before.len() != declared.len()
                                    || after.len() != reference_seq.len()
                                    || before.iter().zip(&after).any(|(a, b)| {
                                        a != b && a != "UNK" && !(nucleotide && a == "ASN")
                                    })
                                {
                                    return Err(invalid(
                                        "canonical one-letter prefix disagrees with reference",
                                    ));
                                }
                                for (mon, symbol) in decoded.iter().zip(&after) {
                                    if (mon == "UNK" && symbol != "UNK")
                                        || (nucleotide
                                            && matches!(mon.as_str(), "N" | "DN")
                                            && !matches!(symbol.as_str(), "UNK" | "ASN"))
                                    {
                                        return Err(invalid("canonical reference claims known identity for an unknown monomer"));
                                    }
                                }
                            } else {
                                let before =
                                    crate::entity::one_letter_sequence(&row[i].value, Some(&kind))?;
                                if before.iter().ne(declared.values().map(|m| m.as_str())) {
                                    return Err(invalid("target one-letter sequence disagrees with monomer metadata"));
                                }
                            }
                            edits.push((
                                row[i].span.clone(),
                                format!("\n{}\n", &reference[replacement.span.clone()]),
                            ));
                        }
                    }
                }
            }
        }
        if !entity
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_')
        {
            return Err(invalid("unsafe entity token for sequence insertion"));
        }
        for row in &ref_table.rows {
            if ref_table.value(row, "_entity_poly_seq.entity_id")? != ref_entity {
                continue;
            }
            let position = ref_table
                .value(row, "_entity_poly_seq.num")?
                .parse::<u32>()
                .map_err(|_| invalid("invalid sequence position"))?;
            if position <= declared.len() as u32 {
                continue;
            }
            let mut values = Vec::new();
            for column in &old_table.columns {
                if column == "_entity_poly_seq.entity_id" {
                    values.push(entity.clone());
                } else if let Ok(i) = ref_table.col(column) {
                    values.push(reference[row[i].span.clone()].to_owned());
                } else {
                    values.push("?".into());
                }
            }
            added.push_str(&format!("\n{}", values.join(" ")));
            rows_added += 1;
        }
    }
    let offset = old_table.span.end;
    edits.push((offset..offset, added));
    edits.sort_by_key(|(s, _)| (s.start, s.end));
    for pair in edits.windows(2) {
        if pair[0].0.end > pair[1].0.start {
            return Err(invalid("overlapping terminal metadata edits"));
        }
    }
    let changes = edits
        .iter()
        .map(|(s, v)| (s.start, s.end, redo[s.clone()].to_owned(), v.clone()))
        .collect();
    let mut cif = redo.to_owned();
    for (span, value) in edits.into_iter().rev() {
        cif.replace_range(span, &value);
    }
    let after_tables = tables(&cif)?;
    let after = table(&after_tables, "_atom_site")?;
    if atoms.columns != after.columns || atoms.rows.len() != after.rows.len() {
        return Err(invalid("terminal repair changed atom layout"));
    }
    for (a, b) in atoms.rows.iter().zip(&after.rows) {
        for (old, new) in a.iter().zip(b) {
            if redo[old.span.clone()] != cif[new.span.clone()] {
                return Err(invalid("terminal repair changed protected atom token"));
            }
        }
    }
    super::super::read_cif_str(&cif)?;
    Ok(TerminalSequenceRepair {
        cif,
        sequence_rows_added: rows_added,
        chain_mapping: mapping,
        changes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn doc(entity: &str, sequence: &str, position: u32) -> String {
        format!("data_1abc\n_entry.id 1abc\n_entity.id {entity}\n_entity.type polymer\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n{sequence}\n#\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 C CA GLY X {entity} {position} A 10 ? 1.234500 2.000 3.0\n")
    }
    #[test]
    fn restores_suffix_with_independent_entities_and_unchanged_atoms() {
        let target = doc("5", "5 1 GLY", 1);
        let reference = doc("3", "3 1 GLY\n3 2 ALA\n3 3 MET", 1);
        let result =
            repair_anchored_terminal_sequences(&target, &reference, "1abc", &["5".into()], false)
                .unwrap();
        assert_eq!(result.sequence_rows_added, 2);
        assert!(result.cif.contains("5 1 GLY\n5 2 ALA\n5 3 MET"));
        assert!(result.cif.contains("1.234500 2.000 3.0"));
        let mut replay = target.clone();
        for (start, end, old, new) in result.changes.iter().rev() {
            assert_eq!(&target[*start..*end], old);
            replay.replace_range(*start..*end, new);
        }
        assert_eq!(replay, result.cif);
    }
    #[test]
    fn refuses_shift_conflict_wrong_identity_and_one_letter_metadata() {
        let target = doc("5", "5 1 GLY", 1);
        let reference = doc("3", "3 1 GLY\n3 2 GLY", 2);
        assert!(repair_anchored_terminal_sequences(
            &target,
            &reference,
            "1abc",
            &["5".into()],
            false
        )
        .is_err());
        let reference = doc("3", "3 1 GLY\n3 2 ALA", 1);
        assert!(repair_anchored_terminal_sequences(
            &target.replace("5 1 GLY", "5 1 ALA"),
            &reference,
            "1abc",
            &["5".into()],
            false
        )
        .is_err());
        assert!(repair_anchored_terminal_sequences(
            &target,
            &reference,
            "2abc",
            &["5".into()],
            false
        )
        .is_err());
        let target = format!(
            "{target}\n_entity_poly.entity_id 5\n_entity_poly.pdbx_seq_one_letter_code G\n"
        );
        assert!(repair_anchored_terminal_sequences(
            &target,
            &reference,
            "1abc",
            &["5".into()],
            false
        )
        .is_err());
    }

    #[test]
    fn synchronizes_one_letter_metadata_only_with_explicit_opt_in() {
        let target = format!("{}\n_entity_poly.entity_id 5\n_entity_poly.type 'polypeptide(L)'\n_entity_poly.pdbx_seq_one_letter_code G\n_entity_poly.pdbx_seq_one_letter_code_can G\n",doc("5","5 1 GLY",1));
        let reference = format!("{}\n_entity_poly.entity_id 3\n_entity_poly.type 'polypeptide(L)'\n_entity_poly.pdbx_seq_one_letter_code\n;GA\n;\n_entity_poly.pdbx_seq_one_letter_code_can GA\n",doc("3","3 1 GLY\n3 2 ALA",1));
        let result =
            repair_anchored_terminal_sequences(&target, &reference, "1abc", &["5".into()], true)
                .unwrap();
        assert_eq!(result.changes.len(), 3);
        assert_eq!(result.sequence_rows_added, 1);
        let mut replay = target.clone();
        for (start, end, old, new) in result.changes.iter().rev() {
            assert_eq!(&target[*start..*end], old);
            replay.replace_range(*start..*end, new);
        }
        assert_eq!(replay, result.cif);
        assert!(repair_anchored_terminal_sequences(
            &target,
            &reference.replace(";GA", ";AA"),
            "1abc",
            &["5".into()],
            true
        )
        .is_err());
        assert!(repair_anchored_terminal_sequences(
            &target.replace("code G\n", "code A\n"),
            &reference,
            "1abc",
            &["5".into()],
            true
        )
        .is_err());
    }

    #[test]
    fn canonical_mask_is_not_confused_with_primary_chemical_identity() {
        let target = format!("{}\n_entity_poly.entity_id 5\n_entity_poly.pdbx_seq_one_letter_code '(DAL)'\n_entity_poly.pdbx_seq_one_letter_code_can X\n",doc("5","5 1 DAL",1).replace("CA GLY","CA DAL"));
        let reference = format!("{}\n_entity_poly.entity_id 3\n_entity_poly.type 'polypeptide(L)'\n_entity_poly.pdbx_seq_one_letter_code '(DAL)(ZGL)'\n_entity_poly.pdbx_seq_one_letter_code_can AX\n",doc("3","3 1 DAL\n3 2 ZGL",1).replace("CA GLY","CA DAL"));
        assert!(repair_anchored_terminal_sequences(
            &target,
            &reference,
            "1abc",
            &["5".into()],
            true
        )
        .is_ok());
        assert!(repair_anchored_terminal_sequences(
            &target.replace("code_can X", "code_can G"),
            &reference,
            "1abc",
            &["5".into()],
            true
        )
        .is_err());
    }
}
