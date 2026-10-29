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
/// and simultaneous one-letter sequence editing are outside this API.
pub fn repair_anchored_terminal_sequences(
    redo: &str,
    reference: &str,
    expected_entry_id: &str,
    entities: &[String],
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
                            return Err(invalid("selected entity has one-letter sequence metadata; explicit synchronization required"));
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
    let mut cif = redo.to_owned();
    cif.insert_str(offset, &added);
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
        changes: vec![(offset, offset, String::new(), added)],
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
            repair_anchored_terminal_sequences(&target, &reference, "1abc", &["5".into()]).unwrap();
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
        assert!(
            repair_anchored_terminal_sequences(&target, &reference, "1abc", &["5".into()]).is_err()
        );
        let reference = doc("3", "3 1 GLY\n3 2 ALA", 1);
        assert!(repair_anchored_terminal_sequences(
            &target.replace("5 1 GLY", "5 1 ALA"),
            &reference,
            "1abc",
            &["5".into()]
        )
        .is_err());
        assert!(
            repair_anchored_terminal_sequences(&target, &reference, "2abc", &["5".into()]).is_err()
        );
        let target = format!(
            "{target}\n_entity_poly.entity_id 5\n_entity_poly.pdbx_seq_one_letter_code G\n"
        );
        assert!(
            repair_anchored_terminal_sequences(&target, &reference, "1abc", &["5".into()]).is_err()
        );
    }
}
