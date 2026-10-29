//! Opt-in metadata edits with author-residue/atom anchors, never entity-ID joins.

use super::*;

#[derive(Debug)]
pub struct IdentityTypeRepair {
    pub cif: String,
    /// Target chain, target entity, reference chain, reference entity.
    pub chain_mapping: Vec<(String, String, String, String)>,
    pub changes: Vec<(usize, usize, String, String)>,
}

fn canonical_id(value: &str) -> String {
    let value = value.to_ascii_lowercase();
    value.strip_prefix("pdb_0000").unwrap_or(&value).to_owned()
}

fn entry_id(tables: &HashMap<String, Table>) -> IoResult<&str> {
    let entry = table(tables, "_entry")?;
    if entry.rows.len() != 1 {
        return Err(invalid("identity repair requires one _entry.id"));
    }
    entry.value(&entry.rows[0], "_entry.id")
}

fn polymers(tables: &HashMap<String, Table>) -> IoResult<BTreeSet<String>> {
    let entities = table(tables, "_entity")?;
    let mut seen = BTreeSet::new();
    let mut result = BTreeSet::new();
    for row in &entities.rows {
        let id = entities.value(row, "_entity.id")?;
        if missing(id) || !seen.insert(id.to_owned()) {
            return Err(invalid("missing/duplicate entity ID"));
        }
        if entities.value(row, "_entity.type")? == "polymer" {
            result.insert(id.to_owned());
        }
    }
    Ok(result)
}

fn atom_entity<'a>(
    tables: &'a HashMap<String, Table>,
    atoms: &Table,
    row: &'a [Cell],
) -> IoResult<&'a str> {
    let raw = atoms.value(row, "_atom_site.label_entity_id")?;
    let chain = atoms.value(row, "_atom_site.label_asym_id")?;
    let mut mapped = None;
    if let Some(asym) = tables.get("_struct_asym") {
        for entry in &asym.rows {
            if asym.value(entry, "_struct_asym.id")? != chain {
                continue;
            }
            let entity = asym.value(entry, "_struct_asym.entity_id")?;
            if missing(entity) {
                continue;
            }
            if mapped.is_some_and(|previous| previous != entity) {
                return Err(invalid("conflicting struct_asym entity anchors"));
            }
            mapped = Some(entity);
        }
    }
    if !missing(raw) {
        if mapped.is_some_and(|entity| entity != raw) {
            return Err(invalid("atom/struct_asym entity conflict"));
        }
        Ok(raw)
    } else {
        Ok(mapped.unwrap_or(raw))
    }
}

/// Prove every observed polymer residue using at least one exact atom key.
/// Alternate conformers/models may repeat keys only within the same chain/entity.
/// No coordinates, sequence alignment, or matching entity numbers are used.
fn anchored_chains(
    target: &HashMap<String, Table>,
    source: &HashMap<String, Table>,
    selected: &BTreeSet<String>,
) -> IoResult<Vec<(String, String, String, String)>> {
    let atoms = table(target, "_atom_site")?;
    let reference = table(source, "_atom_site")?;
    let source_polymers = polymers(source)?;
    let mut index: BTreeMap<Vec<String>, BTreeSet<(String, String)>> = BTreeMap::new();
    for row in &reference.rows {
        let entity = atom_entity(source, reference, row)?;
        if !source_polymers.contains(entity) {
            continue;
        }
        let chain = reference.value(row, "_atom_site.label_asym_id")?;
        if missing(chain) {
            return Err(invalid("reference polymer chain is missing"));
        }
        index
            .entry(atom_key(reference, row)?)
            .or_default()
            .insert((chain.to_owned(), entity.to_owned()));
    }
    // Per residue, collect candidates from matching atoms. Conflicting atom
    // anchors are rejected rather than resolved by majority vote.
    type ResidueAnchor = (String, String, String, String, String, String);
    let mut residues: BTreeMap<ResidueAnchor, BTreeSet<(String, String)>> = BTreeMap::new();
    for row in &atoms.rows {
        let entity = atom_entity(target, atoms, row)?;
        if !selected.contains(entity) {
            continue;
        }
        let chain = atoms.value(row, "_atom_site.label_asym_id")?;
        let key = atom_key(atoms, row)?;
        if missing(chain)
            || key[0].is_empty()
            || key[1].is_empty()
            || key[2].is_empty()
            || key[3].is_empty()
        {
            return Err(invalid(
                "polymer atom lacks author/chain/component/atom anchor",
            ));
        }
        let residue = (
            chain.to_owned(),
            entity.to_owned(),
            key[0].clone(),
            key[1].clone(),
            key[2].clone(),
            key[4].clone(),
        );
        let candidates = residues.entry(residue).or_default();
        if let Some(matches) = index.get(&key) {
            // Union over atoms; differing candidate chains remain ambiguous.
            candidates.extend(matches.iter().cloned());
        }
    }
    let mut chains: BTreeMap<(String, String), (String, String)> = BTreeMap::new();
    let mut observed = BTreeSet::new();
    for ((chain, entity, auth, number, component, insertion), candidates) in residues {
        if candidates.len() != 1 {
            return Err(invalid(format!("unanchored/ambiguous residue {chain}/{auth}:{number}{insertion}/{component} ({} reference chains)", candidates.len())));
        }
        let mapped = candidates.into_iter().next().unwrap();
        if let Some(previous) = chains.insert((chain.clone(), entity.clone()), mapped.clone()) {
            if previous != mapped {
                return Err(invalid(format!(
                    "chain {chain} maps to multiple reference chains/entities"
                )));
            }
        }
        observed.insert(entity);
    }
    if observed != *selected || chains.is_empty() {
        return Err(invalid(
            "selected polymer entity has no observed atom anchors",
        ));
    }
    let mut chain_entities = BTreeMap::new();
    for (chain, entity) in chains.keys() {
        if let Some(old) = chain_entities.insert(chain, entity) {
            if old != entity {
                return Err(invalid("target chain belongs to multiple entities"));
            }
        }
    }
    Ok(chains
        .into_iter()
        .map(|((a, b), (c, d))| (a, b, c, d))
        .collect())
}

fn polymer_type<'a>(
    tables: &'a HashMap<String, Table>,
    entity: &str,
) -> IoResult<Option<&'a Cell>> {
    let Some(poly) = tables.get("_entity_poly") else {
        return Ok(None);
    };
    let mut rows = poly.rows.iter().filter(|row| {
        poly.value(row, "_entity_poly.entity_id")
            .is_ok_and(|id| id == entity)
    });
    let first = rows.next();
    if rows.next().is_some() {
        return Err(invalid("duplicate polymer type entity"));
    }
    first
        .map(|row| poly.col("_entity_poly.type").map(|i| &row[i]))
        .transpose()
}

/// Repair only the explicitly named missing/`other` polymer types and/or entry
/// identity. The caller must pin both inputs independently. All observed
/// polymer residues are anchored for identity repair; type repair anchors the
/// selected entities and additionally requires their full declared sequences
/// and observed label-sequence assignments to agree. No sequence is replaced.
pub fn repair_identity_and_polymer_types(
    redo: &str,
    reference: &str,
    expected_entry_id: &str,
    identity: bool,
    type_entities: &[String],
) -> IoResult<IdentityTypeRepair> {
    let expected = canonical_id(expected_entry_id);
    if expected.len() != 4
        || !expected.as_bytes()[0].is_ascii_digit()
        || !expected.bytes().all(|b| b.is_ascii_alphanumeric())
    {
        return Err(invalid("expected entry ID must be a four-character PDB ID"));
    }
    let target = tables(redo)?;
    let source = tables(reference)?;
    if canonical_id(entry_id(&source)?) != expected {
        return Err(invalid("reference entry ID differs from expected ID"));
    }
    super::super::read_cif_str(reference)?;
    if !identity && canonical_id(entry_id(&target)?) != expected {
        return Err(invalid(
            "target entry ID differs; explicit identity repair required",
        ));
    }
    let target_polymers = polymers(&target)?;
    let selected: BTreeSet<_> = type_entities.iter().cloned().collect();
    if selected.len() != type_entities.len() || !selected.is_subset(&target_polymers) {
        return Err(invalid("duplicate/non-polymer selected type entity"));
    }
    if !identity && selected.is_empty() {
        return Err(invalid("no identity/type repair requested"));
    }
    let mapping = anchored_chains(
        &target,
        &source,
        if identity {
            &target_polymers
        } else {
            &selected
        },
    )?;
    let mut edits: Vec<(Range<usize>, String)> = Vec::new();
    if identity {
        let blocks = tokenize(redo);
        let name = blocks
            .iter()
            .find_map(|t| {
                if let Token::DataBlock(s) = t {
                    Some(*s)
                } else {
                    None
                }
            })
            .unwrap();
        if canonical_id(name) != expected {
            let start = name.as_ptr() as usize - redo.as_ptr() as usize;
            edits.push((start..start + name.len(), expected.clone()));
        }
        match target.get("_entry") {
            Some(entry) => {
                entry_id(&target)?;
                let cell = &entry.rows[0][entry.col("_entry.id")?];
                if canonical_id(&cell.value) != expected {
                    edits.push((cell.span.clone(), expected.clone()));
                }
            }
            None => edits.push((redo.len()..redo.len(), format!("\n_entry.id {expected}\n"))),
        }
    }
    let mut additions = Vec::new();
    if !selected.is_empty() {
        let old_seq = sequence(table(&target, "_entity_poly_seq")?)?;
        let ref_seq = sequence(table(&source, "_entity_poly_seq")?)?;
        for entity in &selected {
            let mapped: BTreeSet<_> = mapping
                .iter()
                .filter(|m| &m.1 == entity)
                .map(|m| m.3.as_str())
                .collect();
            if mapped.len() != 1 {
                return Err(invalid("entity maps to multiple reference entities"));
            }
            let source_entity = mapped.into_iter().next().unwrap();
            let declared: BTreeMap<_, _> = old_seq
                .iter()
                .filter(|((e, _), _)| e == entity)
                .map(|((_, n), m)| (*n, m))
                .collect();
            let reference_seq: BTreeMap<_, _> = ref_seq
                .iter()
                .filter(|((e, _), _)| e == source_entity)
                .map(|((_, n), m)| (*n, m))
                .collect();
            if declared.is_empty() || declared != reference_seq {
                return Err(invalid(format!(
                    "full sequence disagreement for target entity {entity}"
                )));
            }
            let atoms = table(&target, "_atom_site")?;
            for row in &atoms.rows {
                if atom_entity(&target, atoms, row)? != entity {
                    continue;
                }
                let n = atoms
                    .value(row, "_atom_site.label_seq_id")?
                    .parse::<u32>()
                    .map_err(|_| invalid("missing observed label sequence ID"))?;
                if declared.get(&n).map(|m| m.as_str())
                    != Some(atoms.value(row, "_atom_site.label_comp_id")?)
                {
                    return Err(invalid(format!(
                        "observed sequence conflict for entity {entity}, position {n}"
                    )));
                }
            }
            let kind = polymer_type(&source, source_entity)?
                .ok_or_else(|| invalid("reference polymer type missing"))?;
            if missing(&kind.value) || kind.value == "other" {
                return Err(invalid("reference polymer type is not informative"));
            }
            match polymer_type(&target, entity)? {
                Some(old) if missing(&old.value) || old.value == "other" => {
                    edits.push((old.span.clone(), reference[kind.span.clone()].to_owned()))
                }
                Some(_) => return Err(invalid("selected polymer type is already specified")),
                None => {
                    if !entity
                        .bytes()
                        .all(|b| b.is_ascii_alphanumeric() || b == b'_')
                    {
                        return Err(invalid("unsafe entity token for metadata insertion"));
                    }
                    additions.push((entity.clone(), reference[kind.span.clone()].to_owned()));
                }
            }
        }
    }
    if !additions.is_empty() {
        if let Some(poly) = target.get("_entity_poly") {
            // Appending to a loop preserves all its existing lexical metadata.
            // Scalar-to-loop conversion is deliberately unresolved in v1.
            if !poly.is_loop {
                return Err(invalid(
                    "adding polymer rows to scalar metadata is not supported",
                ));
            }
            let entity_col = poly.col("_entity_poly.entity_id")?;
            let type_col = poly.col("_entity_poly.type")?;
            let mut rows = String::new();
            for (entity, kind) in additions {
                let mut values = vec!["?".to_owned(); poly.columns.len()];
                values[entity_col] = entity;
                values[type_col] = kind;
                rows.push_str(&format!("\n{}", values.join(" ")));
            }
            edits.push((poly.span.end..poly.span.end, rows));
        } else {
            let mut text = "\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n".to_owned();
            for (entity, kind) in additions {
                text.push_str(&format!("{entity} {kind}\n"));
            }
            edits.push((redo.len()..redo.len(), text));
        }
    }
    edits.sort_by_key(|(span, _)| (span.start, span.end));
    for pair in edits.windows(2) {
        if pair[0].0.end > pair[1].0.start {
            return Err(invalid("overlapping metadata edits"));
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
    let output = tables(&cif)?;
    let before = table(&target, "_atom_site")?;
    let after = table(&output, "_atom_site")?;
    if before.columns != after.columns || before.rows.len() != after.rows.len() {
        return Err(invalid("atom-site layout changed"));
    }
    for (old, new) in before.rows.iter().zip(&after.rows) {
        for (a, b) in old.iter().zip(new) {
            if redo[a.span.clone()] != cif[b.span.clone()] {
                return Err(invalid("protected atom token changed"));
            }
        }
    }
    if canonical_id(entry_id(&output)?) != expected {
        return Err(invalid("output entry identity mismatch"));
    }
    super::super::read_cif_str(&cif)?;
    Ok(IdentityTypeRepair {
        cif,
        chain_mapping: mapping,
        changes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn doc(id: &str, entity: &str, kind: &str) -> String {
        format!("data_{id}\n_entry.id {id}\n_entity.id {entity}\n_entity.type polymer\n_entity_poly.entity_id {entity}\n_entity_poly.type {kind}\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n{entity} 1 GLY\n#\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 C CA GLY X {entity} 1 A 10 ? 1.234500 2.000 3.0\n")
    }

    #[test]
    fn independent_entities_and_identity_are_explicit() {
        let target = doc("100", "5", "?");
        let reference = doc("1e29", "3", "'polypeptide(L)'");
        let result =
            repair_identity_and_polymer_types(&target, &reference, "1e29", true, &["5".into()])
                .unwrap();
        assert!(result.cif.starts_with("data_1e29\n_entry.id 1e29"));
        assert_eq!(
            result.chain_mapping,
            vec![("X".into(), "5".into(), "X".into(), "3".into())]
        );
        assert_eq!(result.changes.len(), 3);
        assert!(repair_identity_and_polymer_types(
            &target,
            &reference,
            "1e29",
            false,
            &["5".into()]
        )
        .is_err());
        assert!(repair_identity_and_polymer_types(&target, &reference, "2abc", true, &[]).is_err());
    }

    #[test]
    fn fail_closed_on_unanchored_and_conflicting_sequence() {
        let target = doc("1e29", "5", "other");
        let reference = doc("1e29", "3", "'polypeptide(L)'");
        assert!(repair_identity_and_polymer_types(
            &target.replace("A 10 ?", "A 11 ?"),
            &reference,
            "1e29",
            false,
            &["5".into()]
        )
        .is_err());
        assert!(repair_identity_and_polymer_types(
            &target.replace("5 1 GLY", "5 1 ALA"),
            &reference,
            "1e29",
            false,
            &["5".into()]
        )
        .is_err());
        let ambiguous = reference.replace("1 C CA GLY X 3 1 A 10 ? 1.234500 2.000 3.0", "1 C CA GLY X 3 1 A 10 ? 1.234500 2.000 3.0\n2 C CA GLY Y 3 1 A 10 ? 1.234500 2.000 3.0");
        assert!(repair_identity_and_polymer_types(
            &target,
            &ambiguous,
            "1e29",
            false,
            &["5".into()]
        )
        .is_err());
    }

    #[test]
    fn missing_polymer_table_and_loop_rows_are_explicit_additions() {
        let reference = doc("1e29", "3", "'polypeptide(L)'");
        let target =
            doc("1e29", "5", "?").replace("_entity_poly.entity_id 5\n_entity_poly.type ?\n", "");
        let result =
            repair_identity_and_polymer_types(&target, &reference, "1e29", false, &["5".into()])
                .unwrap();
        assert_eq!(result.changes.len(), 1);
        assert!(result.cif.contains("5 'polypeptide(L)'"));
        let loop_target = target
            .replace(
                "_entity.id 5",
                "loop_\n_entity.id\n_entity.type\n6 polymer\n5",
            )
            .replace("_entity.type polymer", "polymer");
        let loop_target = format!("{loop_target}\nloop_\n_entity_poly.entity_id\n_entity_poly.type\n6 'polypeptide(L)'\n#\n");
        let result = repair_identity_and_polymer_types(
            &loop_target,
            &reference,
            "1e29",
            false,
            &["5".into()],
        )
        .unwrap();
        assert!(result
            .cif
            .contains("6 'polypeptide(L)'\n5 'polypeptide(L)'"));
    }

    #[test]
    fn journal_reconstructs_output_and_coordinates_do_not_anchor() {
        let target = doc("100", "5", "other");
        let reference =
            doc("1e29", "3", "'polypeptide(L)'").replace("1.234500 2.000 3.0", "90.0 80.0 70.0");
        let result =
            repair_identity_and_polymer_types(&target, &reference, "1e29", true, &["5".into()])
                .unwrap();
        let mut replay = target.clone();
        for (start, end, old, new) in result.changes.iter().rev() {
            assert_eq!(&target[*start..*end], old);
            replay.replace_range(*start..*end, new);
        }
        assert_eq!(replay, result.cif);
        assert!(result.cif.contains("1.234500 2.000 3.0"));
        assert!(!result.cif.contains("90.0"));
        assert!(repair_identity_and_polymer_types(
            &target.replace("A 10 ?", "A 10 B"),
            &reference,
            "1e29",
            true,
            &[]
        )
        .is_err());
    }

    #[test]
    fn every_residue_and_selected_entity_must_be_anchored() {
        let reference = doc("1e29", "3", "'polypeptide(L)'");
        let target = doc("1e29", "5", "other");
        let extra = format!("{target}2 C CA GLY X 5 1 A 11 ? 1.0 2.0 3.0\n");
        assert!(repair_identity_and_polymer_types(
            &extra,
            &reference,
            "1e29",
            false,
            &["5".into()]
        )
        .is_err());
        let absent = target.replace(
            "_entity.id 5\n_entity.type polymer",
            "loop_\n_entity.id\n_entity.type\n5 polymer\n6 polymer\n#",
        );
        assert!(repair_identity_and_polymer_types(&absent, &reference, "1e29", true, &[]).is_err());
        assert!(
            repair_identity_and_polymer_types(&target, &reference, "1e29", false, &[]).is_err()
        );
        assert!(repair_identity_and_polymer_types(
            &target.replace("other", "'polypeptide(L)'"),
            &reference,
            "1e29",
            false,
            &["5".into()]
        )
        .is_err());
    }

    #[test]
    fn missing_atom_entity_uses_proven_struct_asym_without_token_edits() {
        let target = doc("1e29", "5", "other").replace("X 5 1 A", "X ? 1 A");
        let target = format!("{target}\n_struct_asym.id X\n_struct_asym.entity_id 5\n");
        let reference = doc("1e29", "3", "'polypeptide(L)'");
        let result =
            repair_identity_and_polymer_types(&target, &reference, "1e29", false, &["5".into()])
                .unwrap();
        assert!(result.cif.contains("X ? 1 A"));
        let conflict = target.replace("X ? 1 A", "X 6 1 A");
        assert!(repair_identity_and_polymer_types(
            &conflict,
            &reference,
            "1e29",
            false,
            &["5".into()]
        )
        .is_err());
    }
}
