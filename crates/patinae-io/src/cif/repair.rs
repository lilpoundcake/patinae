//! Explicit, reference-backed metadata repair. Normal readers remain strict.
//!
//! Only known missing water labels and missing terminal sequence rows are
//! restored. Coordinates and other atom-site tokens are never reserialized.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::ops::Range;

use crate::{IoError, IoResult};

use super::lexer::{tokenize, Token};

mod identity;
pub use identity::{repair_identity_and_polymer_types, IdentityTypeRepair};
mod terminal;
pub use terminal::{repair_anchored_terminal_sequences, TerminalSequenceRepair};

#[derive(Debug)]
pub struct MetadataRepair {
    pub cif: String,
    pub water_atoms: usize,
    pub sequence_rows_added: usize,
    /// Byte offsets in the original decompressed CIF, old text, new text.
    pub changes: Vec<(usize, usize, String, String)>,
}

#[derive(Clone)]
struct Cell {
    value: String,
    span: Range<usize>,
}

struct Table {
    columns: Vec<String>,
    rows: Vec<Vec<Cell>>,
    span: Range<usize>,
    is_loop: bool,
}

fn invalid(message: impl Into<String>) -> IoError {
    IoError::parse_msg(message.into())
}

impl Table {
    fn col(&self, field: &str) -> IoResult<usize> {
        self.columns
            .iter()
            .position(|name| name == field)
            .ok_or_else(|| invalid(format!("repair requires column {field}")))
    }

    fn value<'a>(&self, row: &'a [Cell], field: &str) -> IoResult<&'a str> {
        Ok(&row[self.col(field)?].value)
    }
}

fn missing(value: &str) -> bool {
    matches!(value, "." | "?" | "")
}

/// Locate original token spans using the same lexer as the strict reader.
fn tables(text: &str) -> IoResult<HashMap<String, Table>> {
    let tokens = tokenize(text);
    let bytes = text.as_bytes();
    let mut cursor = 0;
    let mut spans = Vec::with_capacity(tokens.len());
    let mut blocks = 0;
    for token in &tokens {
        loop {
            while cursor < bytes.len() && bytes[cursor].is_ascii_whitespace() {
                cursor += 1;
            }
            if bytes.get(cursor) != Some(&b'#') {
                break;
            }
            while cursor < bytes.len() && bytes[cursor] != b'\n' {
                cursor += 1;
            }
        }
        let start = cursor;
        cursor = match token {
            Token::Loop => cursor + 5,
            Token::Missing | Token::Unknown => cursor + 1,
            Token::Eof => bytes.len(),
            Token::SingleQuoted(s) | Token::DoubleQuoted(s) => {
                s.as_ptr() as usize - text.as_ptr() as usize + s.len() + 1
            }
            Token::TextField(s) => s.as_ptr() as usize - text.as_ptr() as usize + s.len() + 2,
            _ => {
                let s = token
                    .as_str()
                    .ok_or_else(|| invalid("unsupported CIF token"))?;
                s.as_ptr() as usize - text.as_ptr() as usize + s.len()
            }
        };
        if cursor > bytes.len() {
            return Err(invalid("unterminated CIF value"));
        }
        if matches!(token, Token::DataBlock(_)) {
            blocks += 1;
        }
        spans.push(start..cursor);
    }
    if blocks != 1 {
        return Err(invalid("metadata repair requires exactly one data block"));
    }
    let relevant = |category: &str| {
        matches!(
            category,
            "_atom_site"
                | "_entry"
                | "_entity_poly"
                | "_entity"
                | "_entity_poly_seq"
                | "_struct_asym"
                | "_pdbx_struct_assembly_gen"
                | "_pdbx_struct_oper_list"
        )
    };
    let cell = |index: usize| -> IoResult<Cell> {
        let value = match &tokens[index] {
            Token::Missing => ".",
            Token::Unknown => "?",
            Token::Value(s)
            | Token::SingleQuoted(s)
            | Token::DoubleQuoted(s)
            | Token::TextField(s) => s,
            _ => return Err(invalid("expected CIF value")),
        };
        Ok(Cell {
            value: value.to_string(),
            span: spans[index].clone(),
        })
    };
    let mut result: HashMap<String, Table> = HashMap::new();
    let mut i = 0;
    while i < tokens.len() {
        if matches!(tokens[i], Token::Loop) {
            let start = i;
            i += 1;
            let mut columns = Vec::new();
            while let Some(Token::DataName(name)) = tokens.get(i) {
                columns.push(name.to_string());
                i += 1;
            }
            let values = i;
            while i < tokens.len()
                && !matches!(
                    tokens[i],
                    Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof
                )
            {
                i += 1;
            }
            let category = columns.first().and_then(|s| s.split_once('.')).map(|s| s.0);
            if let Some(category) = category.filter(|s| relevant(s)) {
                if columns
                    .iter()
                    .any(|s| s.split_once('.').map(|s| s.0) != Some(category))
                    || (i - values) % columns.len() != 0
                {
                    return Err(invalid("incomplete/mixed repair table"));
                }
                let rows = (values..i)
                    .step_by(columns.len())
                    .map(|offset| (offset..offset + columns.len()).map(&cell).collect())
                    .collect::<IoResult<Vec<Vec<Cell>>>>()?;
                let end = if i > values {
                    spans[i - 1].end
                } else {
                    spans[values - 1].end
                };
                let key = category.to_owned();
                if result
                    .insert(
                        key,
                        Table {
                            columns,
                            rows,
                            span: spans[start].start..end,
                            is_loop: true,
                        },
                    )
                    .is_some()
                {
                    return Err(invalid("duplicate repair category"));
                }
            }
        } else if let Token::DataName(name) = &tokens[i] {
            if let Some((category, _)) = name.split_once('.') {
                if relevant(category) {
                    let value = cell(i + 1)?;
                    let table = result.entry(category.to_owned()).or_insert_with(|| Table {
                        columns: Vec::new(),
                        rows: vec![Vec::new()],
                        span: spans[i].start..value.span.end,
                        is_loop: false,
                    });
                    if table.is_loop || table.columns.iter().any(|s| s == name) {
                        return Err(invalid("duplicate single-item repair category"));
                    }
                    table.span.end = value.span.end;
                    table.columns.push(name.to_string());
                    table.rows[0].push(value);
                }
            }
            i += 2;
        } else {
            i += 1;
        }
    }
    Ok(result)
}

fn table<'a>(tables: &'a HashMap<String, Table>, name: &str) -> IoResult<&'a Table> {
    tables
        .get(name)
        .ok_or_else(|| invalid(format!("repair requires {name}")))
}

fn atom_key(table: &Table, row: &[Cell]) -> IoResult<Vec<String>> {
    [
        "auth_asym_id",
        "auth_seq_id",
        "label_comp_id",
        "label_atom_id",
        "pdbx_PDB_ins_code",
    ]
    .iter()
    .map(|field| {
        let value = table.value(row, &format!("_atom_site.{field}"))?;
        Ok(if missing(value) {
            String::new()
        } else {
            value.to_owned()
        })
    })
    .collect()
}

fn sequence(table: &Table) -> IoResult<BTreeMap<(String, u32), String>> {
    let mut result = BTreeMap::new();
    for row in &table.rows {
        let entity = table.value(row, "_entity_poly_seq.entity_id")?.to_owned();
        let num = table
            .value(row, "_entity_poly_seq.num")?
            .parse::<u32>()
            .map_err(|_| invalid("invalid reference sequence position"))?;
        let mon = table.value(row, "_entity_poly_seq.mon_id")?.to_owned();
        if result.insert((entity, num), mon).is_some() {
            return Err(invalid("ambiguous reference sequence position"));
        }
    }
    Ok(result)
}

/// Restore only explicitly requested metadata, using an independently pinned
/// reference CIF. Fails closed on ambiguous matches, sequence disagreements,
/// non-terminal holes, non-identity operations or changed polymer membership.
pub fn repair_metadata(
    redo: &str,
    reference: &str,
    water: bool,
    terminal_sequence: bool,
) -> IoResult<MetadataRepair> {
    repair_metadata_with_assemblies(redo, reference, water, terminal_sequence, false)
}

/// Explicitly opt in to rebuilding assembly lists from reference-backed chain
/// identities. Every target chain must have unambiguous reference membership;
/// unmatched atoms may inherit only an otherwise anchored chain's membership.
pub fn repair_metadata_with_assemblies(
    redo: &str,
    reference: &str,
    water: bool,
    terminal_sequence: bool,
    assembly_membership: bool,
) -> IoResult<MetadataRepair> {
    if assembly_membership && !water {
        return Err(invalid("assembly membership repair requires water repair"));
    }
    let block_id = |text: &str| -> IoResult<String> {
        tokenize(text)
            .into_iter()
            .find_map(|t| match t {
                Token::DataBlock(name) => Some(name.to_ascii_lowercase()),
                _ => None,
            })
            .ok_or_else(|| invalid("missing reference block ID"))
    };
    if block_id(redo)? != block_id(reference)? {
        return Err(invalid("reference block ID differs from target"));
    }
    let target = tables(redo)?;
    let source = tables(reference)?;
    let mut edits: Vec<(Range<usize>, String)> = Vec::new();
    let mut sequence_rows_added = 0;
    if terminal_sequence {
        let old_table = table(&target, "_entity_poly_seq")?;
        let new_table = table(&source, "_entity_poly_seq")?;
        let old = sequence(old_table)?;
        let full = sequence(new_table)?;
        if old.is_empty()
            || old.keys().map(|k| &k.0).collect::<BTreeSet<_>>()
                != full.keys().map(|k| &k.0).collect::<BTreeSet<_>>()
        {
            return Err(invalid("reference sequence entities do not match"));
        }
        for (key, mon) in &old {
            if full.get(key) != Some(mon) {
                return Err(invalid("reference sequence disagreement"));
            }
        }
        for entity in old.keys().map(|k| &k.0).collect::<BTreeSet<_>>() {
            let positions: Vec<_> = old.keys().filter(|k| &k.0 == entity).map(|k| k.1).collect();
            if positions.windows(2).any(|w| w[1] != w[0] + 1) {
                return Err(invalid(
                    "only missing terminal sequence positions can be restored",
                ));
            }
            let positions: Vec<_> = full
                .keys()
                .filter(|k| &k.0 == entity)
                .map(|k| k.1)
                .collect();
            if positions.first() != Some(&1) || positions.windows(2).any(|w| w[1] != w[0] + 1) {
                return Err(invalid("reference sequence is not complete"));
            }
        }
        sequence_rows_added = full.len() - old.len();
        edits.push((
            old_table.span.clone(),
            reference[new_table.span.clone()].to_owned(),
        ));
    }
    let mut water_atoms = 0;
    if water {
        let atoms = table(&target, "_atom_site")?;
        let ref_atoms = table(&source, "_atom_site")?;
        let entities = table(&target, "_entity")?;
        let waters = entities
            .rows
            .iter()
            .filter_map(|row| {
                (entities.value(row, "_entity.type").ok()? == "water")
                    .then(|| entities.value(row, "_entity.id").ok().map(str::to_owned))
                    .flatten()
            })
            .collect::<Vec<_>>();
        if waters.len() != 1 {
            return Err(invalid("requires one unambiguous water entity"));
        }
        let water_entity = &waters[0];
        let mut index: BTreeMap<Vec<String>, BTreeSet<(String, String)>> = BTreeMap::new();
        for row in &ref_atoms.rows {
            index.entry(atom_key(ref_atoms, row)?).or_default().insert((
                ref_atoms.value(row, "_atom_site.label_asym_id")?.to_owned(),
                ref_atoms
                    .value(row, "_atom_site.label_entity_id")?
                    .to_owned(),
            ));
        }
        let asym = table(&target, "_struct_asym")?;
        let mut used: BTreeSet<String> = asym
            .rows
            .iter()
            .map(|r| asym.value(r, "_struct_asym.id").map(str::to_owned))
            .collect::<IoResult<_>>()?;
        for group in &table(&target, "_pdbx_struct_assembly_gen")?.rows {
            used.extend(
                table(&target, "_pdbx_struct_assembly_gen")?
                    .value(group, "_pdbx_struct_assembly_gen.asym_id_list")?
                    .split(',')
                    .map(str::to_owned),
            );
        }
        let mut author_water: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        for row in &atoms.rows {
            let label = atoms.value(row, "_atom_site.label_asym_id")?;
            used.insert(label.to_owned());
            if atoms.value(row, "_atom_site.label_entity_id")? == water_entity && !missing(label) {
                author_water
                    .entry(atoms.value(row, "_atom_site.auth_asym_id")?.to_owned())
                    .or_default()
                    .insert(label.to_owned());
            }
        }
        let mut added = BTreeMap::new();
        let mut source_to_target = BTreeMap::new();
        for row in &atoms.rows {
            let label = atoms.value(row, "_atom_site.label_asym_id")?;
            let entity = atoms.value(row, "_atom_site.label_entity_id")?;
            // Other incomplete metadata is outside this explicitly scoped
            // repair. Only anonymous water caused the assembly failure.
            if !missing(label) {
                continue;
            }
            if !missing(label)
                || !missing(entity)
                || atoms.value(row, "_atom_site.label_comp_id")? != "HOH"
            {
                return Err(invalid(
                    "repair supports only water with both labels missing",
                ));
            }
            let matches = index
                .get(&atom_key(atoms, row)?)
                .ok_or_else(|| invalid("water absent from reference"))?;
            if matches.len() != 1 {
                return Err(invalid("ambiguous reference water atom"));
            }
            let (source_label, source_entity) = matches.iter().next().unwrap();
            let ref_entities = table(&source, "_entity")?;
            if !ref_entities.rows.iter().any(|r| {
                ref_entities.value(r, "_entity.id").ok() == Some(source_entity)
                    && ref_entities.value(r, "_entity.type").ok() == Some("water")
            }) || missing(source_label)
            {
                return Err(invalid("reference match is not labelled water"));
            }
            let author = atoms.value(row, "_atom_site.auth_asym_id")?;
            let labels = author_water.entry(author.to_owned()).or_default();
            if labels.len() > 1 {
                return Err(invalid("ambiguous target water chain"));
            }
            let target_label = if let Some(label) = labels.iter().next() {
                label.clone()
            } else {
                let label = (0..1000)
                    .map(|i| format!("R{i}"))
                    .find(|s| !used.contains(s))
                    .ok_or_else(|| invalid("no free repair chain identifier"))?;
                used.insert(label.clone());
                labels.insert(label.clone());
                added.insert(label.clone(), water_entity.clone());
                label
            };
            if let Some(previous) =
                source_to_target.insert(source_label.clone(), target_label.clone())
            {
                if previous != target_label {
                    return Err(invalid("conflicting water mapping"));
                }
            }
            for (field, value) in [
                ("label_asym_id", target_label),
                ("label_entity_id", water_entity.clone()),
            ] {
                edits.push((
                    row[atoms.col(&format!("_atom_site.{field}"))?].span.clone(),
                    value,
                ));
            }
            water_atoms += 1;
        }
        for (label, entity) in &added {
            if !asym.is_loop {
                return Err(invalid("repair requires loop struct_asym"));
            }
            let mut values = vec!["?".to_owned(); asym.columns.len()];
            values[asym.col("_struct_asym.id")?] = label.clone();
            values[asym.col("_struct_asym.entity_id")?] = entity.clone();
            // Combine all insertions at this offset below.
            edits.push((
                asym.span.end..asym.span.end,
                format!("\n{}", values.join(" ")),
            ));
        }
        let groups = table(&target, "_pdbx_struct_assembly_gen")?;
        let ref_groups = table(&source, "_pdbx_struct_assembly_gen")?;
        let mut chain_sources: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
        if assembly_membership {
            let signatures = |t: &Table| -> IoResult<BTreeSet<(String, String)>> {
                let keys: Vec<_> = t
                    .rows
                    .iter()
                    .map(|r| {
                        Ok((
                            t.value(r, "_pdbx_struct_assembly_gen.assembly_id")?
                                .to_owned(),
                            t.value(r, "_pdbx_struct_assembly_gen.oper_expression")?
                                .to_owned(),
                        ))
                    })
                    .collect::<IoResult<_>>()?;
                let unique: BTreeSet<_> = keys.iter().cloned().collect();
                if unique.len() != keys.len() {
                    return Err(invalid("duplicate assembly group signature"));
                }
                Ok(unique)
            };
            if signatures(groups)? != signatures(ref_groups)? {
                return Err(invalid("reference assembly signatures differ"));
            }
            for row in &atoms.rows {
                let label = atoms.value(row, "_atom_site.label_asym_id")?;
                let matches = index.get(&atom_key(atoms, row)?);
                let label = if missing(label) {
                    let matches = matches.ok_or_else(|| invalid("unanchored water chain"))?;
                    let source_label = &matches.iter().next().unwrap().0;
                    source_to_target
                        .get(source_label)
                        .ok_or_else(|| invalid("missing repaired water chain"))?
                        .as_str()
                } else {
                    label
                };
                let sources = chain_sources.entry(label.to_owned()).or_default();
                if let Some(matches) = matches {
                    if matches.len() != 1 {
                        return Err(invalid("ambiguous reference atom identity"));
                    }
                    let source_label = &matches.iter().next().unwrap().0;
                    if missing(source_label) {
                        return Err(invalid("unlabelled reference chain"));
                    }
                    sources.insert(source_label.clone());
                }
            }
            if chain_sources.values().any(BTreeSet::is_empty) {
                return Err(invalid("target chain has no reference anchor"));
            }
        }
        for group in &groups.rows {
            let list_col = groups.col("_pdbx_struct_assembly_gen.asym_id_list")?;
            let list = &group[list_col].value;
            if !assembly_membership && !list.split(',').any(|s| s.trim().is_empty()) {
                continue;
            }
            let id = groups.value(group, "_pdbx_struct_assembly_gen.assembly_id")?;
            let expression = groups.value(group, "_pdbx_struct_assembly_gen.oper_expression")?;
            if expression != "1" {
                return Err(invalid("water repair requires identity assembly operation"));
            }
            for document in [&target, &source] {
                let operations = table(document, "_pdbx_struct_oper_list")?;
                let op = operations
                    .rows
                    .iter()
                    .find(|r| operations.value(r, "_pdbx_struct_oper_list.id").ok() == Some("1"))
                    .ok_or_else(|| invalid("missing identity operator"))?;
                for i in 1..=3 {
                    for j in 1..=3 {
                        let value: f64 = operations
                            .value(op, &format!("_pdbx_struct_oper_list.matrix[{i}][{j}]"))?
                            .parse()
                            .map_err(|_| invalid("invalid operator"))?;
                        if value != if i == j { 1.0 } else { 0.0 } {
                            return Err(invalid("non-identity matrix"));
                        }
                    }
                    let value: f64 = operations
                        .value(op, &format!("_pdbx_struct_oper_list.vector[{i}]"))?
                        .parse()
                        .map_err(|_| invalid("invalid operator"))?;
                    if value != 0.0 {
                        return Err(invalid("nonzero assembly translation"));
                    }
                }
            }
            let matches: Vec<_> = ref_groups
                .rows
                .iter()
                .filter(|r| {
                    ref_groups
                        .value(r, "_pdbx_struct_assembly_gen.assembly_id")
                        .ok()
                        == Some(id)
                        && ref_groups
                            .value(r, "_pdbx_struct_assembly_gen.oper_expression")
                            .ok()
                            == Some(expression)
                })
                .collect();
            if matches.len() != 1 {
                return Err(invalid("ambiguous reference assembly group"));
            }
            let ref_list: BTreeSet<_> = ref_groups
                .value(matches[0], "_pdbx_struct_assembly_gen.asym_id_list")?
                .split(',')
                .collect();
            let mut labels: Vec<String> = list
                .split(',')
                .filter(|s| !s.trim().is_empty())
                .map(str::to_owned)
                .collect();
            if assembly_membership {
                let mut expected = BTreeSet::new();
                for (label, sources) in &chain_sources {
                    let membership: BTreeSet<_> = sources
                        .iter()
                        .map(|s| ref_list.contains(s.as_str()))
                        .collect();
                    if membership.len() != 1 {
                        return Err(invalid(
                            "target chain spans conflicting reference assemblies",
                        ));
                    }
                    if membership.contains(&true) {
                        expected.insert(label.clone());
                    }
                }
                labels.retain(|label| expected.contains(label));
                for label in expected {
                    if !labels.contains(&label) {
                        labels.push(label);
                    }
                }
                if labels.is_empty() {
                    return Err(invalid("repaired assembly has no members"));
                }
            }
            let old_labels: BTreeSet<_> = labels.iter().map(String::as_str).collect();
            let authors = |t: &Table, labels: &BTreeSet<&str>| -> IoResult<BTreeSet<String>> {
                let mut authors = BTreeSet::new();
                for atom in &t.rows {
                    if labels.contains(t.value(atom, "_atom_site.label_asym_id")?)
                        && !missing(t.value(atom, "_atom_site.label_seq_id")?)
                    {
                        authors.insert(t.value(atom, "_atom_site.auth_asym_id")?.to_owned());
                    }
                }
                Ok(authors)
            };
            if authors(atoms, &old_labels)? != authors(ref_atoms, &ref_list)? {
                return Err(invalid(
                    "polymer assembly membership differs from reference",
                ));
            }
            let mut mapped = false;
            for (source_label, target_label) in &source_to_target {
                if ref_list.contains(source_label.as_str()) {
                    mapped = true;
                    if !labels.contains(target_label) {
                        labels.push(target_label.clone());
                    }
                }
            }
            if !mapped {
                return Err(invalid(
                    "empty assembly member has no reference water mapping",
                ));
            }
            edits.push((group[list_col].span.clone(), labels.join(",")));
        }
    }
    // Reverse-position edits preserve original lexical bytes outside metadata.
    // Merge insertions at a common offset to avoid ambiguous edit ordering.
    edits.sort_by_key(|edit| (edit.0.start, edit.0.end));
    for pair in edits.windows(2) {
        if pair[0].0.end > pair[1].0.start {
            return Err(invalid("overlapping metadata edits"));
        }
    }
    let changes = edits
        .iter()
        .map(|(span, value)| {
            (
                span.start,
                span.end,
                redo[span.clone()].to_owned(),
                value.clone(),
            )
        })
        .collect();
    let mut cif = redo.to_owned();
    for (span, value) in edits.into_iter().rev() {
        cif.replace_range(span, &value);
    }
    let repaired_tables = tables(&cif)?;
    let before = table(&target, "_atom_site")?;
    let after = table(&repaired_tables, "_atom_site")?;
    if before.columns != after.columns || before.rows.len() != after.rows.len() {
        return Err(invalid("repair changed atom-site layout"));
    }
    for (old, new) in before.rows.iter().zip(&after.rows) {
        for (i, field) in before.columns.iter().enumerate() {
            let permitted = water
                && missing(&old[i].value)
                && matches!(
                    field.as_str(),
                    "_atom_site.label_asym_id" | "_atom_site.label_entity_id"
                );
            // Compare the original lexical bytes too: no rounding, reordering,
            // coordinate replacement, or silent changes to auth/label_seq IDs.
            if !permitted && redo[old[i].span.clone()] != cif[new[i].span.clone()] {
                return Err(invalid(format!(
                    "repair changed protected atom field {field}"
                )));
            }
        }
    }
    // The output must satisfy the normal strict reader, without repair flags.
    super::read_cif_str(&cif)?;
    Ok(MetadataRepair {
        cif,
        water_atoms,
        sequence_rows_added,
        changes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const ATOMS: &str = "loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 C CA GLY A 1 2 A 10 ? 1.234500 2.000 3.0\n";

    fn document(sequence: &str) -> String {
        format!("data_test\nloop_\n_entity.id\n_entity.type\n1 polymer\n#\n_entity_poly.entity_id 1\n_entity_poly.type 'polypeptide(L)'\n#\nloop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n{sequence}\n#\n{ATOMS}")
    }

    #[test]
    fn restores_terminal_sequence_without_changing_atom_bytes() {
        let old = document("1 2 GLY");
        let full = document("1 1 MET\n1 2 GLY\n1 3 ALA");
        let repaired = repair_metadata(&old, &full, false, true).unwrap();
        assert_eq!(repaired.sequence_rows_added, 2);
        assert_eq!(repaired.water_atoms, 0);
        assert!(repaired.cif.contains(ATOMS));
        assert_eq!(
            repair_metadata(&repaired.cif, &full, false, true)
                .unwrap()
                .cif,
            repaired.cif
        );
    }

    #[test]
    fn rejects_sequence_disagreement_internal_holes_and_wrong_entry() {
        let full = document("1 1 MET\n1 2 GLY\n1 3 ALA");
        assert!(repair_metadata(&document("1 2 VAL"), &full, false, true).is_err());
        assert!(repair_metadata(&document("1 1 MET\n1 3 ALA"), &full, false, true).is_err());
        assert!(repair_metadata(
            &document("1 2 GLY"),
            &full.replace("data_test", "data_other"),
            false,
            true
        )
        .is_err());
    }

    fn water_documents() -> (String, String) {
        let header = "data_test\nloop_\n_entity.id\n_entity.type\n1 polymer\n2 water\n#\nloop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\n#\nloop_\n_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.matrix[1][1]\n_pdbx_struct_oper_list.matrix[1][2]\n_pdbx_struct_oper_list.matrix[1][3]\n_pdbx_struct_oper_list.vector[1]\n_pdbx_struct_oper_list.matrix[2][1]\n_pdbx_struct_oper_list.matrix[2][2]\n_pdbx_struct_oper_list.matrix[2][3]\n_pdbx_struct_oper_list.vector[2]\n_pdbx_struct_oper_list.matrix[3][1]\n_pdbx_struct_oper_list.matrix[3][2]\n_pdbx_struct_oper_list.matrix[3][3]\n_pdbx_struct_oper_list.vector[3]\n1 1 0 0 0 0 1 0 0 0 0 1 0\n#\nloop_\n_pdbx_struct_assembly.id\n1\n#\nloop_\n_pdbx_struct_assembly_gen.assembly_id\n_pdbx_struct_assembly_gen.oper_expression\n_pdbx_struct_assembly_gen.asym_id_list\n";
        let old =
            format!("{header}1 1 ,A\n#\n{ATOMS}2 O O HOH ? ? . A 501 ? 4.0000 5.0000 6.0000\n");
        let source =
            format!("{header}1 1 A,W\n#\n{ATOMS}2 O O HOH W 2 . A 501 ? 9.0000 9.0000 9.0000\n");
        (old, source)
    }

    #[test]
    fn repairs_water_from_identity_not_reference_coordinates() {
        let (old, source) = water_documents();
        let repaired = repair_metadata(&old, &source, true, false).unwrap();
        assert_eq!(repaired.water_atoms, 1);
        assert!(repaired
            .cif
            .contains("2 O O HOH R0 2 . A 501 ? 4.0000 5.0000 6.0000"));
        assert!(repaired.cif.contains("A,R0"));
        assert!(repaired.cif.contains("R0 2"));
    }

    #[test]
    fn rejects_unknown_or_ambiguous_water_and_changed_membership() {
        let (old, source) = water_documents();
        assert!(repair_metadata(&old, &source.replace("A 501", "A 999"), true, false).is_err());
        let ambiguous = format!("{source}3 O O HOH Z 2 . A 501 ? 9 9 9\n");
        assert!(repair_metadata(&old, &ambiguous, true, false).is_err());
        assert!(repair_metadata(&old, &source.replace("1 1 A,W", "1 1 W"), true, false).is_err());
    }

    #[test]
    fn explicit_membership_mode_separates_protein_ligands_and_water() {
        let (old, source) = water_documents();
        let old = old.replace("1 1 ,A", "1 1 ,A,U\n2 1 ,A,B,U,V")
            .replace("2 O O HOH ? ? . A 501 ? 4.0000 5.0000 6.0000", "2 O O HOH ? ? . A 501 ? 4 5 6\n3 O O HOH ? ? . B 502 ? 7 8 9\n4 C CA GLY B 1 2 B 10 ? 10 11 12\n5 C C1 EDO U 3 . A 101 ? 13 14 15\n6 C C1 EDO V 3 . B 102 ? 16 17 18");
        let source = source.replace("1 1 A,W", "1 1 A,X,WA\n2 1 B,Y,WB")
            .replace("2 O O HOH W 2 . A 501 ? 9.0000 9.0000 9.0000", "2 O O HOH WA 2 . A 501 ? 90 90 90\n3 O O HOH WB 2 . B 502 ? 90 90 90\n4 C CA GLY B 1 2 B 10 ? 90 90 90\n5 C C1 EDO X 3 . A 101 ? 90 90 90\n6 C C1 EDO Y 3 . B 102 ? 90 90 90");
        let prepare = |s: String| {
            s.replace("1 polymer\n2 water", "1 polymer\n2 water\n3 non-polymer")
                .replace(
                    "_pdbx_struct_assembly.id\n1\n",
                    "_pdbx_struct_assembly.id\n1\n2\n",
                )
        };
        let old = prepare(old);
        let source = prepare(source);
        // The old/default mode must still refuse this conflict.
        assert!(repair_metadata(&old, &source, true, false).is_err());
        let fixed = repair_metadata_with_assemblies(&old, &source, true, false, true).unwrap();
        assert!(fixed.cif.contains("1 1 A,U,R0"));
        assert!(fixed.cif.contains("2 1 B,V,R1"));
        assert!(fixed.cif.contains("B 102 ? 16 17 18"));
        assert_eq!(fixed.water_atoms, 2);
        assert!(repair_metadata_with_assemblies(
            &old,
            &source.replace("B 102", "B 999"),
            true,
            false,
            true
        )
        .is_err());
        // One target chain must not merge ligands with different memberships.
        assert!(repair_metadata_with_assemblies(
            &old.replace("EDO V 3", "EDO U 3"),
            &source,
            true,
            false,
            true
        )
        .is_err());
    }
}
