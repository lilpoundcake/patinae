//! Explicit offline repair, not a permissive structure reader.
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.len() != 5
        || !matches!(
            args[1].as_str(),
            "--water" | "--water-assemblies" | "--sequence"
        )
    {
        return Err(
            "usage: repair_cif_metadata --water|--water-assemblies|--sequence REDO.cif REFERENCE.cif OUTPUT.cif"
                .into(),
        );
    }
    let redo = std::fs::read_to_string(&args[2])?;
    let reference = std::fs::read_to_string(&args[3])?;
    let repaired = patinae_io::cif::repair::repair_metadata_with_assemblies(
        &redo,
        &reference,
        matches!(args[1].as_str(), "--water" | "--water-assemblies"),
        args[1] == "--sequence",
        args[1] == "--water-assemblies",
    )?;
    std::fs::write(&args[4], repaired.cif)?;
    println!(
        "water_atoms={} sequence_rows_added={}",
        repaired.water_atoms, repaired.sequence_rows_added
    );
    Ok(())
}
