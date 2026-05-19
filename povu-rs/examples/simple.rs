//! Simple example: Load a GFA file and print basic statistics

use povu::PovuGraph;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <input.gfa>", args[0]);
        std::process::exit(1);
    }

    let gfa_path = &args[1];

    println!("Loading graph from {}...", gfa_path);
    let graph = PovuGraph::load(gfa_path)?;

    println!("\nGraph Statistics:");
    println!("  Vertices: {}", graph.vertex_count());
    println!("  Edges:    {}", graph.edge_count());
    println!("  Paths:    {}", graph.path_count());

    // Show vertices
    println!("\nVertices:");
    let vertices = graph.vertices()?;
    for (i, vertex) in vertices.iter().take(10).enumerate() {
        println!("  {}: id={}, len={} bp",
                 i, vertex.id, vertex.len());
    }
    if vertices.len() > 10 {
        println!("  ... and {} more", vertices.len() - 10);
    }

    // Show edges
    println!("\nEdges:");
    let edges = graph.edges()?;
    for (i, edge) in edges.iter().take(10).enumerate() {
        println!("  {}: {}{} -> {}{}",
                 i,
                 edge.from_id, edge.from_orientation,
                 edge.to_id, edge.to_orientation);
    }
    if edges.len() > 10 {
        println!("  ... and {} more", edges.len() - 10);
    }

    // Show paths
    println!("\nPaths:");
    let paths = graph.paths()?;
    for path in &paths {
        println!("  '{}': {} steps", path.name, path.len());
        if let Some(sample) = path.sample_name() {
            println!("    Sample: {}", sample);
            if let Some(hap) = path.haplotype() {
                println!("    Haplotype: {}", hap);
            }
            if let Some(contig) = path.contig() {
                println!("    Contig: {}", contig);
            }
        }
    }

    Ok(())
}
