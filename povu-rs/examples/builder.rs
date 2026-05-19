//! Example: Build a graph in memory without using GFA files

use povu::{PovuGraph, Orientation};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Building a simple graph in memory...\n");

    // Create an empty graph with capacity hints
    let mut graph = PovuGraph::new(
        10,  // vertices
        12,  // edges
        0,   // paths (not using paths in this example)
    );

    // Build a simple diamond graph:
    //
    //      2 (GGGG)
    //     / \
    //  1 +   + 4
    //     \ /
    //      3 (TTTT)
    //
    // This creates a bubble with two alternative paths

    println!("Adding vertices...");
    graph.add_vertex(1, "AAAA")?;
    graph.add_vertex(2, "GGGG")?; // Alternative 1
    graph.add_vertex(3, "TTTT")?; // Alternative 2
    graph.add_vertex(4, "CCCC")?;

    println!("  Vertex 1: AAAA");
    println!("  Vertex 2: GGGG (alt 1)");
    println!("  Vertex 3: TTTT (alt 2)");
    println!("  Vertex 4: CCCC");

    println!("\nAdding edges...");
    // 1 -> 2
    graph.add_edge(1, Orientation::Forward, 2, Orientation::Forward)?;
    println!("  1+ -> 2+");

    // 1 -> 3
    graph.add_edge(1, Orientation::Forward, 3, Orientation::Forward)?;
    println!("  1+ -> 3+");

    // 2 -> 4
    graph.add_edge(2, Orientation::Forward, 4, Orientation::Forward)?;
    println!("  2+ -> 4+");

    // 3 -> 4
    graph.add_edge(3, Orientation::Forward, 4, Orientation::Forward)?;
    println!("  3+ -> 4+");

    println!("\nFinalizing graph...");
    graph.finalize();

    println!("\n=== Graph Statistics ===");
    println!("Vertices: {}", graph.vertex_count());
    println!("Edges:    {}", graph.edge_count());

    // Query the graph
    println!("\n=== Querying Graph ===");
    let vertices = graph.vertices()?;
    println!("All vertices:");
    for v in &vertices {
        println!("  {}: {} ({}bp)", v.id, v.sequence, v.len());
    }

    let edges = graph.edges()?;
    println!("\nAll edges:");
    for e in &edges {
        println!("  {}{} -> {}{}",
                 e.from_id, e.from_orientation,
                 e.to_id, e.to_orientation);
    }

    println!("\n=== Analyzing for Variation ===");
    match graph.analyze() {
        Ok(analysis) => {
            println!("Successfully analyzed graph!");
            println!("Found {} variation regions", analysis.flubble_count());

            let pvst = analysis.pvst_tree();
            println!("PVST tree has {} vertices", pvst.vertex_count());
        }
        Err(e) => {
            println!("Analysis failed: {}", e);
            println!("(This is expected for very simple graphs)");
        }
    }

    Ok(())
}
