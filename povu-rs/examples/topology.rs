//! Topology example: Analyze graph structure and variation regions

use povu::PovuGraph;
use std::collections::{HashMap, HashSet};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <input.gfa>", args[0]);
        std::process::exit(1);
    }

    let gfa_path = &args[1];

    println!("Loading graph from {}...", gfa_path);
    let graph = PovuGraph::load(gfa_path)?;

    println!("\n=== Basic Statistics ===");
    println!("Vertices: {}", graph.vertex_count());
    println!("Edges:    {}", graph.edge_count());
    println!("Paths:    {}", graph.path_count());

    // Analyze vertex degree distribution
    println!("\n=== Degree Distribution ===");
    let edges = graph.edges()?;
    let mut degrees: HashMap<u64, usize> = HashMap::new();

    for edge in &edges {
        *degrees.entry(edge.from_id).or_insert(0) += 1;
        *degrees.entry(edge.to_id).or_insert(0) += 1;
    }

    let max_degree = degrees.values().max().copied().unwrap_or(0);
    let min_degree = degrees.values().min().copied().unwrap_or(0);
    let avg_degree = if !degrees.is_empty() {
        degrees.values().sum::<usize>() as f64 / degrees.len() as f64
    } else {
        0.0
    };

    println!("  Min degree: {}", min_degree);
    println!("  Max degree: {}", max_degree);
    println!("  Avg degree: {:.2}", avg_degree);

    // Find tips (vertices with degree 1)
    let tips: Vec<u64> = degrees.iter()
        .filter(|(_, &d)| d == 1)
        .map(|(&id, _)| id)
        .collect();
    println!("  Tips: {}", tips.len());

    // Analyze path coverage
    println!("\n=== Path Coverage ===");
    let paths = graph.paths()?;
    let vertices = graph.vertices()?;

    let total_vertices = vertices.len();
    let mut covered_vertices = HashSet::new();

    for path in &paths {
        for step in &path.steps {
            covered_vertices.insert(step.vertex_id);
        }
    }

    let coverage_pct = if total_vertices > 0 {
        (covered_vertices.len() as f64 / total_vertices as f64) * 100.0
    } else {
        0.0
    };

    println!("  Vertices covered by paths: {}/{} ({:.1}%)",
             covered_vertices.len(), total_vertices, coverage_pct);

    // Find vertices not on any path
    let uncovered: Vec<&u64> = vertices.iter()
        .map(|v| &v.id)
        .filter(|id| !covered_vertices.contains(id))
        .take(10)
        .collect();

    if !uncovered.is_empty() {
        println!("  Uncovered vertices (first 10): {:?}", uncovered);
    }

    // Analyze path lengths
    println!("\n=== Path Lengths ===");
    for path in &paths {
        let total_bp: usize = path.steps.iter()
            .filter_map(|step| {
                vertices.iter()
                    .find(|v| v.id == step.vertex_id)
                    .map(|v| v.len())
            })
            .sum();

        println!("  '{}': {} steps, ~{} bp",
                 path.name, path.len(), total_bp);
    }

    // Try to analyze flubbles
    println!("\n=== Variation Analysis ===");
    match graph.analyze() {
        Ok(analysis) => {
            println!("  Flubbles detected: {}", analysis.flubble_count());

            let pvst = analysis.pvst_tree();
            println!("  PVST vertices: {}", pvst.vertex_count());

            println!("\nNote: Use GraphAnalysis for detailed variant calling");
        }
        Err(e) => {
            println!("  Analysis failed: {}", e);
            println!("  (This may be expected for simple or linear graphs)");
        }
    }

    // Detect simple bubbles manually
    println!("\n=== Simple Bubble Detection ===");
    detect_simple_bubbles(&edges, &vertices);

    Ok(())
}

/// Detect simple 2-vertex bubbles in the graph
fn detect_simple_bubbles(
    edges: &[povu::Edge],
    vertices: &[povu::Vertex],
) {
    use std::collections::HashMap;

    // Build adjacency list
    let mut adj: HashMap<u64, Vec<u64>> = HashMap::new();

    for edge in edges {
        adj.entry(edge.from_id).or_default().push(edge.to_id);
    }

    let mut bubble_count = 0;

    // Look for vertices with exactly 2 outgoing edges
    for (source, targets) in &adj {
        if targets.len() != 2 {
            continue;
        }

        let t1 = targets[0];
        let t2 = targets[1];

        // Check if both targets converge to the same vertex
        if let (Some(t1_targets), Some(t2_targets)) = (adj.get(&t1), adj.get(&t2)) {
            let common: HashSet<_> = t1_targets.iter()
                .filter(|t| t2_targets.contains(t))
                .collect();

            if !common.is_empty() {
                bubble_count += 1;
                if bubble_count <= 5 {
                    let v1_len = vertices.iter().find(|v| v.id == t1).map(|v| v.len()).unwrap_or(0);
                    let v2_len = vertices.iter().find(|v| v.id == t2).map(|v| v.len()).unwrap_or(0);

                    println!("  Bubble: {} -> [{} ({} bp), {} ({} bp)] -> {}",
                             source, t1, v1_len, t2, v2_len,
                             common.iter().next().unwrap());
                }
            }
        }
    }

    println!("  Total simple bubbles found: {}", bubble_count);
}
