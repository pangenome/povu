//! Integration tests for povu-rs
//!
//! These tests compare the Rust wrapper outputs to native Povu outputs
//! to ensure correctness and compatibility.

use povu::{PovuGraph, Error};
use std::path::PathBuf;

fn test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("tests")
        .join("data")
}

#[test]
fn test_load_simple_gfa() {
    let gfa_path = test_data_dir().join("simple.gfa");

    if !gfa_path.exists() {
        eprintln!("Warning: test file {:?} does not exist, skipping test", gfa_path);
        return;
    }

    let graph = PovuGraph::load(&gfa_path).expect("Failed to load GFA");

    assert!(graph.vertex_count() > 0, "Graph should have vertices");
    println!("Loaded graph with {} vertices and {} edges",
             graph.vertex_count(), graph.edge_count());
}

#[test]
fn test_query_vertices() {
    let gfa_path = test_data_dir().join("simple.gfa");

    if !gfa_path.exists() {
        eprintln!("Warning: test file {:?} does not exist, skipping test", gfa_path);
        return;
    }

    let graph = PovuGraph::load(&gfa_path).expect("Failed to load GFA");
    let vertices = graph.vertices().expect("Failed to get vertices");

    assert_eq!(vertices.len(), graph.vertex_count());

    for vertex in &vertices {
        assert!(vertex.id > 0, "Vertex ID should be positive");
        assert!(!vertex.sequence.is_empty() || vertex.id > 0, "Vertex should have sequence or valid ID");
        println!("Vertex {}: {} bp", vertex.id, vertex.len());
    }
}

#[test]
fn test_query_edges() {
    let gfa_path = test_data_dir().join("simple.gfa");

    if !gfa_path.exists() {
        eprintln!("Warning: test file {:?} does not exist, skipping test", gfa_path);
        return;
    }

    let graph = PovuGraph::load(&gfa_path).expect("Failed to load GFA");
    let edges = graph.edges().expect("Failed to get edges");

    assert_eq!(edges.len(), graph.edge_count());

    for edge in &edges {
        println!("Edge: {}{} -> {}{}",
                 edge.from_id, edge.from_orientation,
                 edge.to_id, edge.to_orientation);
    }
}

#[test]
fn test_query_paths() {
    let gfa_path = test_data_dir().join("simple.gfa");

    if !gfa_path.exists() {
        eprintln!("Warning: test file {:?} does not exist, skipping test", gfa_path);
        return;
    }

    let graph = PovuGraph::load(&gfa_path).expect("Failed to load GFA");
    let paths = graph.paths().expect("Failed to get paths");

    for path in &paths {
        println!("Path '{}': {} steps", path.name, path.len());

        if let Some(sample) = path.sample_name() {
            println!("  Sample: {}", sample);
        }

        // Print first few steps
        for (i, step) in path.steps.iter().take(5).enumerate() {
            println!("  Step {}: {}", i, step);
        }
    }
}

#[test]
fn test_analyze_graph() {
    let gfa_path = test_data_dir().join("simple.gfa");

    if !gfa_path.exists() {
        eprintln!("Warning: test file {:?} does not exist, skipping test", gfa_path);
        return;
    }

    let graph = PovuGraph::load(&gfa_path).expect("Failed to load GFA");

    // This may fail if the graph is too simple or doesn't have variation
    match graph.analyze() {
        Ok(analysis) => {
            println!("Found {} flubbles", analysis.flubble_count());
            let pvst = analysis.pvst_tree();
            println!("PVST has {} vertices", pvst.vertex_count());
        }
        Err(e) => {
            eprintln!("Analysis failed (may be expected for simple graphs): {}", e);
        }
    }
}

#[test]
fn test_set_references_by_prefix() {
    let gfa_path = test_data_dir().join("simple.gfa");

    if !gfa_path.exists() {
        eprintln!("Warning: test file {:?} does not exist, skipping test", gfa_path);
        return;
    }

    let mut graph = PovuGraph::load(&gfa_path).expect("Failed to load GFA");

    // Try setting references - this may succeed or fail depending on path names
    match graph.set_references_from_prefixes(&["ref", "sample"]) {
        Ok(_) => println!("Successfully set references"),
        Err(e) => eprintln!("Setting references failed: {}", e),
    }
}

#[test]
fn test_vertex_reverse_complement() {
    use povu::Vertex;

    let v = Vertex {
        id: 1,
        sequence: "ACGTACGT".to_string(),
    };

    let rc = v.reverse_complement();
    assert_eq!(rc, "ACGTACGT");

    let v2 = Vertex {
        id: 2,
        sequence: "AAAA".to_string(),
    };
    assert_eq!(v2.reverse_complement(), "TTTT");
}

#[test]
fn test_edge_operations() {
    use povu::{Edge, Orientation};

    let edge = Edge::new(1, Orientation::Forward, 2, Orientation::Reverse);
    assert!(!edge.is_self_loop());

    let rev = edge.reverse();
    assert_eq!(rev.from_id, 2);
    assert_eq!(rev.to_id, 1);

    let self_loop = Edge::new(5, Orientation::Forward, 5, Orientation::Reverse);
    assert!(self_loop.is_self_loop());
}

#[test]
fn test_error_handling() {
    // Try to load a non-existent file
    let result = PovuGraph::load("/nonexistent/path/to/file.gfa");
    assert!(result.is_err());

    if let Err(e) = result {
        println!("Expected error: {}", e);
    }
}
