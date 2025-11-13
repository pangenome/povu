//! Tests for in-memory graph building API

use povu::{PovuGraph, Orientation};

#[test]
fn test_create_empty_graph() {
    let graph = PovuGraph::new(10, 10, 0);
    assert_eq!(graph.vertex_count(), 0);
    assert_eq!(graph.edge_count(), 0);
}

#[test]
fn test_add_single_vertex() {
    let mut graph = PovuGraph::new(10, 10, 0);

    let idx = graph.add_vertex(1, "ACGT").expect("Failed to add vertex");
    assert_eq!(idx, 0, "First vertex should have index 0");

    assert_eq!(graph.vertex_count(), 1);
}

#[test]
fn test_add_multiple_vertices() {
    let mut graph = PovuGraph::new(10, 10, 0);

    graph.add_vertex(1, "AAAA").expect("Failed to add vertex 1");
    graph.add_vertex(2, "CCCC").expect("Failed to add vertex 2");
    graph.add_vertex(3, "GGGG").expect("Failed to add vertex 3");
    graph.add_vertex(4, "TTTT").expect("Failed to add vertex 4");

    assert_eq!(graph.vertex_count(), 4);
}

#[test]
fn test_add_edge() {
    let mut graph = PovuGraph::new(10, 10, 0);

    graph.add_vertex(1, "AAAA").unwrap();
    graph.add_vertex(2, "CCCC").unwrap();

    let edge_idx = graph.add_edge(
        1, Orientation::Forward,
        2, Orientation::Forward
    ).expect("Failed to add edge");

    assert_eq!(edge_idx, 0, "First edge should have index 0");
    assert_eq!(graph.edge_count(), 1);
}

#[test]
fn test_build_simple_path() {
    let mut graph = PovuGraph::new(10, 10, 0);

    // Build a simple linear path: 1 -> 2 -> 3
    graph.add_vertex(1, "AAA").unwrap();
    graph.add_vertex(2, "CCC").unwrap();
    graph.add_vertex(3, "GGG").unwrap();

    graph.add_edge(1, Orientation::Forward, 2, Orientation::Forward).unwrap();
    graph.add_edge(2, Orientation::Forward, 3, Orientation::Forward).unwrap();

    assert_eq!(graph.vertex_count(), 3);
    assert_eq!(graph.edge_count(), 2);
}

#[test]
fn test_build_diamond_graph() {
    let mut graph = PovuGraph::new(10, 10, 0);

    // Build diamond: 1 -> 2 -> 4
    //                  \-> 3 ->/
    graph.add_vertex(1, "AAAA").unwrap();
    graph.add_vertex(2, "GGGG").unwrap();
    graph.add_vertex(3, "TTTT").unwrap();
    graph.add_vertex(4, "CCCC").unwrap();

    graph.add_edge(1, Orientation::Forward, 2, Orientation::Forward).unwrap();
    graph.add_edge(1, Orientation::Forward, 3, Orientation::Forward).unwrap();
    graph.add_edge(2, Orientation::Forward, 4, Orientation::Forward).unwrap();
    graph.add_edge(3, Orientation::Forward, 4, Orientation::Forward).unwrap();

    assert_eq!(graph.vertex_count(), 4);
    assert_eq!(graph.edge_count(), 4);
}

#[test]
fn test_finalize_graph() {
    let mut graph = PovuGraph::new(10, 10, 0);

    graph.add_vertex(1, "ACGT").unwrap();
    graph.add_vertex(2, "TGCA").unwrap();
    graph.add_edge(1, Orientation::Forward, 2, Orientation::Forward).unwrap();

    // Should not panic
    graph.finalize();

    // Graph should still be queryable
    assert_eq!(graph.vertex_count(), 2);
    assert_eq!(graph.edge_count(), 1);
}

#[test]
fn test_query_built_vertices() {
    let mut graph = PovuGraph::new(10, 10, 0);

    graph.add_vertex(10, "AAAA").unwrap();
    graph.add_vertex(20, "CCCC").unwrap();
    graph.add_vertex(30, "GGGG").unwrap();

    graph.finalize();

    let vertices = graph.vertices().expect("Failed to get vertices");

    assert_eq!(vertices.len(), 3);

    // Check IDs
    let ids: Vec<u64> = vertices.iter().map(|v| v.id).collect();
    assert!(ids.contains(&10));
    assert!(ids.contains(&20));
    assert!(ids.contains(&30));

    // Check sequences
    for v in &vertices {
        match v.id {
            10 => assert_eq!(v.sequence, "AAAA"),
            20 => assert_eq!(v.sequence, "CCCC"),
            30 => assert_eq!(v.sequence, "GGGG"),
            _ => panic!("Unexpected vertex ID: {}", v.id),
        }
    }
}

#[test]
fn test_query_built_edges() {
    let mut graph = PovuGraph::new(10, 10, 0);

    graph.add_vertex(1, "A").unwrap();
    graph.add_vertex(2, "C").unwrap();
    graph.add_vertex(3, "G").unwrap();

    graph.add_edge(1, Orientation::Forward, 2, Orientation::Forward).unwrap();
    graph.add_edge(2, Orientation::Reverse, 3, Orientation::Forward).unwrap();

    graph.finalize();

    let edges = graph.edges().expect("Failed to get edges");

    assert_eq!(edges.len(), 2);

    // Check edge 1 -> 2
    let edge1 = edges.iter().find(|e| e.from_id == 1).expect("Edge 1->2 not found");
    assert_eq!(edge1.from_id, 1);
    assert_eq!(edge1.from_orientation, Orientation::Forward);
    assert_eq!(edge1.to_id, 2);
    assert_eq!(edge1.to_orientation, Orientation::Forward);

    // Check edge 2 -> 3
    let edge2 = edges.iter().find(|e| e.from_id == 2).expect("Edge 2->3 not found");
    assert_eq!(edge2.from_id, 2);
    assert_eq!(edge2.from_orientation, Orientation::Reverse);
    assert_eq!(edge2.to_id, 3);
    assert_eq!(edge2.to_orientation, Orientation::Forward);
}

#[test]
fn test_bidirected_orientations() {
    let mut graph = PovuGraph::new(10, 10, 0);

    graph.add_vertex(1, "ACGT").unwrap();
    graph.add_vertex(2, "TGCA").unwrap();

    // Test all orientation combinations
    graph.add_edge(1, Orientation::Forward, 2, Orientation::Forward).unwrap();
    graph.add_edge(1, Orientation::Forward, 2, Orientation::Reverse).unwrap();
    graph.add_edge(1, Orientation::Reverse, 2, Orientation::Forward).unwrap();
    graph.add_edge(1, Orientation::Reverse, 2, Orientation::Reverse).unwrap();

    assert_eq!(graph.edge_count(), 4);
}

#[test]
fn test_empty_sequence() {
    let mut graph = PovuGraph::new(10, 10, 0);

    // Should handle empty sequences
    let result = graph.add_vertex(1, "");
    assert!(result.is_ok(), "Should accept empty sequence");
}

#[test]
fn test_long_sequence() {
    let mut graph = PovuGraph::new(10, 10, 0);

    // Test with a longer sequence
    let long_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let result = graph.add_vertex(1, long_seq);
    assert!(result.is_ok(), "Should handle long sequences");

    graph.finalize();

    let vertices = graph.vertices().unwrap();
    assert_eq!(vertices[0].sequence, long_seq);
    assert_eq!(vertices[0].len(), 40);
}

#[test]
fn test_self_loop() {
    let mut graph = PovuGraph::new(10, 10, 0);

    graph.add_vertex(1, "ACGT").unwrap();

    // Self-loop: vertex connects to itself
    let result = graph.add_edge(1, Orientation::Forward, 1, Orientation::Reverse);
    assert!(result.is_ok(), "Should allow self-loops");

    assert_eq!(graph.edge_count(), 1);
}

#[test]
fn test_comparison_with_loaded_graph() {
    use std::path::PathBuf;

    let test_data = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("tests")
        .join("data")
        .join("simple.gfa");

    if !test_data.exists() {
        eprintln!("Warning: Test data not found, skipping comparison test");
        return;
    }

    // Load graph from file
    let loaded_graph = PovuGraph::load(&test_data).expect("Failed to load GFA");

    // Build equivalent graph in memory
    let mut built_graph = PovuGraph::new(
        loaded_graph.vertex_count(),
        loaded_graph.edge_count(),
        0
    );

    // Copy vertices
    let loaded_vertices = loaded_graph.vertices().unwrap();
    for v in &loaded_vertices {
        built_graph.add_vertex(v.id, &v.sequence).unwrap();
    }

    // Copy edges
    let loaded_edges = loaded_graph.edges().unwrap();
    for e in &loaded_edges {
        built_graph.add_edge(
            e.from_id, e.from_orientation,
            e.to_id, e.to_orientation
        ).unwrap();
    }

    built_graph.finalize();

    // Compare
    assert_eq!(built_graph.vertex_count(), loaded_graph.vertex_count());
    assert_eq!(built_graph.edge_count(), loaded_graph.edge_count());
}
