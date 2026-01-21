//! Edge representation in the bidirected pangenome graph

use crate::Orientation;

/// An edge in the bidirected pangenome graph
///
/// Each edge connects two vertices at specific orientations (ends).
/// In a bidirected graph, each segment has two ends: forward (left/5') and reverse (right/3').
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Edge {
    /// Source vertex ID
    pub from_id: u64,
    /// Source vertex orientation/end
    pub from_orientation: Orientation,
    /// Target vertex ID
    pub to_id: u64,
    /// Target vertex orientation/end
    pub to_orientation: Orientation,
}

impl Edge {
    /// Create a new edge
    pub fn new(
        from_id: u64,
        from_orientation: Orientation,
        to_id: u64,
        to_orientation: Orientation,
    ) -> Self {
        Edge {
            from_id,
            from_orientation,
            to_id,
            to_orientation,
        }
    }

    /// Check if this is a self-loop (edge from a vertex to itself)
    pub fn is_self_loop(&self) -> bool {
        self.from_id == self.to_id
    }

    /// Get the reverse of this edge
    pub fn reverse(&self) -> Self {
        Edge {
            from_id: self.to_id,
            from_orientation: self.to_orientation.flip(),
            to_id: self.from_id,
            to_orientation: self.from_orientation.flip(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edge_creation() {
        let edge = Edge::new(1, Orientation::Forward, 2, Orientation::Reverse);
        assert_eq!(edge.from_id, 1);
        assert_eq!(edge.to_id, 2);
        assert!(!edge.is_self_loop());
    }

    #[test]
    fn test_self_loop() {
        let edge = Edge::new(1, Orientation::Forward, 1, Orientation::Reverse);
        assert!(edge.is_self_loop());
    }

    #[test]
    fn test_edge_reverse() {
        let edge = Edge::new(1, Orientation::Forward, 2, Orientation::Reverse);
        let rev = edge.reverse();
        assert_eq!(rev.from_id, 2);
        assert_eq!(rev.to_id, 1);
        assert_eq!(rev.from_orientation, Orientation::Forward);
        assert_eq!(rev.to_orientation, Orientation::Reverse);
    }
}
