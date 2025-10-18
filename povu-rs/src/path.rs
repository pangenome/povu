//! Path (reference walk) representation in the pangenome graph

/// Orientation of a vertex in a path
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(C)]
pub enum Orientation {
    /// Forward orientation (left/5' end)
    Forward = 0,
    /// Reverse orientation (right/3' end)
    Reverse = 1,
}

impl Orientation {
    /// Flip the orientation
    pub fn flip(&self) -> Self {
        match self {
            Orientation::Forward => Orientation::Reverse,
            Orientation::Reverse => Orientation::Forward,
        }
    }

    /// Convert to GFA orientation character
    pub fn to_gfa_char(&self) -> char {
        match self {
            Orientation::Forward => '+',
            Orientation::Reverse => '-',
        }
    }

    /// Parse from GFA orientation character
    pub fn from_gfa_char(c: char) -> Option<Self> {
        match c {
            '+' => Some(Orientation::Forward),
            '-' => Some(Orientation::Reverse),
            _ => None,
        }
    }
}

impl std::fmt::Display for Orientation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_gfa_char())
    }
}

/// A step in a path through the graph
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Step {
    /// Vertex ID at this step
    pub vertex_id: u64,
    /// Orientation of the vertex
    pub orientation: Orientation,
}

impl Step {
    /// Create a new step
    pub fn new(vertex_id: u64, orientation: Orientation) -> Self {
        Step {
            vertex_id,
            orientation,
        }
    }

    /// Get the reverse of this step
    pub fn reverse(&self) -> Self {
        Step {
            vertex_id: self.vertex_id,
            orientation: self.orientation.flip(),
        }
    }
}

impl std::fmt::Display for Step {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.vertex_id, self.orientation)
    }
}

/// A path (reference walk) through the pangenome graph
///
/// Represents a sequence of oriented vertices that forms a reference path.
/// In GFA terms, this corresponds to a P-line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Path {
    /// Path name (e.g., "sample#0#chr1" in PanSN format)
    pub name: String,
    /// Sequence of steps through the graph
    pub steps: Vec<Step>,
}

impl Path {
    /// Create a new path
    pub fn new(name: String, steps: Vec<Step>) -> Self {
        Path { name, steps }
    }

    /// Get the number of steps in the path
    pub fn len(&self) -> usize {
        self.steps.len()
    }

    /// Check if the path is empty
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }

    /// Parse sample name from PanSN format (sample#haplotype#contig)
    ///
    /// Returns the sample name portion before the first '#'
    pub fn sample_name(&self) -> Option<&str> {
        self.name.split('#').next()
    }

    /// Parse haplotype from PanSN format
    pub fn haplotype(&self) -> Option<&str> {
        let mut parts = self.name.split('#');
        parts.next()?; // skip sample
        parts.next()
    }

    /// Parse contig from PanSN format
    pub fn contig(&self) -> Option<&str> {
        let mut parts = self.name.split('#');
        parts.next()?; // skip sample
        parts.next()?; // skip haplotype
        parts.next()
    }

    /// Get an iterator over the vertex IDs in this path
    pub fn vertex_ids(&self) -> impl Iterator<Item = u64> + '_ {
        self.steps.iter().map(|s| s.vertex_id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orientation_flip() {
        assert_eq!(Orientation::Forward.flip(), Orientation::Reverse);
        assert_eq!(Orientation::Reverse.flip(), Orientation::Forward);
    }

    #[test]
    fn test_orientation_gfa() {
        assert_eq!(Orientation::Forward.to_gfa_char(), '+');
        assert_eq!(Orientation::Reverse.to_gfa_char(), '-');
        assert_eq!(Orientation::from_gfa_char('+'), Some(Orientation::Forward));
        assert_eq!(Orientation::from_gfa_char('-'), Some(Orientation::Reverse));
        assert_eq!(Orientation::from_gfa_char('x'), None);
    }

    #[test]
    fn test_step_creation() {
        let step = Step::new(42, Orientation::Forward);
        assert_eq!(step.vertex_id, 42);
        assert_eq!(step.orientation, Orientation::Forward);
        assert_eq!(format!("{}", step), "42+");
    }

    #[test]
    fn test_step_reverse() {
        let step = Step::new(42, Orientation::Forward);
        let rev = step.reverse();
        assert_eq!(rev.vertex_id, 42);
        assert_eq!(rev.orientation, Orientation::Reverse);
    }

    #[test]
    fn test_path_pansn_parsing() {
        let path = Path::new(
            "HG002#1#chr1".to_string(),
            vec![
                Step::new(1, Orientation::Forward),
                Step::new(2, Orientation::Reverse),
            ],
        );

        assert_eq!(path.sample_name(), Some("HG002"));
        assert_eq!(path.haplotype(), Some("1"));
        assert_eq!(path.contig(), Some("chr1"));
        assert_eq!(path.len(), 2);
    }

    #[test]
    fn test_path_vertex_ids() {
        let path = Path::new(
            "ref".to_string(),
            vec![
                Step::new(1, Orientation::Forward),
                Step::new(2, Orientation::Reverse),
                Step::new(3, Orientation::Forward),
            ],
        );

        let ids: Vec<u64> = path.vertex_ids().collect();
        assert_eq!(ids, vec![1, 2, 3]);
    }
}
