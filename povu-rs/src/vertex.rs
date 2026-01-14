//! Vertex (segment) representation in the pangenome graph

/// A vertex in the bidirected pangenome graph
///
/// Each vertex represents a segment in the GFA file, with an ID and sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Vertex {
    /// Unique vertex ID
    pub id: u64,
    /// DNA sequence for this vertex
    pub sequence: String,
}

impl Vertex {
    /// Get the length of the sequence
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if the vertex has an empty sequence
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get the reverse complement of this vertex's sequence
    pub fn reverse_complement(&self) -> String {
        reverse_complement(&self.sequence)
    }
}

/// Compute reverse complement of a DNA sequence
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => c,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ACGT"), "ACGT");
        assert_eq!(reverse_complement("AAAA"), "TTTT");
        assert_eq!(reverse_complement("GCGC"), "GCGC");
        assert_eq!(reverse_complement("ATCG"), "CGAT");
    }

    #[test]
    fn test_vertex_rc() {
        let v = Vertex {
            id: 1,
            sequence: "ACGT".to_string(),
        };
        assert_eq!(v.reverse_complement(), "ACGT");
    }
}
