//! Graph analysis and variant calling functionality

use crate::{ffi, Error, Result};
use std::path::Path;

/// Result of analyzing a graph for variation regions
///
/// This structure provides access to the detected flubbles (variation regions)
/// and the hierarchical PVST tree. It can also generate VCF output.
pub struct GraphAnalysis {
    pub(crate) inner: *mut ffi::PovuFlubbles,
}

impl GraphAnalysis {
    /// Get the number of flubbles (variation regions) detected
    pub fn flubble_count(&self) -> usize {
        unsafe { ffi::povu_flubbles_count(self.inner) }
    }

    /// Get the PVST tree representing the hierarchical structure
    pub fn pvst_tree(&self) -> PvstTree {
        let tree_ptr = unsafe { ffi::povu_flubbles_get_pvst_tree(self.inner) };
        PvstTree { inner: tree_ptr }
    }

    /// Write VCF output to a file
    ///
    /// # Arguments
    ///
    /// * `path` - Output path for the VCF file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use povu::PovuGraph;
    ///
    /// let graph = PovuGraph::load("test.gfa")?;
    /// let analysis = graph.analyze()?;
    /// analysis.write_vcf("output.vcf")?;
    /// # Ok::<(), povu::Error>(())
    /// ```
    pub fn write_vcf(&self, _path: impl AsRef<Path>) -> Result<()> {
        // TODO: Implement VCF generation when FFI is complete
        Err(Error::Povu {
            code: 1,
            message: "VCF generation not yet implemented".to_string(),
        })
    }
}

impl Drop for GraphAnalysis {
    fn drop(&mut self) {
        if !self.inner.is_null() {
            unsafe {
                ffi::povu_flubbles_free(self.inner);
            }
        }
    }
}

/// Pangenome Variation Structure Tree (PVST)
///
/// Represents the hierarchical nesting relationships between variation regions.
pub struct PvstTree {
    inner: *mut ffi::PovuPvstTree,
}

impl PvstTree {
    /// Get the number of vertices in the PVST tree
    pub fn vertex_count(&self) -> usize {
        if self.inner.is_null() {
            return 0;
        }
        unsafe { ffi::povu_pvst_tree_vertex_count(self.inner) }
    }

    // TODO: Add methods for tree traversal and querying specific nodes
}

impl Drop for PvstTree {
    fn drop(&mut self) {
        if !self.inner.is_null() {
            unsafe {
                ffi::povu_pvst_tree_free(self.inner);
            }
        }
    }
}

// Safety: GraphAnalysis can be sent between threads
unsafe impl Send for GraphAnalysis {}
unsafe impl Send for PvstTree {}
