//! # povu-rs: Rust bindings for Povu
//!
//! This crate provides both low-level FFI bindings and a high-level idiomatic Rust API
//! for Povu, a toolkit for exploring regions of genomic variation in pangenomes.
//!
//! ## Quick Start
//!
//! ```no_run
//! use povu::PovuGraph;
//!
//! // Simple one-shot conversion
//! povu::gfa_to_vcf("input.gfa", "output.vcf", None)?;
//!
//! // Or with more control
//! let graph = PovuGraph::load("input.gfa")?;
//! println!("Graph has {} vertices and {} edges",
//!          graph.vertex_count(), graph.edge_count());
//!
//! let analysis = graph.analyze()?;
//! analysis.write_vcf("output.vcf")?;
//! # Ok::<(), povu::Error>(())
//! ```
//!
//! ## Features
//!
//! - Load GFA files into bidirected pangenome graphs
//! - Query graph topology (vertices, edges, reference paths)
//! - Detect flubbles (regions of variation)
//! - Generate VCF variant calls
//! - Access hierarchical PVST (Pangenome Variation Structure Tree)

mod ffi;
mod error;
mod graph;
mod analysis;
mod vertex;
mod edge;
mod path;

pub use error::{Error, Result};
pub use graph::PovuGraph;
pub use analysis::GraphAnalysis;
pub use vertex::Vertex;
pub use edge::Edge;
pub use path::{Path, Step, Orientation};

/// One-shot convenience function to convert GFA to VCF
///
/// # Arguments
///
/// * `gfa_path` - Path to input GFA file
/// * `vcf_path` - Path to output VCF file
/// * `ref_file` - Optional path to file containing reference path names
///
/// # Example
///
/// ```no_run
/// povu::gfa_to_vcf("input.gfa", "output.vcf", None)?;
/// # Ok::<(), povu::Error>(())
/// ```
pub fn gfa_to_vcf(
    gfa_path: impl AsRef<std::path::Path>,
    vcf_path: impl AsRef<std::path::Path>,
    ref_file: Option<impl AsRef<std::path::Path>>,
) -> Result<()> {
    use std::ffi::CString;

    let gfa_path_c = CString::new(gfa_path.as_ref().to_string_lossy().as_ref())?;
    let vcf_path_c = CString::new(vcf_path.as_ref().to_string_lossy().as_ref())?;
    let ref_file_c = ref_file
        .as_ref()
        .map(|p| CString::new(p.as_ref().to_string_lossy().as_ref()))
        .transpose()?;

    let mut error = ffi::PovuError::default();

    let success = unsafe {
        ffi::povu_gfa_to_vcf(
            gfa_path_c.as_ptr(),
            vcf_path_c.as_ptr(),
            ref_file_c.as_ref().map_or(std::ptr::null(), |c| c.as_ptr()),
            &mut error as *mut _,
        )
    };

    if !success {
        return Err(Error::from_ffi_error(error));
    }

    Ok(())
}
