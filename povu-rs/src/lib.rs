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
//! povu::gfa_to_vcf("input.gfa", "output.vcf", None::<&std::path::Path>)?;
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

mod analysis;
mod edge;
mod error;
mod ffi;
mod graph;
pub mod native_gfa;
mod path;
pub mod vcf;
mod vertex;

pub use analysis::GraphAnalysis;
pub use edge::Edge;
pub use error::{Error, Result};
pub use graph::PovuGraph;
pub use native_gfa::{
    detect_flubble_stack, gfa_to_vcf_document, FlubbleBoundary, FlubbleCandidate, NativeGfa,
};
pub use path::{Orientation, Path, Step};
pub use vcf::{
    AlleleConstruction, AlleleFrequency, AlternateAllele, Contig, GenotypeAllele, GenotypeColumn,
    Info as VcfInfo, OrderKey as VcfOrderKey, Record as VcfRecord, VariantCall, VariantSource,
    VariantType, VcfDocument,
};
pub use vertex::Vertex;

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
/// povu::gfa_to_vcf("input.gfa", "output.vcf", None::<&std::path::Path>)?;
/// # Ok::<(), povu::Error>(())
/// ```
pub fn gfa_to_vcf(
    gfa_path: impl AsRef<std::path::Path>,
    vcf_path: impl AsRef<std::path::Path>,
    ref_file: Option<impl AsRef<std::path::Path>>,
) -> Result<()> {
    let document = native_gfa::gfa_to_vcf_document(gfa_path, ref_file)?;
    document.write_path(vcf_path)
}
