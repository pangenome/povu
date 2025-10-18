//! Main graph API for loading and querying pangenome graphs

use crate::{ffi, Error, GraphAnalysis, Result, Vertex, Edge, Path, Step, Orientation};
use std::ffi::{CStr, CString};
use std::path::Path as StdPath;

/// A pangenome variation graph loaded from GFA
///
/// This is the main entry point for working with Povu graphs.
/// It provides methods to query the graph topology and perform analysis.
pub struct PovuGraph {
    pub(crate) inner: *mut ffi::PovuGraph,
}

impl PovuGraph {
    /// Load a pangenome graph from a GFA file
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the GFA file
    ///
    /// # Example
    ///
    /// ```no_run
    /// use povu::PovuGraph;
    ///
    /// let graph = PovuGraph::load("test.gfa")?;
    /// println!("Loaded graph with {} vertices", graph.vertex_count());
    /// # Ok::<(), povu::Error>(())
    /// ```
    pub fn load(path: impl AsRef<StdPath>) -> Result<Self> {
        let path_str = path.as_ref().to_string_lossy();
        let c_path = CString::new(path_str.as_ref())?;
        let mut error = ffi::PovuError::default();

        let inner = unsafe { ffi::povu_graph_from_gfa(c_path.as_ptr(), &mut error as *mut _) };

        if inner.is_null() {
            return Err(Error::from_ffi_error(error));
        }

        Ok(PovuGraph { inner })
    }

    /// Get the number of vertices in the graph
    pub fn vertex_count(&self) -> usize {
        unsafe { ffi::povu_graph_vertex_count(self.inner) }
    }

    /// Get the number of edges in the graph
    pub fn edge_count(&self) -> usize {
        unsafe { ffi::povu_graph_edge_count(self.inner) }
    }

    /// Get the number of reference paths in the graph
    pub fn path_count(&self) -> usize {
        unsafe { ffi::povu_graph_path_count(self.inner) }
    }

    /// Get all vertices in the graph
    ///
    /// Returns a vector of vertices with their IDs and sequences.
    pub fn vertices(&self) -> Result<Vec<Vertex>> {
        let mut count: usize = 0;
        let vertices_ptr = unsafe { ffi::povu_graph_get_vertices(self.inner, &mut count as *mut _) };

        if vertices_ptr.is_null() {
            return Ok(Vec::new());
        }

        let mut vertices = Vec::with_capacity(count);
        for i in 0..count {
            let v = unsafe { &*vertices_ptr.add(i) };
            let sequence = if v.sequence.is_null() {
                String::new()
            } else {
                unsafe {
                    CStr::from_ptr(v.sequence)
                        .to_string_lossy()
                        .into_owned()
                }
            };

            vertices.push(Vertex {
                id: v.id,
                sequence,
            });
        }

        unsafe { ffi::povu_vertices_free(vertices_ptr, count) };

        Ok(vertices)
    }

    /// Get all edges in the graph
    ///
    /// Returns a vector of edges with their endpoints and orientations.
    pub fn edges(&self) -> Result<Vec<Edge>> {
        let mut count: usize = 0;
        let edges_ptr = unsafe { ffi::povu_graph_get_edges(self.inner, &mut count as *mut _) };

        if edges_ptr.is_null() {
            return Ok(Vec::new());
        }

        let mut edges = Vec::with_capacity(count);
        for i in 0..count {
            let e = unsafe { &*edges_ptr.add(i) };
            edges.push(Edge {
                from_id: e.from_id,
                from_orientation: unsafe { std::mem::transmute(e.from_orientation) },
                to_id: e.to_id,
                to_orientation: unsafe { std::mem::transmute(e.to_orientation) },
            });
        }

        unsafe { ffi::povu_edges_free(edges_ptr, count) };

        Ok(edges)
    }

    /// Get all reference paths in the graph
    ///
    /// Returns a vector of paths with their names and steps.
    pub fn paths(&self) -> Result<Vec<Path>> {
        let mut count: usize = 0;
        let paths_ptr = unsafe { ffi::povu_graph_get_paths(self.inner, &mut count as *mut _) };

        if paths_ptr.is_null() {
            return Ok(Vec::new());
        }

        let mut paths = Vec::with_capacity(count);
        for i in 0..count {
            let p = unsafe { &*paths_ptr.add(i) };

            let name = if p.name.is_null() {
                String::new()
            } else {
                unsafe { CStr::from_ptr(p.name).to_string_lossy().into_owned() }
            };

            let mut steps = Vec::with_capacity(p.steps_count);
            if !p.steps.is_null() {
                for j in 0..p.steps_count {
                    let s = unsafe { &*p.steps.add(j) };
                    steps.push(Step {
                        vertex_id: s.vertex_id,
                        orientation: unsafe { std::mem::transmute(s.orientation) },
                    });
                }
            }

            paths.push(Path { name, steps });
        }

        unsafe { ffi::povu_paths_free(paths_ptr, count) };

        Ok(paths)
    }

    /// Set reference paths from a file
    ///
    /// The file should contain one reference path name per line.
    ///
    /// # Arguments
    ///
    /// * `path` - Path to the file containing reference names
    pub fn set_references_from_file(&mut self, path: impl AsRef<StdPath>) -> Result<()> {
        let path_str = path.as_ref().to_string_lossy();
        let c_path = CString::new(path_str.as_ref())?;
        let mut error = ffi::PovuError::default();

        let success = unsafe {
            ffi::povu_graph_set_references_from_file(self.inner, c_path.as_ptr(), &mut error as *mut _)
        };

        if !success {
            return Err(Error::from_ffi_error(error));
        }

        Ok(())
    }

    /// Set reference paths by name prefixes
    ///
    /// Selects all paths whose names start with any of the given prefixes.
    ///
    /// # Arguments
    ///
    /// * `prefixes` - Slice of prefix strings (e.g., `&["HG", "NA"]`)
    ///
    /// # Example
    ///
    /// ```no_run
    /// use povu::PovuGraph;
    ///
    /// let mut graph = PovuGraph::load("test.gfa")?;
    /// graph.set_references_from_prefixes(&["HG002", "HG003"])?;
    /// # Ok::<(), povu::Error>(())
    /// ```
    pub fn set_references_from_prefixes(&mut self, prefixes: &[impl AsRef<str>]) -> Result<()> {
        let c_strings: Result<Vec<CString>> = prefixes
            .iter()
            .map(|s| CString::new(s.as_ref()).map_err(Error::from))
            .collect();
        let c_strings = c_strings?;

        let c_ptrs: Vec<*const libc::c_char> = c_strings.iter().map(|s| s.as_ptr()).collect();

        let mut error = ffi::PovuError::default();

        let success = unsafe {
            ffi::povu_graph_set_references_from_prefixes(
                self.inner,
                c_ptrs.as_ptr(),
                c_ptrs.len(),
                &mut error as *mut _,
            )
        };

        if !success {
            return Err(Error::from_ffi_error(error));
        }

        Ok(())
    }

    /// Analyze the graph to find flubbles (regions of variation)
    ///
    /// This runs the flubble detection algorithm and returns an analysis
    /// object that can be used to query the variation structure and generate VCF output.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use povu::PovuGraph;
    ///
    /// let graph = PovuGraph::load("test.gfa")?;
    /// let analysis = graph.analyze()?;
    /// println!("Found {} variation regions", analysis.flubble_count());
    /// # Ok::<(), povu::Error>(())
    /// ```
    pub fn analyze(&self) -> Result<GraphAnalysis> {
        let mut error = ffi::PovuError::default();

        let flubbles = unsafe { ffi::povu_graph_find_flubbles(self.inner, &mut error as *mut _) };

        if flubbles.is_null() {
            return Err(Error::from_ffi_error(error));
        }

        Ok(GraphAnalysis { inner: flubbles })
    }
}

impl Drop for PovuGraph {
    fn drop(&mut self) {
        if !self.inner.is_null() {
            unsafe {
                ffi::povu_graph_free(self.inner);
            }
        }
    }
}

// Safety: PovuGraph can be sent between threads
unsafe impl Send for PovuGraph {}
