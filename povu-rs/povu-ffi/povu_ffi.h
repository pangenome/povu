/**
 * @file povu_ffi.h
 * @brief C FFI interface for Povu pangenome variation toolkit
 *
 * This header provides a C-compatible interface to the Povu C++ library,
 * primarily intended for Rust FFI bindings but usable from any C/C++ code.
 *
 * @author Povu Contributors
 */

#ifndef POVU_FFI_H
#define POVU_FFI_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * =============================================================================
 * Opaque Types
 * =============================================================================
 */

/** @brief Opaque handle to a pangenome variation graph */
typedef struct PovuGraph PovuGraph;

/** @brief Opaque handle to flubble (variation region) detection results */
typedef struct PovuFlubbles PovuFlubbles;

/** @brief Opaque handle to a Pangenome Variation Structure Tree */
typedef struct PovuPvstTree PovuPvstTree;

/** @brief Opaque handle to VCF output data */
typedef struct PovuVcfOutput PovuVcfOutput;

/*
 * =============================================================================
 * Enums and Data Structures
 * =============================================================================
 */

/**
 * @brief Orientation of a vertex in a bidirected graph
 *
 * In a bidirected graph, each vertex can be traversed in two orientations:
 * - FORWARD: Left-to-right / 5' to 3' direction
 * - REVERSE: Right-to-left / 3' to 5' direction (reverse complement)
 */
typedef enum {
    POVU_ORIENTATION_FORWARD = 0,  /**< Forward orientation (left/5') */
    POVU_ORIENTATION_REVERSE = 1   /**< Reverse orientation (right/3') */
} PovuOrientation;

/**
 * @brief Vertex (segment) in the pangenome graph
 *
 * Represents a node in the graph with its ID and DNA sequence.
 */
typedef struct {
    uint64_t id;           /**< Unique vertex identifier */
    const char* sequence;  /**< DNA sequence (null-terminated) */
    size_t sequence_len;   /**< Length of the sequence in base pairs */
} PovuVertex;

/**
 * @brief Edge connecting two vertices in the bidirected graph
 *
 * Each edge connects two vertices at specific orientations.
 */
typedef struct {
    uint64_t from_id;                /**< Source vertex ID */
    PovuOrientation from_orientation; /**< Source vertex orientation */
    uint64_t to_id;                   /**< Target vertex ID */
    PovuOrientation to_orientation;   /**< Target vertex orientation */
} PovuEdge;

/**
 * @brief A single step in a path through the graph
 *
 * Represents visiting a specific vertex in a specific orientation.
 */
typedef struct {
    uint64_t vertex_id;         /**< Vertex being visited */
    PovuOrientation orientation; /**< Orientation of traversal */
} PovuStep;

/**
 * @brief A reference path through the pangenome graph
 *
 * Represents a complete walk through the graph, typically corresponding
 * to a reference genome or haplotype.
 */
typedef struct {
    const char* name;      /**< Path name (e.g., "sample#1#chr1") */
    size_t name_len;       /**< Length of the name string */
    PovuStep* steps;       /**< Array of steps in the path */
    size_t steps_count;    /**< Number of steps in the path */
} PovuPath;

/**
 * @brief A flubble (region of variation/bubble)
 *
 * Represents a detected region of variation with multiple possible
 * paths through the graph.
 */
typedef struct {
    uint64_t id;                  /**< Flubble identifier */
    const char* type_name;        /**< Type: "flubble", "tiny", "parallel", "concealed", etc. */
    uint64_t start_vertex_id;     /**< Starting vertex of the region */
    uint64_t end_vertex_id;       /**< Ending vertex of the region */
    PovuStep** walks;             /**< Array of walks (each walk is array of steps) */
    size_t* walk_lengths;         /**< Length of each walk */
    size_t walks_count;           /**< Number of alternative walks */
} PovuFlubble;

/**
 * @brief Error information from FFI calls
 *
 * When an FFI function fails, this structure contains the error details.
 * The caller must free the message using povu_error_free().
 */
typedef struct {
    int code;         /**< Error code (non-zero on error) */
    char* message;    /**< Error message (must be freed) */
} PovuError;

/*
 * =============================================================================
 * Graph Loading and Management
 * =============================================================================
 */

/**
 * @brief Load a pangenome graph from a GFA file
 *
 * @param gfa_path Path to the GFA file to load
 * @param error Error structure to receive error information (can be NULL)
 * @return PovuGraph* Graph handle, or NULL on error
 *
 * @note The returned graph must be freed using povu_graph_free()
 */
PovuGraph* povu_graph_from_gfa(const char* gfa_path, PovuError* error);

/**
 * @brief Free a graph and all associated resources
 *
 * @param graph Graph to free (can be NULL)
 */
void povu_graph_free(PovuGraph* graph);

/*
 * =============================================================================
 * Graph Topology Queries
 * =============================================================================
 */

/** @brief Get the number of vertices in the graph */
size_t povu_graph_vertex_count(const PovuGraph* graph);

/** @brief Get the number of edges in the graph */
size_t povu_graph_edge_count(const PovuGraph* graph);

/** @brief Get the number of reference paths in the graph */
size_t povu_graph_path_count(const PovuGraph* graph);

/**
 * @brief Get all vertices from the graph
 *
 * @param graph Graph to query
 * @param count Output: number of vertices returned
 * @return Array of vertices (must be freed with povu_vertices_free)
 */
PovuVertex* povu_graph_get_vertices(const PovuGraph* graph, size_t* count);

/**
 * @brief Get all edges from the graph
 *
 * @param graph Graph to query
 * @param count Output: number of edges returned
 * @return Array of edges (must be freed with povu_edges_free)
 */
PovuEdge* povu_graph_get_edges(const PovuGraph* graph, size_t* count);

/**
 * @brief Get all reference paths from the graph
 *
 * @param graph Graph to query
 * @param count Output: number of paths returned
 * @return Array of paths (must be freed with povu_paths_free)
 */
PovuPath* povu_graph_get_paths(const PovuGraph* graph, size_t* count);

/** @brief Free vertex array returned by povu_graph_get_vertices */
void povu_vertices_free(PovuVertex* vertices, size_t count);

/** @brief Free edge array returned by povu_graph_get_edges */
void povu_edges_free(PovuEdge* edges, size_t count);

/** @brief Free path array returned by povu_graph_get_paths */
void povu_paths_free(PovuPath* paths, size_t count);

/*
 * =============================================================================
 * Reference Selection
 * =============================================================================
 */

/**
 * @brief Set reference paths from a file
 *
 * @param graph Graph to configure
 * @param ref_file Path to file containing reference path names (one per line)
 * @param error Error structure (can be NULL)
 * @return true on success, false on error
 */
bool povu_graph_set_references_from_file(PovuGraph* graph, const char* ref_file, PovuError* error);

/**
 * @brief Set reference paths by name prefixes
 *
 * Selects all paths whose names start with any of the given prefixes.
 *
 * @param graph Graph to configure
 * @param prefixes Array of prefix strings
 * @param count Number of prefixes in the array
 * @param error Error structure (can be NULL)
 * @return true on success, false on error
 */
bool povu_graph_set_references_from_prefixes(PovuGraph* graph, const char** prefixes, size_t count, PovuError* error);

/*
 * =============================================================================
 * Flubble Detection
 * =============================================================================
 */

/**
 * @brief Detect flubbles (regions of variation) in the graph
 *
 * @param graph Graph to analyze
 * @param error Error structure (can be NULL)
 * @return Flubbles handle, or NULL on error (must be freed with povu_flubbles_free)
 */
PovuFlubbles* povu_graph_find_flubbles(PovuGraph* graph, PovuError* error);

/** @brief Free flubbles detection results */
void povu_flubbles_free(PovuFlubbles* flubbles);

/** @brief Get the number of flubbles detected */
size_t povu_flubbles_count(const PovuFlubbles* flubbles);

/**
 * @brief Get a specific flubble by index
 *
 * @param flubbles Flubbles handle
 * @param index Index of the flubble to retrieve
 * @return Flubble structure (must be freed with povu_flubble_free), or NULL
 *
 * @note This function is not yet fully implemented
 */
PovuFlubble* povu_flubbles_get(const PovuFlubbles* flubbles, size_t index);

/** @brief Free a single flubble structure */
void povu_flubble_free(PovuFlubble* flubble);

/*
 * =============================================================================
 * PVST Tree Access
 * =============================================================================
 */

/**
 * @brief Get the PVST (Pangenome Variation Structure Tree) from flubbles
 *
 * The PVST represents the hierarchical nesting relationships between
 * variation regions.
 *
 * @param flubbles Flubbles handle
 * @return PVST tree handle (must be freed with povu_pvst_tree_free)
 */
PovuPvstTree* povu_flubbles_get_pvst_tree(const PovuFlubbles* flubbles);

/** @brief Free PVST tree */
void povu_pvst_tree_free(PovuPvstTree* tree);

/** @brief Get the number of vertices in the PVST tree */
size_t povu_pvst_tree_vertex_count(const PovuPvstTree* tree);

/*
 * =============================================================================
 * VCF Output
 * =============================================================================
 */

/**
 * @brief Call variants from detected flubbles
 *
 * @param flubbles Flubbles to call variants from
 * @param error Error structure (can be NULL)
 * @return VCF output handle, or NULL on error (must be freed with povu_vcf_free)
 *
 * @note This function is not yet fully implemented
 */
PovuVcfOutput* povu_flubbles_call_variants(PovuFlubbles* flubbles, PovuError* error);

/**
 * @brief Write VCF output to a file
 *
 * @param vcf VCF output handle
 * @param path Output file path
 * @param error Error structure (can be NULL)
 * @return true on success, false on error
 */
bool povu_vcf_write_to_file(const PovuVcfOutput* vcf, const char* path, PovuError* error);

/**
 * @brief Get VCF output as a string
 *
 * @param vcf VCF output handle
 * @param length Output: length of the returned string
 * @return VCF string (must be freed with povu_string_free), or NULL on error
 */
char* povu_vcf_to_string(const PovuVcfOutput* vcf, size_t* length);

/** @brief Free VCF output handle */
void povu_vcf_free(PovuVcfOutput* vcf);

/** @brief Free a string returned by this API */
void povu_string_free(char* str);

/*
 * =============================================================================
 * Convenience Functions
 * =============================================================================
 */

/**
 * @brief One-shot GFA to VCF conversion
 *
 * Loads a GFA file, detects flubbles, and writes VCF output in a single call.
 *
 * @param gfa_path Path to input GFA file
 * @param vcf_path Path to output VCF file
 * @param ref_file Optional path to reference paths file (can be NULL)
 * @param error Error structure (can be NULL)
 * @return true on success, false on error
 *
 * @note This function is not yet fully implemented
 */
bool povu_gfa_to_vcf(const char* gfa_path, const char* vcf_path,
                     const char* ref_file, PovuError* error);

/*
 * =============================================================================
 * Error Handling
 * =============================================================================
 */

/**
 * @brief Free an error structure
 *
 * Frees the message string in the error structure. Safe to call on
 * already-freed or NULL errors.
 *
 * @param error Error structure to free (can be NULL)
 */
void povu_error_free(PovuError* error);

#ifdef __cplusplus
}
#endif

#endif // POVU_FFI_H
