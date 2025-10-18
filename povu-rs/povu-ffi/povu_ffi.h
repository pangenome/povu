#ifndef POVU_FFI_H
#define POVU_FFI_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Opaque types
typedef struct PovuGraph PovuGraph;
typedef struct PovuFlubbles PovuFlubbles;
typedef struct PovuPvstTree PovuPvstTree;
typedef struct PovuVcfOutput PovuVcfOutput;

// Graph orientation
typedef enum {
    POVU_ORIENTATION_FORWARD = 0,
    POVU_ORIENTATION_REVERSE = 1
} PovuOrientation;

// Vertex/segment information
typedef struct {
    uint64_t id;
    const char* sequence;
    size_t sequence_len;
} PovuVertex;

// Edge information
typedef struct {
    uint64_t from_id;
    PovuOrientation from_orientation;
    uint64_t to_id;
    PovuOrientation to_orientation;
} PovuEdge;

// Path/walk step
typedef struct {
    uint64_t vertex_id;
    PovuOrientation orientation;
} PovuStep;

// Reference path
typedef struct {
    const char* name;
    size_t name_len;
    PovuStep* steps;
    size_t steps_count;
} PovuPath;

// Flubble/bubble information
typedef struct {
    uint64_t id;
    const char* type_name;  // "flubble", "tiny", "parallel", "concealed", etc.
    uint64_t start_vertex_id;
    uint64_t end_vertex_id;
    PovuStep** walks;       // Array of walks (each walk is array of steps)
    size_t* walk_lengths;   // Length of each walk
    size_t walks_count;
} PovuFlubble;

// Error handling
typedef struct {
    int code;
    char* message;
} PovuError;

// Graph loading and management
PovuGraph* povu_graph_from_gfa(const char* gfa_path, PovuError* error);
void povu_graph_free(PovuGraph* graph);

// Graph topology queries
size_t povu_graph_vertex_count(const PovuGraph* graph);
size_t povu_graph_edge_count(const PovuGraph* graph);
size_t povu_graph_path_count(const PovuGraph* graph);

PovuVertex* povu_graph_get_vertices(const PovuGraph* graph, size_t* count);
PovuEdge* povu_graph_get_edges(const PovuGraph* graph, size_t* count);
PovuPath* povu_graph_get_paths(const PovuGraph* graph, size_t* count);

void povu_vertices_free(PovuVertex* vertices, size_t count);
void povu_edges_free(PovuEdge* edges, size_t count);
void povu_paths_free(PovuPath* paths, size_t count);

// Reference selection
bool povu_graph_set_references_from_file(PovuGraph* graph, const char* ref_file, PovuError* error);
bool povu_graph_set_references_from_prefixes(PovuGraph* graph, const char** prefixes, size_t count, PovuError* error);

// Flubble detection
PovuFlubbles* povu_graph_find_flubbles(PovuGraph* graph, PovuError* error);
void povu_flubbles_free(PovuFlubbles* flubbles);

size_t povu_flubbles_count(const PovuFlubbles* flubbles);
PovuFlubble* povu_flubbles_get(const PovuFlubbles* flubbles, size_t index);
void povu_flubble_free(PovuFlubble* flubble);

// PVST tree access
PovuPvstTree* povu_flubbles_get_pvst_tree(const PovuFlubbles* flubbles);
void povu_pvst_tree_free(PovuPvstTree* tree);

size_t povu_pvst_tree_vertex_count(const PovuPvstTree* tree);
// TODO: More PVST tree traversal functions

// VCF output
PovuVcfOutput* povu_flubbles_call_variants(PovuFlubbles* flubbles, PovuError* error);
bool povu_vcf_write_to_file(const PovuVcfOutput* vcf, const char* path, PovuError* error);
char* povu_vcf_to_string(const PovuVcfOutput* vcf, size_t* length);
void povu_vcf_free(PovuVcfOutput* vcf);
void povu_string_free(char* str);

// One-shot convenience function
bool povu_gfa_to_vcf(const char* gfa_path, const char* vcf_path,
                     const char* ref_file, PovuError* error);

// Error handling
void povu_error_free(PovuError* error);

#ifdef __cplusplus
}
#endif

#endif // POVU_FFI_H
