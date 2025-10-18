#include "povu_ffi.h"

#include <cstring>
#include <exception>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Povu headers
#include "povu/algorithms/flubbles.hpp"
#include "povu/common/app.hpp"
#include "povu/graph/bidirected.hpp"
#include "povu/graph/pvst.hpp"
#include "povu/graph/spanning_tree.hpp"
#include "povu/io/from_gfa.hpp"
#include "povu/io/to_vcf.hpp"
#include "povu/genomics/vcf.hpp"

// Internal wrapper structures
struct PovuGraph {
    povu::bidirected::VG* vg;
    core::config config;

    PovuGraph(povu::bidirected::VG* graph) : vg(graph) {
        config.set_inc_refs(true);
        config.set_inc_vtx_labels(true);
    }

    ~PovuGraph() {
        delete vg;
    }
};

struct PovuFlubbles {
    povu::pvst::Tree pvst_tree;
    PovuGraph* graph; // non-owning reference

    PovuFlubbles(povu::pvst::Tree&& tree, PovuGraph* g)
        : pvst_tree(std::move(tree)), graph(g) {}
};

struct PovuPvstTree {
    const povu::pvst::Tree* tree; // non-owning reference

    PovuPvstTree(const povu::pvst::Tree* t) : tree(t) {}
};

struct PovuVcfOutput {
    std::string vcf_content;
};

// Helper function to set error
static void set_error(PovuError* error, int code, const char* message) {
    if (error) {
        error->code = code;
        error->message = strdup(message);
    }
}

// Graph loading and management
PovuGraph* povu_graph_new(size_t vertex_capacity, size_t edge_capacity, size_t path_capacity) {
    try {
        povu::bidirected::VG* vg = new povu::bidirected::VG(
            vertex_capacity,
            edge_capacity,
            path_capacity
        );
        return new PovuGraph(vg);
    } catch (const std::exception& e) {
        return nullptr;
    }
}

PovuGraph* povu_graph_from_gfa(const char* gfa_path, PovuError* error) {
    try {
        core::config config;
        config.set_input_gfa(gfa_path);
        config.set_inc_refs(true);
        config.set_inc_vtx_labels(true);

        povu::bidirected::VG* vg = povu::io::from_gfa::to_bd(config);
        if (!vg) {
            set_error(error, 1, "Failed to load GFA file");
            return nullptr;
        }

        return new PovuGraph(vg);
    } catch (const std::exception& e) {
        set_error(error, 1, e.what());
        return nullptr;
    }
}

void povu_graph_free(PovuGraph* graph) {
    delete graph;
}

// Graph building
size_t povu_graph_add_vertex(PovuGraph* graph, uint64_t id, const char* sequence) {
    if (!graph || !sequence) {
        return static_cast<size_t>(-1);
    }

    try {
        std::string seq_str(sequence);
        pt::idx_t idx = graph->vg->add_vertex(id, seq_str);
        return static_cast<size_t>(idx);
    } catch (const std::exception& e) {
        return static_cast<size_t>(-1);
    }
}

size_t povu_graph_add_edge(PovuGraph* graph,
                           uint64_t from_id, PovuOrientation from_orientation,
                           uint64_t to_id, PovuOrientation to_orientation) {
    if (!graph) {
        return static_cast<size_t>(-1);
    }

    try {
        // Convert orientations to v_end_e enum
        povu::types::graph::v_end_e from_end =
            (from_orientation == POVU_ORIENTATION_FORWARD)
                ? povu::types::graph::v_end_e::L
                : povu::types::graph::v_end_e::R;

        povu::types::graph::v_end_e to_end =
            (to_orientation == POVU_ORIENTATION_FORWARD)
                ? povu::types::graph::v_end_e::L
                : povu::types::graph::v_end_e::R;

        pt::idx_t idx = graph->vg->add_edge(from_id, from_end, to_id, to_end);
        return static_cast<size_t>(idx);
    } catch (const std::exception& e) {
        return static_cast<size_t>(-1);
    }
}

bool povu_graph_add_path(PovuGraph* graph, const char* name,
                         const PovuStep* steps, size_t steps_count) {
    if (!graph || !name || !steps) {
        return false;
    }

    try {
        // For now, path addition requires going through the liteseq API
        // or manually constructing the path data structures.
        // This is a TODO - needs more integration with the ref/path system

        // Temporary stub - full implementation requires ref management
        return false;
    } catch (const std::exception& e) {
        return false;
    }
}

void povu_graph_finalize(PovuGraph* graph) {
    if (!graph) {
        return;
    }

    try {
        graph->vg->shrink_to_fit();
        graph->vg->gen_genotype_metadata();
    } catch (const std::exception& e) {
        // Silently fail - not critical
    }
}

// Graph topology queries
size_t povu_graph_vertex_count(const PovuGraph* graph) {
    return graph ? graph->vg->vtx_count() : 0;
}

size_t povu_graph_edge_count(const PovuGraph* graph) {
    return graph ? graph->vg->edge_count() : 0;
}

size_t povu_graph_path_count(const PovuGraph* graph) {
    return graph ? graph->vg->get_ref_count() : 0;
}

PovuVertex* povu_graph_get_vertices(const PovuGraph* graph, size_t* count) {
    if (!graph || !count) {
        if (count) *count = 0;
        return nullptr;
    }

    size_t vtx_count = graph->vg->vtx_count();
    *count = vtx_count;

    PovuVertex* vertices = new PovuVertex[vtx_count];

    for (size_t i = 0; i < vtx_count; i++) {
        const auto& v = graph->vg->get_vertex_by_idx(i);
        vertices[i].id = v.id();

        // Allocate and copy sequence
        const std::string& label = v.get_label();
        char* seq = new char[label.size() + 1];
        std::strcpy(seq, label.c_str());
        vertices[i].sequence = seq;
        vertices[i].sequence_len = label.size();
    }

    return vertices;
}

PovuEdge* povu_graph_get_edges(const PovuGraph* graph, size_t* count) {
    if (!graph || !count) {
        if (count) *count = 0;
        return nullptr;
    }

    size_t edge_count = graph->vg->edge_count();
    *count = edge_count;

    PovuEdge* edges = new PovuEdge[edge_count];

    for (size_t i = 0; i < edge_count; i++) {
        const auto& e = graph->vg->get_edge(i);

        edges[i].from_id = graph->vg->v_idx_to_id(e.get_v1_idx());
        edges[i].from_orientation = static_cast<PovuOrientation>(e.get_v1_end());
        edges[i].to_id = graph->vg->v_idx_to_id(e.get_v2_idx());
        edges[i].to_orientation = static_cast<PovuOrientation>(e.get_v2_end());
    }

    return edges;
}

PovuPath* povu_graph_get_paths(const PovuGraph* graph, size_t* count) {
    if (!graph || !count) {
        if (count) *count = 0;
        return nullptr;
    }

    size_t path_count = graph->vg->get_ref_count();
    *count = path_count;

    if (path_count == 0) {
        return nullptr;
    }

    PovuPath* paths = new PovuPath[path_count];

    for (size_t ref_id = 0; ref_id < path_count; ref_id++) {
        const auto* ref_vec = graph->vg->get_ref_vec(ref_id);
        const auto& ref = graph->vg->get_ref_by_id(ref_id);

        // Get path name
        std::string path_name = ref.get_label();
        char* name = new char[path_name.size() + 1];
        std::strcpy(name, path_name.c_str());
        paths[ref_id].name = name;
        paths[ref_id].name_len = path_name.size();

        // Get steps - access through the ref vector from liteseq
        if (ref_vec && ref_vec->len > 0) {
            paths[ref_id].steps_count = ref_vec->len;
            paths[ref_id].steps = new PovuStep[ref_vec->len];

            for (size_t j = 0; j < ref_vec->len; j++) {
                paths[ref_id].steps[j].vertex_id = ref_vec->ids[j];
                paths[ref_id].steps[j].orientation = static_cast<PovuOrientation>(ref_vec->ors[j]);
            }
        } else {
            paths[ref_id].steps_count = 0;
            paths[ref_id].steps = nullptr;
        }
    }

    return paths;
}

void povu_vertices_free(PovuVertex* vertices, size_t count) {
    if (vertices) {
        for (size_t i = 0; i < count; i++) {
            delete[] vertices[i].sequence;
        }
        delete[] vertices;
    }
}

void povu_edges_free(PovuEdge* edges, size_t /*count*/) {
    delete[] edges;
}

void povu_paths_free(PovuPath* paths, size_t count) {
    if (paths) {
        for (size_t i = 0; i < count; i++) {
            delete[] paths[i].name;
            delete[] paths[i].steps;
        }
        delete[] paths;
    }
}

// Reference selection
bool povu_graph_set_references_from_file(PovuGraph* graph, const char* ref_file, PovuError* error) {
    if (!graph || !ref_file) {
        set_error(error, 1, "Invalid arguments");
        return false;
    }

    try {
        graph->config.set_references_txt(ref_file);
        graph->config.set_ref_input_format(core::input_format_e::file_path);
        return true;
    } catch (const std::exception& e) {
        set_error(error, 1, e.what());
        return false;
    }
}

bool povu_graph_set_references_from_prefixes(PovuGraph* graph, const char** prefixes,
                                              size_t count, PovuError* error) {
    if (!graph || !prefixes) {
        set_error(error, 1, "Invalid arguments");
        return false;
    }

    try {
        std::vector<std::string> prefix_vec;
        for (size_t i = 0; i < count; i++) {
            prefix_vec.push_back(prefixes[i]);
        }
        graph->config.set_ref_name_prefixes(std::move(prefix_vec));
        graph->config.set_ref_input_format(core::input_format_e::params);
        return true;
    } catch (const std::exception& e) {
        set_error(error, 1, e.what());
        return false;
    }
}

// Flubble detection
PovuFlubbles* povu_graph_find_flubbles(PovuGraph* graph, PovuError* error) {
    if (!graph) {
        set_error(error, 1, "Invalid graph");
        return nullptr;
    }

    try {
        // Build spanning tree
        povu::spanning_tree::Tree st(*graph->vg, graph->config);

        // Find flubbles
        povu::pvst::Tree pvst_tree = povu::flubbles::find_flubbles(st, graph->config);

        return new PovuFlubbles(std::move(pvst_tree), graph);
    } catch (const std::exception& e) {
        set_error(error, 1, e.what());
        return nullptr;
    }
}

void povu_flubbles_free(PovuFlubbles* flubbles) {
    delete flubbles;
}

size_t povu_flubbles_count(const PovuFlubbles* flubbles) {
    if (!flubbles) return 0;
    return flubbles->pvst_tree.vertex_count();
}

PovuFlubble* povu_flubbles_get(const PovuFlubbles* flubbles, size_t index) {
    if (!flubbles || index >= flubbles->pvst_tree.vertex_count()) {
        return nullptr;
    }

    // TODO: Implement detailed flubble access
    // This requires more complex iteration through the PVST tree
    // For now, return nullptr - this can be expanded later
    return nullptr;
}

void povu_flubble_free(PovuFlubble* flubble) {
    if (flubble) {
        if (flubble->walks) {
            for (size_t i = 0; i < flubble->walks_count; i++) {
                delete[] flubble->walks[i];
            }
            delete[] flubble->walks;
        }
        delete[] flubble->walk_lengths;
        delete flubble;
    }
}

// PVST tree access
PovuPvstTree* povu_flubbles_get_pvst_tree(const PovuFlubbles* flubbles) {
    if (!flubbles) return nullptr;
    return new PovuPvstTree(&flubbles->pvst_tree);
}

void povu_pvst_tree_free(PovuPvstTree* tree) {
    delete tree;
}

size_t povu_pvst_tree_vertex_count(const PovuPvstTree* tree) {
    return tree ? tree->tree->vertex_count() : 0;
}

// VCF output
PovuVcfOutput* povu_flubbles_call_variants(PovuFlubbles* flubbles, PovuError* error) {
    if (!flubbles || !flubbles->graph) {
        set_error(error, 1, "Invalid flubbles or graph");
        return nullptr;
    }

    try {
        // Create VCF output to string stream
        std::ostringstream oss;

        // TODO: Implement actual variant calling
        // This requires integrating the variant calling pipeline
        // For now, create a placeholder

        set_error(error, 1, "VCF generation not yet implemented in FFI layer");
        return nullptr;

    } catch (const std::exception& e) {
        set_error(error, 1, e.what());
        return nullptr;
    }
}

bool povu_vcf_write_to_file(const PovuVcfOutput* vcf, const char* path, PovuError* error) {
    if (!vcf || !path) {
        set_error(error, 1, "Invalid arguments");
        return false;
    }

    try {
        std::ofstream out(path);
        if (!out) {
            set_error(error, 1, "Failed to open output file");
            return false;
        }
        out << vcf->vcf_content;
        return true;
    } catch (const std::exception& e) {
        set_error(error, 1, e.what());
        return false;
    }
}

char* povu_vcf_to_string(const PovuVcfOutput* vcf, size_t* length) {
    if (!vcf || !length) {
        if (length) *length = 0;
        return nullptr;
    }

    *length = vcf->vcf_content.size();
    char* result = new char[*length + 1];
    std::strcpy(result, vcf->vcf_content.c_str());
    return result;
}

void povu_vcf_free(PovuVcfOutput* vcf) {
    delete vcf;
}

void povu_string_free(char* str) {
    delete[] str;
}

// One-shot convenience function
bool povu_gfa_to_vcf(const char* gfa_path, const char* vcf_path,
                     const char* ref_file, PovuError* error) {
    try {
        core::config config;
        config.set_task(core::task_e::gfa2vcf);
        config.set_input_gfa(gfa_path);
        config.set_inc_refs(true);
        config.set_inc_vtx_labels(true);

        if (ref_file) {
            config.set_references_txt(ref_file);
            config.set_ref_input_format(core::input_format_e::file_path);
        }

        // TODO: Implement full gfa2vcf pipeline
        set_error(error, 1, "gfa_to_vcf not yet implemented in FFI layer");
        return false;

    } catch (const std::exception& e) {
        set_error(error, 1, e.what());
        return false;
    }
}

// Error handling
void povu_error_free(PovuError* error) {
    if (error && error->message) {
        free(error->message);
        error->message = nullptr;
        error->code = 0;
    }
}
