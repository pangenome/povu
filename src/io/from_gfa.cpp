#include <cassert>
#include <cstddef>
#include <format>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <sys/types.h>
#include <string>
#include <thread>
#include <vector>

#include <utility>
#include <mutex>
#include <atomic>


#include "gfakluge.hpp"

#include "../graph/bidirected.hpp"
#include "./io.hpp"
#include "handlegraph/types.hpp"

namespace hg = handlegraph;

namespace io::from_gfa {
/**
 * Count the number of lines of each type in a GFA file.
 *
 * @param[in] filename The GFA file to read.
 * @return A map from line type to number of lines of that type.
 */
std::map<char, uint64_t> gfa_line_counts(const char* filename) {
  int gfa_fd = -1;
  char* gfa_buf = nullptr;
  std::size_t gfa_filesize = gfak::mmap_open(filename, gfa_buf, gfa_fd);
  if (gfa_fd == -1) {
    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
    exit(1);
  }
  std::string line;
  std::size_t i = 0;
  //bool seen_newline = true;
  std::map<char, uint64_t> counts;
  while (i < gfa_filesize) {
    if (i == 0 || gfa_buf[i-1] == '\n') { counts[gfa_buf[i]]++; }
    ++i;
  }
  gfak::mmap_close(gfa_buf, gfa_fd, gfa_filesize);
  return counts;
}


/**
 * This fn assumes source (src) and sink (snk) are the same value so no need to
 * pass it twice or check.
 * A self loop can be in the forward, reverse, or mixed strand
 *
 * Examples:
 * L 1 + 1 + is a forward self loop
 * L 1 - 1 - is a reverse self loop
 *
 * L 1 + 1 -  and L 1 - 1 + are mixed self loops that aren't
 * representable in a bidirected graph without node duplication
 *
 * @param[in] v_id vertex index
 * @param[in] src_f true if source orientation is in forward strand
 * @param[in] snk_f true if sink orientation is in forward strand
 */
void handle_self_loop(bidirected::VariationGraph &vg, id_t v_id, bool src_f, bool snk_f) {

  if (src_f == snk_f) {
    vg.add_edge(v_id, bidirected::VertexEnd::l,
                v_id, bidirected::VertexEnd::r);
  }

  if (src_f != snk_f) {
    std::cerr << __func__ << " unhandled case\n";
  }

}


/**
 * To a variation graph represented as a bidirected graph
 *
 *
 *
 * @param [in] filename The GFA file to read
 * @return A VariationGraph object from the GFA file
 */
bidirected::VariationGraph to_vg(const char* filename, const core::config& app_config) {
  std::string fn_name = std::format("[povu::io::{}]", __func__);

  //#ifdef DEBUG
  //if (app_config.verbosity() > 2) { std::cout << fn_name << std::endl; }
  //#endif

  gfak::GFAKluge gg = gfak::GFAKluge();

  gg.parse_gfa_file(filename);

  /*
    Preprocess the GFA
    ------------------

    scan over the file to count edges and sequences in parallel
  */
  std::map<char, uint64_t> line_counts;
  std::size_t min_id = std::numeric_limits<uint64_t>::max();
  std::size_t max_id = std::numeric_limits<uint64_t>::min();
  {
    std::thread x([&]() {
      gg.for_each_sequence_line_in_file(
        filename,
        [&](gfak::sequence_elem s) {
          uint64_t id = stol(s.name);
          min_id = std::min(min_id, id);
          max_id = std::max(max_id, id);
        });
    });
    line_counts = gfa_line_counts(filename);
    x.join();
  }

  //std::cout << "min_id: " << min_id << " max_id: " << max_id << std::endl;

  std::size_t node_count = line_counts['S'];
  std::size_t edge_count = line_counts['L'];
  std::size_t path_count = line_counts['P'];

  /*
    Build the digraph
    -----------------

  */
  // compute max nodes by difference between max and min ids
  std::size_t max_nodes = max_id - min_id + 1;

  // the node count and max and min ids are not necessarily the same
  // this is based on ids and node count not being the same
  bidirected::VariationGraph vg(max_nodes, edge_count, path_count);

  // we want to start counting from 0 so we have to have an offset_value which would subtruct from the min_id
  // this is because the input graph is not guaranteed to start from 0
  // so we have to make sure that the graph starts from 0 as is the digraph

  std::size_t offset_value = min_id;
  //std::cout << "offset_value: " << offset_value << std::endl;

  // add nodes
  // ---------
  {
    gg.for_each_sequence_line_in_file(
      filename,
      [&](const gfak::sequence_elem& s) {
        vg.create_handle(s.sequence, std::stoll(s.name));
        //vg.get_vertex_mut(std::stoll(s.name) - offset_value).set_name(s.name);
      });
  }

  assert(vg.size() == node_count);

  //std::cout << "[io::gfa_to_vg]" << "Nodes added Graph size: " << vg.size() << std::endl;

  // add edges
  // ---------
  {
    gg.for_each_edge_line_in_file(
      filename,
      [&](const gfak::edge_elem& e) {
        if (e.source_name.empty()) return;

        if (e.source_name == e.sink_name) {
          handle_self_loop(vg, stoll(e.source_name),
                           e.source_orientation_forward, e.sink_orientation_forward);
          return;
        }

        /*
          The handlegraph create_edge method doesn't make sense in this case
          because the edge is bidirected
          and not the vertex itself so we don't return a handle to a vertex side
          but rather a handle to a vertex
          this would make more sense in the biedged graph
        */
        auto v1_end =
          e.source_orientation_forward ? bidirected::VertexEnd::r : bidirected::VertexEnd::l;
        auto v2_end =
          e.sink_orientation_forward ? bidirected::VertexEnd::l : bidirected::VertexEnd::r;

        //std::cout << "source "<< e.source_name << v1_end << " sink " << e.sink_name << v2_end << "\n";

        vg.add_edge(stoll(e.source_name), v1_end,
                    stoll(e.sink_name), v2_end);
      });
  }

  //std::cout  << "[io::gfa_to_vg]" << "Edges added\n";


  // TODO: use handles?
  //std::map<std::string, std::size_t> path_id_map;

  //std::vector<std::vector<std::pair<std::size_t, std::size_t>>> path_spans;

  std::size_t path_pos {};

  // add paths
  // ---------
  // do this by associating each node with a path
  if (path_count > 0) {
    std::vector<std::vector<bidirected::id_n_orientation_t>> raw_paths;
    std::vector<bidirected::id_n_orientation_t> raw_path;

    gg.for_each_path_line_in_file(
        filename, [&](const gfak::path_elem &path) {
          handlegraph::path_handle_t p_h = vg.create_path_handle(
              path.name,
              *std::begin(path.segment_names) == *std::rbegin(path.segment_names));

          path_pos = 0;

          bidirected::side_n_id_t path_start =
            bidirected::side_n_id_t{path.orientations.front() ? bidirected::VertexEnd::l : bidirected::VertexEnd::r,
                                    std::stoull(*std::begin(path.segment_names)) - offset_value};
          bidirected::side_n_id_t path_end =
            bidirected::side_n_id_t{path.orientations.back() ? bidirected::VertexEnd::r : bidirected::VertexEnd::r,
                                    std::stoull(*std::rbegin(path.segment_names)) - offset_value};

          vg.add_haplotype_start_node(path_start);
          vg.add_haplotype_stop_node(path_end);

          for (std::size_t i{}; i < path.segment_names.size(); ++i) {
            const std::string &s = path.segment_names[i];
            bool orientation = path.orientations[i];

            handlegraph::nid_t id = std::stoull(s) - offset_value;

            bidirected::id_n_orientation_t id_n_orientation =
              bidirected::id_n_orientation_t{
              std::stoull(s) - offset_value,
              orientation ? bidirected::orientation_t::forward : bidirected::orientation_t::reverse};

            raw_path.push_back(id_n_orientation);

            /*
              this can be done through handleGraph's append_step but this is
              preferable in my because it also sets the step value which is
              usable for variant calling
            */

            bidirected::Vertex &v = vg.get_vertex_mut(id);
            v.add_path(std::stoll(p_h.data), path_pos);

            path_pos += v.get_label().length();
          }
          raw_paths.push_back(raw_path);
          raw_path.clear();
        });

     vg.set_raw_paths(raw_paths);
  }

  // populate tips
  // -------------

  for (std::size_t v_idx{}; v_idx < vg.size(); ++v_idx) {
    const bidirected::Vertex &v = vg.get_vertex(v_idx);
    if (v.get_edges_l().empty() && v.get_edges_r().empty()) {
      std::cerr << std::format(" {} WARN isolated node {} \n", fn_name, v.get_name());
      vg.add_tip(v_idx, graph_types::VertexEnd::l);
    }
    else if (v.get_edges_l().empty()) { vg.add_tip(v_idx, graph_types::VertexEnd::l); }
    else if (v.get_edges_l().empty()) { vg.add_tip(v_idx, graph_types::VertexEnd::r); }
  }

  //std::cout << "Paths added " << std::endl;
  return vg;
}

};
