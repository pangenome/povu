#include <format>
#include <sys/types.h>
#include <utility>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <thread>
#include <mutex>
#include <functional>
#include <atomic>
#include <vector>
#include <functional>

#include "gfakluge.hpp"

#include "../graph/digraph.hpp"
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
 * Read a GFA file into a DiGraph.
 *
 *
 * @param[in] filename The GFA file to read.
 * @param[out] dg The DiGraph to read into.
 */
digraph::DiGraph to_digraph(const char* filename) {

  std::cout << "[io::gfa_to_digraph]" << "\n";

  gfak::GFAKluge gg = gfak::GFAKluge();

  std::cout << "he\n";
  gg.parse_gfa_file(filename);
  std::cout << "oe\n";
  /*
    Preprocess the GFA
    ------------------

scan over the file to count edges and sequences in parallel
  */


  std::map<char, uint64_t> line_counts;
  std::size_t min_id = std::numeric_limits<uint64_t>::max();
  std::size_t max_id = std::numeric_limits<uint64_t>::min();
  {
std::thread x(
  [&]() {
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

  std::cout << "min_id: " << min_id << " max_id: " << max_id << std::endl;

  std::size_t node_count = line_counts['S'];
  std::size_t edge_count = line_counts['L'];
  std::size_t path_count = line_counts['P'];


  /*
Build the digraph
-----------------

  */

  // compute max nodes by difference between max and min ids
  std::size_t max_nodes = max_id - min_id + 1;

  digraph::DiGraph dg(max_nodes, path_count);

  // we want to start counting from 0 so we have to have an offset_value which would subtruct from the min_id
  // this is because the input graph is not guaranteed to start from 0
  // so we have to make sure that the graph starts from 0 as is the digraph

  std::size_t offset_value = min_id;
  std::cout << "offset_value: " << offset_value << std::endl;

  // add nodes
  // ---------
  {
gg.for_each_sequence_line_in_file(
  filename,
  [&](const gfak::sequence_elem& s) {
dg.create_handle(s.sequence, std::stoll(s.name) - offset_value);
  });
  }

  assert(dg.size() == node_count);

  std::cout << "Nodes added "
<< " Graph size: " << dg.size() << std::endl;

  // add edges
  // ---------
  {
gg.for_each_edge_line_in_file(
  filename,
  [&](const gfak::edge_elem& e) {
if (e.source_name.empty()) return;
// TODO: is this not pointless computation for the sake of following the libhandlegraph the API?
hg::handle_t a =
  dg.get_handle(stoll(e.source_name) - offset_value,
!e.source_orientation_forward);
hg::handle_t b =
  dg.get_handle(stoll(e.sink_name) - offset_value,
!e.sink_orientation_forward);
dg.create_edge(a, b);
  });
  }

  std::cout << "Edges added "
<< " Graph size: " << dg.size() << std::endl;


  // TODO: use handles?
  //std::map<std::string, std::size_t> path_id_map;

  //std::vector<std::vector<std::pair<std::size_t, std::size_t>>> path_spans;

  std::size_t path_pos{0};

  // add paths
  // ---------
  // do this by associating each node with a path
  if (path_count > 0) {
gg.for_each_path_line_in_file(
  filename,
  [&](const gfak::path_elem& path) {
handlegraph::path_handle_t p_h = dg.create_path_handle(path.name);

//metadata.path_id_map[path.name] = path_counter++;

//auto x = std::vector<std::pair<std::size_t, std::size_t>>();
//x.reserve(path.segment_names.size());

//metadata.path_spans.push_back(x);

//std::vector<std::pair<std::size_t, std::size_t>>& curr_spans =
//metadata.path_spans.back();

path_pos = 0;

for (auto& s : path.segment_names) {

  handlegraph::nid_t id = std::stoull(s) - offset_value;

  // TODO: go through the handle graph API
  // this can be done through handle but this is preferable in my case
  digraph::Vertex& v = dg.get_vertex_mut(id);

  //if (curr_spans.empty()) {
  //curr_spans.push_back(std::make_pair(0, v.get_seq().length()));
  //}
  //else {
  //std::pair<std::size_t, std::size_t> prev =
  //  curr_spans.back();
  //std::pair<std::size_t, std::size_t> curr =
  //  std::make_pair(prev.second + prev.first, v.get_seq().length());

  //curr_spans.push_back(curr);
  //}

  if (v.add_path(std::stoll(p_h.data), path_pos) < 0) {
std::cout << "error setting path" << std::endl;
  }

  path_pos += v.get_seq().length();
  //std::cout << "path_pos: " << path_pos << std::endl;
}

  });
  }

  std::cout << "Paths added "
<< " Graph size: " << dg.size() << std::endl;

  dg.compute_start_nodes();
  dg.compute_stop_nodes();

  return dg;
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

  #ifdef DEBUG
  if (app_config.verbosity() > 2) { std::cout << fn_name << std::endl; }
  #endif

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
        vg.create_handle(s.sequence, std::stoll(s.name) - offset_value);
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
          handle_self_loop(vg, stoll(e.source_name) - offset_value,
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

        vg.add_edge(stoll(e.source_name) - offset_value, v1_end,
                    stoll(e.sink_name) - offset_value, v2_end);
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
    gg.for_each_path_line_in_file(
      filename,
      [&](const gfak::path_elem& path) {
        handlegraph::path_handle_t p_h =
          vg.create_path_handle(path.name,
                                *std::begin(path.segment_names) == *std::rbegin(path.segment_names));
        path_pos = 0;

        handlegraph::nid_t start_id = std::stoull(*std::begin(path.segment_names)) - offset_value;
        handlegraph::nid_t end_id = std::stoull(*std::rbegin(path.segment_names)) - offset_value;

        vg.add_haplotype_start_node(start_id);
        vg.add_haplotype_stop_node(end_id);

        for (auto& s : path.segment_names) {
          handlegraph::nid_t id = std::stoull(s) - offset_value;

          /*
            this can be done through handleGraph's append_step but this is preferable in my
            because it also sets the step value
            which is usable for variant calling
          */

          bidirected::Vertex& v = vg.get_vertex_mut(id);
          v.add_path(std::stoll(p_h.data), path_pos);

          path_pos += v.get_label().length();
        }
      });
  }

  //std::cout << "Paths added " << std::endl;
  return vg;
}

};
