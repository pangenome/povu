#include <cassert>
#include <cstddef>
#include <cstdlib>
// #include <exception>
#include <format>
// #include <functional>
#include <fstream>
#include <iostream>
//#include <map>
#include <set>
#include <sstream>
#include <string>
#include <sys/types.h>
// #include <limits>
// #include <thread>
#include <tuple>
#include <vector>

#include "gfakluge.hpp"

#include "../graph/bidirected.hpp"
#include "../graph/graph.hpp"
#include "./io.hpp"
// #include "../common/types.hpp"
// #include "handlegraph/types.hpp"


namespace io::from_gfa {
namespace bd = povu::bidirected;
namespace pg = povu::graph;
namespace pgt = povu::graph_types;
using namespace povu::graph;


std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}


int foo(const char* filename, std::vector<std::size_t>& v_ids, std::vector<std::tuple<std::size_t, pgt::or_t, std::size_t, pgt::or_t>>& edges) {
  std::ifstream file(filename); // Open the file
  if (!file) {
    std::cerr << "Unable to open file" << std::endl;
    return 1;
  }

  std::string line;
  std::vector<std::string> tokens;

  while (std::getline(file, line)) {
    tokens = splitString(line, '\t');
    if (tokens[0] == "S") {

      v_ids.push_back(stoull(tokens[1]));
    }
    else if (tokens[0] == "L") {

      std::size_t src = stoull(tokens[1]);
      // char src_o = tokens[2][0];
      pgt::or_t src_o = tokens[2][0] == '+' ? pgt::or_t::forward : pgt::or_t::reverse;
      std::size_t snk = stoull(tokens[3]);
      pgt::or_t snk_o = tokens[4][0] == '+' ? pgt::or_t::forward : pgt::or_t::reverse;
      edges.push_back(std::make_tuple(src, src_o, snk, snk_o));
    }
  }

  file.close(); // Close the file
  return 0;
}


/**
 * To a variation graph represented as a bidirected graph
 *
 *
 *
 * @param [in] filename The GFA file to read
 * @return A VariationGraph object from the GFA file
 */
povu::graph::Graph to_pv_graph(const char* filename, const core::config& app_config) {
  std::string fn_name { std::format("[povu::io::{}]", __func__) };

  std::vector<std::size_t> v_ids;

  std::vector<std::tuple<std::size_t, pgt::or_t, std::size_t, pgt::or_t>> edges;

  foo(filename, v_ids, edges);

  pg::Graph g(v_ids.size(), edges.size());


  /*
    add nodes
    ---------

  */

  for (std::size_t v_id : v_ids) { g.add_vertex(v_id); }


  /*
    add edges
    ---------

  */

  for (auto [src, src_or, snk, snk_or] : edges) {

    if (src == snk) {
      if (src_or == snk_or) {
        g.add_edge(src, pgt::VertexEnd::l, src, pgt::VertexEnd::r);
      }

      if (src_or != snk_or) {
        std::cerr << std::format("{} invalid self loop\n", fn_name);
      }
      continue;
    }

    auto v1_end = src_or == pgt::or_t::forward ? pgt::v_end::r : pgt::v_end::l;
    auto v2_end = snk_or == pgt::or_t::forward ? pgt::v_end::l : pgt::v_end::r;

    g.add_edge(src, v1_end, snk, v2_end);
  }

  // populate tips
  // -------------

  for (std::size_t v_idx {}; v_idx < g.size(); ++v_idx) {
    const Vertex& v = g.get_vertex_by_idx(v_idx);
    std::size_t v_id = v.id();
    if (v.get_edges_l().empty() && v.get_edges_r().empty()) {
      if (app_config.verbosity() > 3) {
        std::cerr << std::format(" {} WARN isolated node {} \n", fn_name, v.id());
      }

      g.add_tip(v_id, pgt::VertexEnd::l);
    }
    else if (v.get_edges_l().empty()) { g.add_tip(v_id, pgt::v_end::l); }
    else if (v.get_edges_r().empty()) { g.add_tip(v_id, pgt::v_end::r); }
  }

  return g;
}


/**
 * Count the number of lines of each type in a GFA file.
 *
 * @param[in] filename The GFA file to read.
 * @return A map from line type to number of lines of that type.
 */
void gfa_line_counts(const char* filename, std::map<char, uint64_t>& line_counts) {
  int gfa_fd = -1;
  char* gfa_buf = nullptr;
  std::size_t gfa_filesize = gfak::mmap_open(filename, gfa_buf, gfa_fd);
  if (gfa_fd == -1) {
    std::cerr << "Couldn't open GFA file " << filename << "." << std::endl;
    exit(1);
  }

  for (std::size_t i{}; i < gfa_filesize; ++i) {
    if (i == 0 || gfa_buf[i-1] == '\n') { line_counts[gfa_buf[i]]++; }
    //if (gfa_buf[i] == 'S' || gfa_buf[i] == 'L' || gfa_buf[i] == 'P') {
    //  line_counts[gfa_buf[i]]++;
    //}
  }

  //std::size_t i {};
  //while (i < gfa_filesize) {
  //  if (i == 0 || gfa_buf[i-1] == '\n') { line_counts[gfa_buf[i]]++; }
  //  ++i;
  //}
  gfak::mmap_close(gfa_buf, gfa_fd, gfa_filesize);
  //return counts;
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
void handle_self_loop(bd::VG &vg, id_t v_id, bool src_f, bool snk_f) {
  if (src_f == snk_f) {
    vg.add_edge(v_id, pgt::v_end::l, v_id, pgt::v_end::r);
  }

  if (src_f != snk_f) { std::cerr << __func__ << " unhandled case\n"; }
}


/**
 * To a variation graph represented as a bidirected graph
 *
 *
 *
 * @param [in] filename The GFA file to read
 * @return A VariationGraph object from the GFA file
 */
bd::VG to_bd(const char* filename, const core::config& app_config) {
  std::string fn_name { std::format("[povu::io::{}]", __func__) };

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
  gfa_line_counts(filename, line_counts);

  //std::cout << "min_id: " << min_id << " max_id: " << max_id << std::endl;

  std::size_t node_count = line_counts['S'];
  std::size_t edge_count = line_counts['L'];
  std::size_t path_count = line_counts['P'];

  /*
    Build the digraph
    -----------------

  */
  // compute max nodes by difference between max and min ids
  //std::size_t max_nodes = max_id - min_id + 1;

  // the node count and max and min ids are not necessarily the same
  // this is based on ids and node count not being the same
  bd::VG vg(node_count, edge_count, path_count);

  // we want to start counting from 0 so we have to have an offset_value which would subtruct from the min_id
  // this is because the input graph is not guaranteed to start from 0
  // so we have to make sure that the graph starts from 0 as is the digraph

  //std::size_t offset_value = min_id;
  //std::cout << "offset_value: " << offset_value << std::endl;

  /*
    add nodes
    ---------

  */
  //std::size_t min_id{}; // offset value
  std::size_t curr_id{};
  {
    gg.for_each_sequence_line_in_file(
      filename,
      [&](const gfak::sequence_elem& s) {
        //curr_id = stoull(s.name);
        vg.create_handle(s.sequence, stoull(s.name));
        //min_id = std::min(min_id, curr_id);
      });
  }

  assert(vg.size() == node_count);

  //std::cout << "[io::gfa_to_vg]" << "Nodes added Graph size: " << vg.size() << std::endl;

  /*
    add edges
    ---------

  */
  {
    gg.for_each_edge_line_in_file(filename, [&](const gfak::edge_elem& e) {
      if (e.source_name.empty()) return;

      if (e.source_name == e.sink_name) {
          handle_self_loop(
            vg, stoll(e.source_name), e.source_orientation_forward, e.sink_orientation_forward);
          return;
      }

      /*
        The handlegraph create_edge method doesn't make sense in this case
        because the edge is bidirected
        and not the vertex itself so we don't return a handle to a vertex side
        but rather a handle to a vertex
        this would make more sense in the biedged graph
      */
      auto v1_end = e.source_orientation_forward ? pgt::v_end::r : pgt::v_end::l;
      auto v2_end = e.sink_orientation_forward ? pgt::v_end::l : pgt::v_end::r;

      //if (e.source_name == "94679" ) {
      //std::cerr << "source " << e.source_name << v1_end << " (" << e.source_orientation_forward
      //          << " sink " << e.sink_name << v2_end << " (" << e.sink_orientation_forward << "\n";
      //}


      vg.add_edge(stoll(e.source_name), v1_end, stoll(e.sink_name), v2_end);
    });
  }

  //std::cout  << "[io::gfa_to_vg]" << "Edges added\n";


  // TODO: use handles?
  //std::map<std::string, std::size_t> path_id_map;

  //std::vector<std::vector<std::pair<std::size_t, std::size_t>>> path_spans;


  auto path_elem_idx_to_id_o = [&vg](const gfak::path_elem& path, std::size_t i) {
    const std::string& s = path.segment_names[i];
    std::size_t v_idx = vg.id_to_idx(std::stoull(s));
    bool orientation = path.orientations[i];
    pgt::orientation_t o = orientation ? pgt::orientation_t::forward : pgt::orientation_t::reverse;

    pgt::id_n_orientation_t id_n_orientation = pgt::id_n_orientation_t{v_idx, o};
    return id_n_orientation;
  };

  std::size_t path_pos {}; // the position of a base in a reference path

  /*
    add paths
    ---------

  */
  // do this by associating each node && edge with a reference/color
  if (app_config.get_task() == core::task_t::call && path_count > 0) {
    std::vector<std::vector<pgt::id_n_orientation_t>> raw_paths;
    std::vector<pgt::id_n_orientation_t> raw_path;

    // for each reference path (P line) in the GFA file
    gg.for_each_path_line_in_file(filename, [&](const gfak::path_elem& path) {
      handlegraph::path_handle_t p_h =
        vg.create_path_handle(path.name, *std::begin(path.segment_names) == *std::rbegin(path.segment_names));

      path_pos = 1;

      std::size_t s_v_idx = vg.id_to_idx(std::stoull(*std::begin(path.segment_names)));
      pgt::side_n_id_t path_start = pgt::side_n_id_t{ path.orientations.front() ? pgt::v_end::l : pgt::v_end::r, s_v_idx};

      std::size_t e_v_idx = vg.id_to_idx(std::stoull(*std::rbegin(path.segment_names)));
      pgt::side_n_id_t path_end = pgt::side_n_id_t { path.orientations.back() ? pgt::v_end::r : pgt::v_end::r, e_v_idx };

      // do we need this?
      vg.add_haplotype_start_node(path_start);
      vg.add_haplotype_stop_node(path_end);

      for (std::size_t i{}; i < path.segment_names.size(); ++i) {

        pgt::id_n_orientation_t id_n_orientation = path_elem_idx_to_id_o(path, i);

        raw_path.push_back(id_n_orientation);

        /*
          color the edge
          ...............
        */
        if (i+1 < path.segment_names.size()) {
          pgt::id_n_orientation_t id_n_orientation_next = path_elem_idx_to_id_o(path, i+1);
          bd::Edge &e = vg.get_edge_mut(id_n_orientation, id_n_orientation_next);
          e.add_ref(std::stoll(p_h.data));
        }

        /*
          color the vertex
          ................
         */

        /*
          this can be done through handleGraph's append_step but this is
          preferable in my because it also sets the step value which is
          usable for variant calling
        */

        bd::Vertex& v = vg.get_vertex_mut(id_n_orientation.v_idx);
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
    const bd::Vertex &v = vg.get_vertex(v_idx);
    if (v.get_edges_l().empty() && v.get_edges_r().empty()) {
      std::cerr << std::format(" {} WARN isolated node {} \n", fn_name, v.get_name());
      vg.add_tip(v_idx, pgt::VertexEnd::l);
    }
    else if (v.get_edges_l().empty()) { vg.add_tip(v_idx, pgt::v_end::l); }
    else if (v.get_edges_l().empty()) { vg.add_tip(v_idx, pgt::v_end::r); }
  }


  return vg;
}

};
