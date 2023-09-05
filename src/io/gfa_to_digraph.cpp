
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
#include "./io.hpp"

namespace io {

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

// TODO:
// - handle orientation of edges?
void reader() {
  gfak::GFAKluge gg = gfak::GFAKluge();
  std::string filename = "/home/sluggie/src/phd/domibubble-cpp/deps/gfakluge/data/test.gfa";
  gg.parse_gfa_file(filename);
  std::cout << gg << std::endl;

  // count the number of nodes, edges and paths
  std::map<char, uint64_t> line_counts;
  uint64_t min_id = std::numeric_limits<uint64_t>::max();
  uint64_t max_id = std::numeric_limits<uint64_t>::min();
  {
	std::thread x(
	  [&]() {
		gg.for_each_sequence_line_in_file(
		  filename.c_str(),
		  [&](gfak::sequence_elem s) {
			uint64_t id = stol(s.name);
			min_id = std::min(min_id, id);
			max_id = std::max(max_id, id);
		  });
	  });
	line_counts = gfa_line_counts(filename.c_str());
	x.join();
  }

  uint64_t node_count = line_counts['S'];
  uint64_t edge_count = line_counts['L'];
  uint64_t path_count = line_counts['P'];

  digraph::DiGraph graph = digraph::DiGraph(node_count);
  
  gg.for_each_sequence_line_in_file(
	filename.c_str(),
	[&](const gfak::sequence_elem& s) {
	  // s.name is a string but is a number
	  std::cout << s.name << " " << s.sequence << std::endl;
	  //uint64_t id = stol(s.name);
	  graph.create_handle(s.sequence, std::stoll(s.name));
	  //graph->create_handle(s.sequence, id - id_increment);
	  //if (progress) progress_meter->increment(1);
	});

	
  gg.for_each_edge_line_in_file(
	filename.c_str(),
	[&](const gfak::edge_elem& e) {
	  if (e.source_name.empty()) return;
	  std::cout << "src: " << e.source_name << " sink: " << e.sink_name << " "
				<< e.source_orientation_forward << " " << e.sink_orientation_forward
				<< std::endl;

	  // TODO: is this not pointless computation for the sake of following the libhandlegraph the API?
	  // graph.create_edge(e.source_name, e.sink_name);
	  // if (e.source_name.empty()) return;
	  // handlegraph::handle_t a = graph.get_handle(stol(e.source_name), !e.source_orientation_forward);
	  // handlegraph::handle_t b = graph.get_handle(stol(e.sink_name), !e.sink_orientation_forward);
	  handlegraph::handle_t a = graph.get_handle(stoll(e.source_name), !e.source_orientation_forward);
	  handlegraph::handle_t b = graph.get_handle(stoll(e.sink_name), !e.sink_orientation_forward);
	  graph.create_edge(a, b);
	  //if (progress) progress_meter->increment(1);
	});

  graph.compute_start_nodes();
  graph.compute_stop_nodes();

  graph.print_dot();
  
  //Graph sh;

  //std::cout << sh.foo(8,9) << "\n";
  
  //auto hg = Graph();

  //handlegraph::MutableHandleGraph* graph = new handlegraph::MutableHandleGraph();


  //auto hg = handlegraph::MutableHandleGraph();
  //gfa_to_handle( "../gfakluge/data/test.gfa", &hg, 1, true, false);
  
  //return 0;
}

/**
 * Read a GFA file into a DiGraph.
 *
 *
 * @param[in] filename The GFA file to read.
 * @param[out] dg The DiGraph to read into.
 */
void reader(char* filename, digraph::DiGraph dg) {
  gfak::GFAKluge gg = gfak::GFAKluge();
  //std::string filename = "/home/sluggie/src/phd/domibubble-cpp/deps/gfakluge/data/test.gfa";
  gg.parse_gfa_file(filename);
  //std::cout << gg << std::endl;

  // count the number of nodes, edges and paths
  std::map<char, uint64_t> line_counts;
  uint64_t min_id = std::numeric_limits<uint64_t>::max();
  uint64_t max_id = std::numeric_limits<uint64_t>::min();
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

  uint64_t node_count = line_counts['S'];
  uint64_t edge_count = line_counts['L'];
  uint64_t path_count = line_counts['P'];

  //digraph::DiGraph graph = digraph::DiGraph(node_count);
  
  gg.for_each_sequence_line_in_file(
	filename,
	[&](const gfak::sequence_elem& s) {
	  // s.name is a string but is a number
	  std::cout << s.name << " " << s.sequence << std::endl;
	  //uint64_t id = stol(s.name);
	  dg.create_handle(s.sequence, std::stoll(s.name));
	  //graph->create_handle(s.sequence, id - id_increment);
	  //if (progress) progress_meter->increment(1);
	});

	
  gg.for_each_edge_line_in_file(
	filename,
	[&](const gfak::edge_elem& e) {
	  if (e.source_name.empty()) return;
	  std::cout << "src: " << e.source_name << " sink: " << e.sink_name << " "
				<< e.source_orientation_forward << " " << e.sink_orientation_forward
				<< std::endl;

	  // TODO: is this not pointless computation for the sake of following the libhandlegraph the API?
	  // graph.create_edge(e.source_name, e.sink_name);
	  // if (e.source_name.empty()) return;
	  // handlegraph::handle_t a = graph.get_handle(stol(e.source_name), !e.source_orientation_forward);
	  // handlegraph::handle_t b = graph.get_handle(stol(e.sink_name), !e.sink_orientation_forward);
	  handlegraph::handle_t a = dg.get_handle(stoll(e.source_name), !e.source_orientation_forward);
	  handlegraph::handle_t b = dg.get_handle(stoll(e.sink_name), !e.sink_orientation_forward);
	  dg.create_edge(a, b);
	  //if (progress) progress_meter->increment(1);
	});

  dg.compute_start_nodes();
  dg.compute_stop_nodes();

  dg.print_dot();
  
}

  
}
