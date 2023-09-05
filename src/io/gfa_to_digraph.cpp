
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

// TODO:
// - what about nodes without incoming or outgoing edges?
/**
 * Read a GFA file into a DiGraph.
 *
 *
 * @param[in] filename The GFA file to read.
 * @param[out] dg The DiGraph to read into.
 */
void gfa_to_digraph(char* filename, digraph::DiGraph* dg) {
  gfak::GFAKluge gg = gfak::GFAKluge();
  gg.parse_gfa_file(filename);

  // scan over the file to count edges and sequences in parallel
  // -----------------------------------------------------------
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

  std::size_t node_count = line_counts['S'];
  std::size_t edge_count = line_counts['L'];
  std::size_t path_count = line_counts['P'];

  // add nodes
  // ---------
  {
	gg.for_each_sequence_line_in_file(
	  filename,
	  [&](const gfak::sequence_elem& s) {
		dg->create_handle(s.sequence, std::stoll(s.name));
	  });
  }

  // add edges
  // ---------
  {
	gg.for_each_edge_line_in_file(
	  filename,
	  [&](const gfak::edge_elem& e) {
		if (e.source_name.empty()) return;
		// TODO: is this not pointless computation for the sake of following the libhandlegraph the API?
		handlegraph::handle_t a = dg->get_handle(stoll(e.source_name), !e.source_orientation_forward);
		handlegraph::handle_t b = dg->get_handle(stoll(e.sink_name), !e.sink_orientation_forward);
		dg->create_edge(a, b);
	  });	
  }

  dg->compute_start_nodes();
  dg->compute_stop_nodes();


}
  
} // namespace io
