#ifndef COMMON_HPP
#define COMMON_HPP


#include <cstddef>
#include <iostream>
#include <set>
#include <sys/types.h>

namespace graph_types {

/*
 * ===================
 * Graph related types
 * ===================
 */


/*
 * black edge is default
 * TODO: pick a better default
 * gray edge is a bi-edge
 */
enum color { gray, black };

// implement << operator for color
std::ostream& operator<<(std::ostream& os, const color& c);

// Eq class and node id
struct eq_n_id_t {
  std::size_t eq_class;
  std::size_t v_id;
};



/**
  * l (left) 5' or +
  * r (right) 3' or -
  */
// TODO: replace with struct or class to allow methods like complement
enum class VertexEnd {
  l,
  r
};
typedef VertexEnd v_end_t;
typedef VertexEnd v_end;
std::ostream& operator<<(std::ostream& os, const VertexEnd& vt);

/**
 * l (left) 5' or +
 * r (right) 3' or -
 */
enum class VertexType {
    l,
    r,
    dummy
};
typedef VertexType v_type;
std::ostream& operator<<(std::ostream& os, const VertexType& vt);

// Merge path_t and biedged PathInfo into one namespace
struct path_t {
  std::string name; // name as pertains the GFA file
  std::size_t id; // numerical id to be associated with handle
  bool is_circular; // is the path circular?
};

VertexEnd complement(VertexEnd s);


struct side_n_id_t {
  VertexEnd v_end;
  std::size_t v_idx;

  // prototype < operator
  friend bool operator<(const side_n_id_t& lhs, const side_n_id_t& rhs);

  // method complement
  side_n_id_t complement() const;
};
std::ostream& operator<<(std::ostream& os, const side_n_id_t& x);

struct canonical_sese {
  std::size_t start;
  std::size_t end;
  std::set<std::size_t> in_sese;
};

} // namespace graph_types


namespace common_fns {
/**
 * @brief given a bi-Edged index return its index in the bi-Directed graph
 *
 * @param
 * @param
 * @return
 */
std::size_t to_bidirected_idx(std::size_t be_idx, bool has_dummy=true);

/**
 * @brief given a bi-Directed index return the indexes of the left and right vertices in the bi-Edged graph
 *
 * @param
 * @param
 * @return
 */
std::pair<std::size_t, std::size_t> frm_bidirected_idx(std::size_t bd_idx, bool has_dummy=true);

} // namespace common_fns


namespace genomic_types {

  //typedef std::pair<bidirected::VertexEnd, id_t> side_n_id_t; // TODO: replace with struct
  //typedef std::vector<side_n_id_t> subpath_t;
  //typedef std::vector<subpath_t> subpaths_t;

// TODO: make use of this or delete
enum variant_type {
    SNP,
    DEL,
    INS,
    INV,
    DUP,
    CNV,
    BND
};

// TODO which version of VCF is best?
enum output_format {
    VCF, //  currently outputs v4.2
    PAF, // not yet supported
};

}
#endif
