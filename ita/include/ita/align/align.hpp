#ifndef IT_ALN_HPP
#define IT_ALN_HPP

#include <iostream>    // for char_traits, basic_ostream, ope...
#include <string>      // for basic_string, string, to_string
#include <string_view> // for operator<<, string_view
#include <vector>      // for vector

#include <quilt/constants.hpp>	 // for
#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/types.hpp>	 // for q

#include "ita/genomics/allele.hpp" // for itn_t
				   //

// #include "povu/common/constants.hpp" // for INF, INVALID_IDX
//  #include "povu/graph/types.hpp"	     // for walk_t

namespace ita::align
{
inline constexpr std::string_view MODULE = "povu::align";

namespace pgt = quilt::types::graph;

/**
 * used in untangling, the level of alignment
 */

enum class aln_level_e {
	step,
	at // allele traversal
};

struct aln_scores_t {
	qt::idx_t match;
	qt::idx_t mismatch;
	qt::idx_t gap_open;
	qt::idx_t gap_extend;
};

struct aln_result_t {
	qt::idx_t score;
	std::string et; // edit transcript
};

class Matrix
{
	std::vector<qt::idx_t> data;
	qt::idx_t row_count_;
	qt::idx_t col_count_;

public:
	Matrix(qt::idx_t row_count, qt::idx_t col_count)
	    : data(row_count * col_count, 0), row_count_(row_count),
	      col_count_(col_count)
	{}

	qt::idx_t &operator()(qt::idx_t row, qt::idx_t col)
	{
		return data[row * this->col_count_ + col];
	}

	qt::idx_t operator()(qt::idx_t row, qt::idx_t col) const
	{
		return data[row * this->col_count_ + col];
	}

	// Print the matrix to an output stream (default is std::cout)
	void print(std::ostream &os = std::cout) const
	{
		for (qt::idx_t i = 0; i < this->row_count_; ++i) {
			for (qt::idx_t j = 0; j < this->col_count_; ++j) {
				os << ((*this)(i, j) == pc::INVALID_IDX
					       ? pc::INF
					       : std::to_string((*this)(i, j)))
				   << " ";
			}
			os << "\n";
		}
	}
};

std::string align(const ia::at_itn &iw, const ia::at_itn &jw,
		  aln_level_e level);

} // namespace ita::align

#endif // IT_ALN_HPP
