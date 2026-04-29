#include "povu/refs/refs.hpp"

#include <liteseq/gfa.h>

#include "povu/graph/types.hpp"

namespace oza::refs
{

std::string to_string(ref_format_e r)
{
	if (r == ref_format_e::PANSN)
		return "PANSN";
	else
		return "UNDEFINED";
}

char lq_strand_to_char(liteseq::strand s)
{
	return (s == liteseq::strand::STRAND_FWD) ? '>' : '<';
}

ptg::or_e lq_strand_to_pv_or(liteseq::strand s)
{
	return s == lq::strand::STRAND_FWD ? ptg::or_e::forward
					   : ptg::or_e::reverse;
}

}; // namespace oza::refs
