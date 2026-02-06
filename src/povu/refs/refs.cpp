#include "povu/refs/refs.hpp"
#include <liteseq/gfa.h>

namespace povu::refs
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
}; // namespace povu::refs
