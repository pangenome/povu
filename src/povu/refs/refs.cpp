#include <liteseq/gfa.h>

namespace povu::refs
{
char lq_strand_to_char(liteseq::strand s)
{
	return (s == liteseq::strand::STRAND_FWD) ? '>' : '<';
}
}; // namespace povu::refs
