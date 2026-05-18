#ifndef TO_STRUCTURE_EXPORT_IO_HPP
#define TO_STRUCTURE_EXPORT_IO_HPP

#include <cstddef>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "ita/genomics/vcf.hpp"
#include "povu/common/app.hpp"
#include "povu/graph/bidirected.hpp"
#include "povu/graph/pvst.hpp"

namespace mto::to_structure_export
{
inline constexpr std::string_view MODULE = "povu::io::to_structure_export";

namespace bd = povu::bidirected;
namespace pvst = povu::pvst;

class Writer
{
	std::ofstream out_;
	const bd::VG &graph_;
	const std::vector<pvst::Tree> &pvsts_;
	std::map<const pvst::VertexBase *, std::string> pvst_node_ids_;
	std::size_t next_variant_order_{0};
	bool wrote_variant_{false};
	bool finished_{false};

public:
	Writer(const core::config &app_config, const bd::VG &graph,
	       const std::vector<pvst::Tree> &pvsts);
	Writer(const Writer &) = delete;
	Writer &operator=(const Writer &) = delete;
	Writer(Writer &&) = delete;
	Writer &operator=(Writer &&) = delete;
	~Writer() = default;

	void write_variant_calls(iv::VcfRecIdx &vcf_recs);
	void finish();

private:
	void index_pvst_nodes();
	void write_prefix(const core::config &app_config);
	void write_accepted_gfa(const core::config &app_config);
	void write_pvst_boundary_candidates();
	void write_pvst_nodes();
	void write_variant_call(const std::string &chrom, const iv::VcfRec &record);
	[[nodiscard]]
	std::string node_id(std::size_t component_idx, std::size_t node_idx) const;
	[[nodiscard]]
	std::string node_id(const pvst::VertexBase *vertex) const;
};

} // namespace mto::to_structure_export

#endif // TO_STRUCTURE_EXPORT_IO_HPP
