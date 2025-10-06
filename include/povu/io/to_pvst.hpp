#ifndef TO_PVST_IO_HPP
#define TO_PVST_IO_HPP

#include <string>      // for string
#include <string_view> // for string_view

#include "povu/common/app.hpp" // for config
#include "povu/graph/pvst.hpp" // for Tree

namespace povu::io::to_pvst
{
constexpr std::string_view MODULE = "povu::io::to_pvst";
namespace pvst = povu::pvst;

void write_pvst(const pvst::Tree &bt, const std::string &base_name,
		const core::config &app_config);
} // namespace povu::io::to_pvst

// add namespace alias for povu::compat
// below tells clang-tidy to skip that specific check for the next line.
// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pv_to_pvst = povu::io::to_pvst;

#endif
