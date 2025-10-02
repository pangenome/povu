#ifndef FROM_PVST_IO_HPP
#define FROM_PVST_IO_HPP

#include <string>

#include "../../include/graph/pvst.hpp"

namespace povu::io::from_pvst
{
constexpr std::string_view MODULE = "povu::io::from_pvst";

pvst::Tree read_pvst(const std::string &fp);
} // namespace povu::io::from_pvst

// add namespace alias for povu::compat
// below tells clang-tidy to skip that specific check for the next line.
// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pv_frm_pvst = povu::io::from_pvst;

#endif
