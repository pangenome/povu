#ifndef FROM_PVST_IO_HPP
#define FROM_PVST_IO_HPP

#include <string>      // for string
#include <string_view> // for string_view

#include "povu/graph/pvst.hpp" // for Tree

namespace mto::from_pvst
{
constexpr std::string_view MODULE = "povu::io::from_pvst";

pvst::Tree read_pvst(const std::string &fp);
} // namespace mto::from_pvst

// add namespace alias for povu::compat
// below tells clang-tidy to skip that specific check for the next line.
// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pv_frm_pvst = mto::from_pvst;

#endif
