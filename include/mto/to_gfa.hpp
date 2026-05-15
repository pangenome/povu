#include <filesystem> // for path

#include "povu/graph/bidirected.hpp" // for bd::VG

namespace mto::to_gfa
{
constexpr std::string_view MODULE = "povu::io::to_gfa";
void write_gfa(const bd::VG &g, const std::filesystem::path &fp);
} // namespace mto::to_gfa
