#include <filesystem> // for path

#include <oza/graph/bidirected.hpp> // for bd::VG

namespace mto::to_gfa
{
void write_gfa(const bd::VG &g, const std::filesystem::path &fp);
} // namespace mto::to_gfa
