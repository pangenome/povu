#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp"
#include "povu/graph/types.hpp"
// #include "povu/io/to_gfa.hpp"
#include <filesystem>
#include <fstream> // for basic_ofstream, operator<<, bas...

namespace povu::io::to_gfa
{
constexpr std::string_view MODULE = "povu::io::to_gfa";
void write_gfa(const bd::VG &g, const std::filesystem::path &fp);
} // namespace povu::io::to_gfa
