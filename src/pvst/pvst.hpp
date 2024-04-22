#ifndef PVST_HPP
#define PVST_HPP

#include <filesystem>
#include <vector>

#include "../core/core.hpp"
#include "../common/common.hpp"
#include "../graph/tree.hpp"


namespace pvst {
// TODO create PVST class that inherits from Tree
using namespace graph_types;
tree::Tree compute_pvst(std::vector<eq_n_id_t> v, const core::config& app_config);

void to_text(tree::Tree const& pvst, std::filesystem::path const& output_path);
} // namespace pvst
#endif
