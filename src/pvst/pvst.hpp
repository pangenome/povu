#ifndef PVST_HPP
#define PVST_HPP

#include <filesystem>
#include <cstddef>
#include <vector>
#include <tuple>

#include "../core/core.hpp"
#include "../graph/tree.hpp"


namespace pvst {
// TODO create PVST class that inherits from Tree

tree::Tree compute_pvst(std::vector<core::eq_n_id_t> v, const core::config& app_config);

void to_text(tree::Tree const& pvst, std::filesystem::path const& output_path);

} // namespace pvst
#endif
