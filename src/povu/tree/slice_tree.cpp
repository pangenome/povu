#include "povu/tree/slice_tree.hpp" // for poi

namespace povu::tree::interval_tree
{
std::string to_string(comp_type ct)
{
	switch (ct) {
	case comp_type::NO_OVERLAP:
		return "NO_OVERLAP";
	case comp_type::EXISTS:
		return "EXISTS";
	case comp_type::REPLACE_ALT:
		return "REPLACE_ALT";
	case comp_type::MERGE_EXTEND:
		return "MERGE_EXTEND";
	case comp_type::MERGE_REPLACE:
		return "MERGE_REPLACE";
	case comp_type::INSERT_ALT:
		return "INSERT_ALT";
	case comp_type::CONTAINED:
		return "CONTAINED";
	case comp_type::EXTEND_ALT:
		return "EXTEND_ALT";
	}

	ERR("Unknown comp_type value: {}", static_cast<pt::u8>(ct));
}

std::string to_string(update_type ut)
{
	switch (ut) {
	case update_type::DO_NOTHING:
		return "DO_NOTHING";
	case update_type::INSERT_LEAF:
		return "INSERT_LEAF";
	case update_type::INSERT_ALT:
		return "INSERT_ALT";
	case update_type::REPLACE_ALT:
		return "REPLACE_ALT";
	case update_type::MERGE_EXTEND:
		return "MERGE_EXTEND";
	case update_type::MERGE_REPLACE:
		return "MERGE_REPLACE";
	case update_type::EXTEND_ALT:
		return "EXTEND_ALT";
	}

	ERR("Unknown update_type value: {}", static_cast<pt::u8>(ut));
}

} // namespace povu::tree::interval_tree
