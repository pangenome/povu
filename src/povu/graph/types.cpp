#include "povu/graph/types.hpp"

namespace povu::types::graph
{

std::string_view to_str(v_type_e vt)
{
	switch (vt) {
	case v_type_e::l:
		return "+";
	case v_type_e::r:
		return "-";
	case v_type_e::dummy:
		return "*";
	default:
		return "?";
	}
};

std::ostream &operator<<(std::ostream &os, const v_type_e &vt)
{
	return os << to_str(vt);
}

/*
 * v_end_e
 * ---------
 */
std::string_view to_str(v_end_e ve)
{
	return ve == v_end_e::l ? "+" : "-";
}

std::ostream &operator<<(std::ostream &os, const v_end_e &ve)
{
	return os << to_str(ve);
}

std::string_view to_str(color_e ve)
{
	switch (ve) {
	case color_e::gray:
		return "gray";
	case color_e::black:
		return "black";
	default:
		return "unknown";
	}
}

std::ostream &operator<<(std::ostream &os, const color_e &c)
{
	return os << to_str(c);
}

/*
 * Side and SideID
 * ----
 */
bool operator<(const side_n_id_t &lhs, const side_n_id_t &rhs)
{
	if (lhs.v_idx < rhs.v_idx)
		return true;
	else if (lhs.v_idx == rhs.v_idx)
		return lhs.v_end < rhs.v_end;
	else
		return false;
}

std::ostream &operator<<(std::ostream &os, const side_n_id_t &x)
{
	os << x.v_idx << x.v_end;
	return os;
}

side_n_id_t side_n_id_t::complement() const
{
	return {types::graph::complement(this->v_end), this->v_idx};
}

/*
 * Orientation
 * ------------
 */

std::string_view to_str(or_e o)
{
	return o == or_e::forward ? ">" : "<";
};

// >> and << might be better than + and -
std::ostream &operator<<(std::ostream &os, const or_e &o)
{
	return os << to_str(o);
}

or_e flip(or_e o)
{
	return o == or_e::forward ? or_e::reverse : or_e::forward;
}

/*
 * id and orientation
 * -----------------
 */
std::ostream &operator<<(std::ostream &os, const id_or_t &x)
{
	os << x.orientation << x.v_id;
	return os;
}

bool operator<(const id_or_t &lhs, const id_or_t &rhs)
{
	if (lhs.v_id < rhs.v_id)
		return true;
	else if (lhs.v_id == rhs.v_id)
		return lhs.orientation < rhs.orientation;
	else
		return false;
}

bool operator!=(const id_or_t &lhs, const id_or_t &rhs)
{
	return lhs.v_id != rhs.v_id || lhs.orientation != rhs.orientation;
}

bool operator==(const id_or_t &lhs, const id_or_t &rhs)
{
	return lhs.v_id == rhs.v_id && lhs.orientation == rhs.orientation;
}

std::string to_string(const walk_t &w)
{
	std::string res;
	for (const id_or_t &s : w)
		res += s.as_str();

	return res;
}

} // namespace povu::types::graph
