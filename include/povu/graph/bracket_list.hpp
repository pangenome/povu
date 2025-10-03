#ifndef B_LIST_HPP
#define B_LIST_HPP

#include <cstddef>
#include <list>
#include <unordered_map>

"#include "povu/common/constants.hpp"

namespace povu::bracket_list
{
using namespace povu::constants;

class Bracket;
typedef std::list<Bracket> BracketList;

/*
 * Bracket
 * -------
 *  holds metadata about a back edge
 *
 */
class Bracket
{
	std::size_t back_edge_id_;
	std::size_t recent_size_;
	std::size_t recent_class_; // TODO: rename to class?

public:
	// Bracket(std::size_t backedge_id, std::size_t recent_size, std::size_t
	// recent_class);
	Bracket(std::size_t backedge_id);

	std::size_t back_edge_id();

	std::size_t recent_size() const;
	std::size_t recent_class() const;

	void set_recent_size(std::size_t s);
	void set_recent_class(std::size_t c);
};

/*
 * BracketList wrapper
 * -------------------
 *
 *
 */
class WBracketList
{
	BracketList brackets;
	std::unordered_map<std::size_t, std::list<Bracket>::iterator>
		bracket_map;

public:
	WBracketList();
	void push(Bracket br);
	std::size_t size() const;
	Bracket &top();
	BracketList &get_bracket_list();
	std::unordered_map<std::size_t, std::list<Bracket>::iterator>
	get_bracket_map();
	void del(std::size_t backedge_id);
	void concat(WBracketList *other);
};

} // namespace povu::bracket_list

#endif
