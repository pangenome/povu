#include "./bracket_list.hpp"

namespace povu::bracket_list {

/*
 * Bracket
 * -------
 */

Bracket::Bracket(std::size_t backedge_id)
  : back_edge_id_(backedge_id), recent_size_(UNDEFINED_SIZE_T), recent_class_(UNDEFINED_SIZE_T){}

std::size_t Bracket::back_edge_id() { return this->back_edge_id_; }
std::size_t Bracket::recent_class() const { return this->recent_class_; }
std::size_t Bracket::recent_size() const { return this->recent_size_; }
void Bracket::set_recent_size(std::size_t s) { this->recent_size_ = s; }
void Bracket::set_recent_class(std::size_t c) { this->recent_class_ = c; }


/*
 * BracketList Wrapper
 * -------------------
 */

WBracketList::WBracketList() {
  this->brackets = std::list<Bracket>{};
  this->bracket_map = std::unordered_map<std::size_t, std::list<Bracket>::iterator>{};
}

std::size_t WBracketList::size() const {
  return this->brackets.size();
}

void WBracketList::push(Bracket br) {
  this->brackets.push_front(br);
  this->bracket_map[br.back_edge_id()] = this->brackets.begin();
}

Bracket& WBracketList::top() {
  return this->brackets.front();
}

void WBracketList::del(std::size_t be_id) {
  // check if bracket_map has the backedge_id
  // TODO: the be_id should always be in the bracket_map, why the need to check?
  if (this->bracket_map.find(be_id) == this->bracket_map.end()) {
    return;
    // throw std::runtime_error("Bracket not found");
  }
  this->brackets.erase(this->bracket_map[be_id]);
  this->bracket_map.erase(be_id);
}

BracketList& WBracketList::get_bracket_list() {
  return this->brackets;
}

void WBracketList::concat(WBracketList* child) {
  this->brackets.splice(this->brackets.begin(), child->get_bracket_list());

  // update the bracket_map by iterating through the list
  for (auto it = this->brackets.begin(); it != this->brackets.end(); ++it) {
    this->bracket_map[it->back_edge_id()] = it;
  }
}

} // namespace povu::bracket_list
