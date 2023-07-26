#include <algorithm>
#include <cstddef>
#include <format>
#include <iostream>
#include <limits>
#include <set>
#include <stack>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <queue>
#include <map>

#include "./u_graph.hpp"
#include "./digraph.hpp"
#include "./spanning_tree.hpp"
#include "./tree.hpp"


// undirected graph
namespace u_graph {

/*
 * Vertex
 * ======
 */

// constructor(s)
// --------------
Vertex::Vertex() : edge_idxs(std::set<e_v_t>{}) { }
//Vertex::Vertex (std::size_t edge_idx, std::size_t v_idx) : edge_idxs(std::set<e_v_t>{edge_idx}) { }

// getters
// -------
std::set<e_v_t> const& Vertex::get_adjacent_vertices() const {
  return edge_idxs;
}

std::vector<std::size_t> Vertex::edge_indexes() const {
  std::vector<std::size_t> edge_idxs;
  for (auto e : this->edge_idxs) { edge_idxs.push_back(e.e_idx); }
  return edge_idxs;
}

std::vector<std::size_t> Vertex::adj_vertices() const {
  std::vector<std::size_t> v;
  for (auto e : this->edge_idxs) { v.push_back(e.v_idx); }
  return v;
}

// setters
// -------
void Vertex::add_edge_idx(std::size_t e_idx, std::size_t v_idx) {
  this->edge_idxs.insert(e_v_t{e_idx, v_idx});
}

//void Vertex::del_edge_idx(std::size_t e_idx) { this->edge_idxs.erase(e_idx); }

/*
 * Flow Graph
 * ==========
 */

// at least 2
FlowGraph::FlowGraph(std::size_t initial_len) {
  std::vector<Vertex> d(initial_len, Vertex());
  this->adj_list = std::move(d);
  std::size_t edge_idx = this->edges.size(); // == 0

  this->edges.push_back(u_graph::Edge{0, initial_len - 1});

  this->adj_list[0].add_edge_idx(edge_idx, initial_len - 1);
  this->adj_list[initial_len - 1].add_edge_idx(edge_idx, 0);
}

FlowGraph::FlowGraph(digraph::DiGraph const& di_graph) {
  std::size_t size = di_graph.size();
    /*
    handle the start of the graph
   */

  // initialize the flow graph
  // with a dummy start and stop node and connect them
  std::vector<Vertex> d(size+2, Vertex());
  this->adj_list = std::move(d);

  std::size_t edge_idx = this->edges.size();

  this->edges.push_back(u_graph::Edge{0, size + 1});

  // connect start to end
  this->adj_list[0].add_edge_idx(edge_idx, size + 1);
  this->adj_list[size + 1].add_edge_idx(edge_idx, 0);

  // connect all digraph start nodes to flow graph dummy start node
  for (auto const& start_node : di_graph.starts()) {
    this->set_start_node(start_node);
    continue;

    edge_idx = this->edges.size();
    this->edges.push_back(u_graph::Edge{0, start_node+1});
    this->adj_list[0].add_edge_idx(edge_idx, start_node+1);
    // no inc because zero index
    this->adj_list[start_node].add_edge_idx(edge_idx, 0);
  }

  /*
    handle the ends of the graph
   */
  for (std::size_t i{}; i < size; ++i) {
    for (auto const& edge : di_graph.get_vertex(i).out()) {
      this->add_edge(i, edge.to(), edge.get_color());
    }
  }

  /*
    handle the ends of the graph
   */
  // connect all digraph end nodes to stop node
  for (auto const& stop_node : di_graph.stops()) {
    this->set_stop_node(stop_node);
    continue;

    edge_idx = this->edges.size();

    this->edges.push_back(u_graph::Edge{stop_node+1, this->stop_node_internal_idx()});

    this->adj_list[this->stop_node_internal_idx()].add_edge_idx(edge_idx, stop_node+1);
    this->adj_list[stop_node+1].add_edge_idx(edge_idx, this->stop_node_internal_idx());
    //this->set_stop_node(stop_node);
  }
}

FlowGraph::FlowGraph(spanning_tree::Tree& t) {
  //FlowGraph g = FlowGraph(t.size());

  //std::size_t initial_len = t.size();

  auto foo =[&](std::size_t n1, std::size_t n2, core::color c, int weight){
    std::size_t edge_idx = this->edges.size();

    this->edges.push_back( u_graph::Edge{n1, n2, c, weight});
    this->adj_list[n1].add_edge_idx(edge_idx, n2);
    this->adj_list[n2].add_edge_idx(edge_idx, n1);
  };

  std::vector<Vertex> d(t.size(), Vertex());

  this->adj_list = std::move(d);
  //std::size_t edge_idx = this->edges.size(); // == 0


  //this->edges.push_back(u_graph::Edge{0, initial_len - 1});



  //this->adj_list[0].add_edge_idx(edge_idx, initial_len - 1);
  //this->adj_list[initial_len - 1].add_edge_idx(edge_idx, 0);
  //std::size_t edge_idx = this->edges.size(); // == 0

    for (std::size_t j{}; j < t.size(); ++j) {
      std::size_t i = t.get_sorted(j);

      for (auto edge : t.get_child_edges(i)) {
        foo(i, edge.get_child(), edge.get_color(), edge.get_class());
      }

      for (auto bee : t.get_obe_w_id(i)) {
        spanning_tree::BackEdge be = t.get_backedge_given_id(bee.first);
        if (!be.is_capping_backedge()) {
          foo(i, be.get_tgt(), be.get_color(), be.get_class());
        }
      }
    }

    //std::cout << "internal " << this->size_internal() << std::endl;

}

std::size_t FlowGraph::size_internal() {
  return this->adj_list.size();
}

std::size_t FlowGraph::size() {
  return this->size_internal() - 2;
}

std::size_t FlowGraph::stop_node_internal_idx() {
  return this->size_internal() - 1;
}

Vertex& FlowGraph::start_node_internal() {
  return this->adj_list.at(this->start_node_id);
}

Vertex& FlowGraph::stop_node_internal() {
  return this->adj_list.at(this->stop_node_internal_idx());
}

Vertex& FlowGraph::get_vertex_mut(std::size_t vertex) {
  return this->adj_list.at(++vertex);
}

Vertex const& FlowGraph::get_vertex(std::size_t vertex) const {
  return this->adj_list.at(++vertex);
}

Vertex const& FlowGraph::get_vertex_internal(std::size_t vertex) const {
  return this->adj_list.at(vertex);
}

// will break when called with internal vertex indexes
//std::set<std::size_t>const& FlowGraph::get_adjacent_vertices(std::size_t vertex) const {
//  return this->get_vertex(vertex).get_adjacent_vertices();
//}

// will break when called with internal vertex indexes
std::vector<std::size_t>
FlowGraph::get_edge_indexes_internal(std::size_t vertex) const {
  return this->get_vertex_internal(vertex).edge_indexes();
}

void FlowGraph::set_start_node(std::size_t vertex) {

  std::size_t  edge_idx = this->edges.size();
  this->edges.push_back(u_graph::Edge{0, vertex+1});

  this->start_node_internal().add_edge_idx(edge_idx, vertex+1);
  this->get_vertex_mut(vertex).add_edge_idx(edge_idx, 0);
}

void FlowGraph::set_stop_node(std::size_t vertex) {

  std::size_t edge_idx = this->edges.size();
  this->edges.push_back(u_graph::Edge{vertex+1, this->stop_node_internal_idx()});

  this->get_vertex_mut(vertex).add_edge_idx(edge_idx, this->stop_node_internal_idx());
  this->stop_node_internal().add_edge_idx(edge_idx, vertex+1);
}

void FlowGraph::add_edge(std::size_t n1, std::size_t n2, core::color c, int weight, bool inc) {
    if (inc) {
      ++n1; ++n2;
    }
    //++n1; ++n2;
  // move the last node further right if n2 or n1 reaches it
  // TODO: what about self loops i.e n1 == n2?
  std::size_t curr_size = this->size();
  std::size_t greater = std::max(n1, n2);

  // TODO: update current start and stop nodes
  if (curr_size < greater) {
    std::size_t new_internal_size = this->size_internal() + greater - curr_size;

    // update the edge from dummy start to dummy stop
    // is out of range
    for (auto edge_idx: this->start_node_internal().edge_indexes()) {

      Edge& edge = this->edges.at(edge_idx);

      // TODO: check if stop_node is left
      if (edge.right() == this->stop_node_internal_idx()) {
        edge.set_right(new_internal_size - 1);
        //this->start_node_internal().del_adjacent_vertex(e);
        break;
      }
    }

    this->adj_list.resize(new_internal_size, Vertex());

    std::size_t curr_stop_node_idx = curr_size + 1;
    std::swap(this->adj_list[curr_stop_node_idx],
              this->adj_list[new_internal_size - 1]);
  }

  std::size_t edge_idx = this->edges.size();

  //std::cout  << " n1 " << n1 << " n2 " << n2 << " edis " << edge_idx << std::endl;
  this->edges.push_back(u_graph::Edge{n1, n2, c, weight});
  this->adj_list[n1].add_edge_idx(edge_idx, n2);
  this->adj_list[n2].add_edge_idx(edge_idx, n1);
}

spanning_tree::Tree FlowGraph::compute_spanning_tree() {
  spanning_tree::Tree t = spanning_tree::Tree(this->size_internal());

  //std::set<std::size_t> explored;
  std::set<std::size_t> seen;
  std::stack<std::size_t> visited;

  std::size_t current_vertex{this->start_node_id};
  visited.push(current_vertex);

  std::size_t counter{0};

  while (!visited.empty()) {
    current_vertex = visited.top();

    if (!seen.count(current_vertex)) {
      t.set_dfs_num(current_vertex, counter);
      t.set_sort(counter, current_vertex);
      ++counter;
    }

    seen.insert(current_vertex);

    // TODO: simplify below for loop
    // - replace f with not_explored
    // bool not_explored{false}; // the current vertex has not been explored
    bool f{false};
    Vertex const& v =  this->get_vertex_internal(current_vertex);

    // TODO: better condition here for speedup
    for (auto adj : v.get_adjacent_vertices()) {
      Edge& edge = this->edges.at(adj.e_idx);
      std::size_t a = adj.v_idx;

      //if (a < current_vertex) { continue; }

      if (seen.find(a) == seen.end()) {
        t.add_tree_edge(current_vertex, a, edge.get_color());
        visited.push(a);
        f = true;
        break;
      }
      else if (
        !t.is_root(current_vertex) &&
        t.get_parent(current_vertex) != a &&
        !t.has_child(current_vertex, a) &&
        !t.has_ibe(current_vertex, a) &&
        !t.has_obe(current_vertex, a)
      ) {
        // TODO: why the has child and not parent test?
        //std::cout << "adding back edge: " << current_vertex << " -> " << a << std::endl;
        t.add_be(current_vertex, a, false, edge.get_color());
      }
    }

    if (!f) { visited.pop(); }
  }

  return t;
}

void FlowGraph::compute_pst() {
  //std::set<std::size_t> grey;
  std::set<std::size_t> black;
  std::queue<int> q;

  std::size_t current_vertex{this->start_node_id};

  q.push(current_vertex);

  // implement breath first search
  //grey.insert(current_vertex);

  while (!q.empty()) {
    current_vertex = q.front();
    //black.insert(current_vertex);
    std::cout << "current vertex: " << current_vertex << std::endl;
    q.pop();

    Vertex const& v = this->get_vertex_internal(current_vertex);

    for (auto adj : v.get_adjacent_vertices()) {

      Edge& edge = this->edges.at(adj.e_idx);
      std::size_t a = adj.v_idx;


      if (current_vertex == this->start_node_id && a == this->stop_node_internal_idx()) { continue; }



      if (!black.count(a)) {
        std::cout << "\t" << a << " w: " << edge.get_weight() << std::endl;
        q.push(a);
        black.insert(a);
      }
    }

    black.insert(current_vertex);
  }
  //this->get_vertex_internal(current_vertex).set_color(core::color::grey);

}

void construct_pst(std::vector<std::size_t> const& v) {
    tree::Tree t = tree::Tree(30, false);

    std::map<std::size_t, std::size_t> count_map;
    std::map<std::size_t, std::vector<std::size_t>> pos_map;

    for (auto i : v) {
      if (count_map.find(i) == count_map.end()) {
        count_map[i] = 1;
      }
      else {
        count_map[i] += 1;
      }
    }

    for (std::size_t i{}; i < v.size(); ++i) {
      pos_map[v[i]].push_back(i);
    }

    std::vector<std::size_t> temp{};

    std::size_t counter{};
    std::size_t cr{};
    std::size_t cr_idx{};

    auto find = [&](std::size_t q) {
      for (std::size_t j{}; j < temp.size(); ++j) {
        if (q == temp[j]) { return j; }
      }
    };

    std::vector<bool> nesting(pos_map.size(), false);
    
    //
    for (std::size_t i{}; i < v.size(); ++i) {
      if (pos_map[v[i]].size() == 1) { continue; }
      temp = pos_map[v[i]];

      std::size_t res;

      // find
      res = find(i);

      if (nesting[v[i]]) {
        // the last
        cr_idx = t.get_parent(cr_idx);
        nesting[v[i]] = false;

        if ( i+1 < v.size() && v[i] == v[i+1]) {
          std::size_t p = t.get_parent(cr_idx);
          t.add_vertex(p, counter, v[i]);
          cr_idx = counter;
          ++counter;
        }
      }
      else {
        if (i > 0 && v[i] == v[i-1] && res != temp.size() - 1) {
          std::size_t p = t.get_parent(cr_idx);
          
          t.add_vertex(p, counter, v[i]);
          cr_idx = counter;
          ++counter;

          if (v[i] != v[i+1]) {
            nesting[v[i]] = true;
            continue;
          }
        } else if (i > 0 && v[i] == v[i-1] && res == temp.size() - 1) {
          //std::size_t p = t.get_parent(cr_idx);
          //t.add_vertex(p, counter, v[i]);
          //cr_idx = counter;
          //++counter;
          //continue;
        }


        
        if ( v[i] == v[i+1] && res != temp.size() - 2) {
          // happens twice

          t.add_vertex(cr_idx,counter, v[i]);
          cr_idx = counter;
          ++counter;
                      
        }
        else if (v[i] == v[i+1] && res == temp.size() - 2){
          // first of a possible pair?
          // happens once

          t.add_vertex(cr_idx,counter, v[i]);
          ++counter;
          //cr_idx = counter;
          
        }
        else if (v[i] != v[i+1] && i > 0 && v[i] ==v[i-1]) {
          
        }
        else if (v[i] != v[i+1]) {

          if (counter == 0) {
            t.add_vertex(0,0, v[i]);
            cr_idx = 0;
          } else {
            t.add_vertex(cr_idx,counter, v[i]);
          }

          nesting[v[i]] = true;

          cr = v[i];
          cr_idx = counter;
          ++counter;
        } else {
          
        }
      }

    }

    
}

struct Foo {
  std::size_t src;
  std::size_t tgt;
  std::size_t class_;
};
  
void construct_pvst(std::vector<Foo> const& v) {

  tree::Tree t = tree::Tree(30, false);

  //std::map<std::size_t, std::size_t> count_map;
    std::map<std::size_t, std::vector<std::size_t>> pos_map;

    //for (auto j : v) {
    // std::size_t i = j.class_;
    // if (count_map.find(i) == count_map.end()) {
    //   count_map[i] = 1;
    // }
    // else {
    //   count_map[i] += 1;
    // }
    //}

    for (std::size_t i{}; i < v.size(); ++i) {
      //std::size_t i = j.class_;
      pos_map[v[i].class_].push_back(i);
    }

    std::vector<std::size_t> temp{};

    std::size_t counter{};
    std::size_t pr{};
    std::size_t pr_idx{};

    auto find = [&](std::size_t q) {
      for (std::size_t j{}; j < temp.size(); ++j) {
        if (q == temp[j]) { return j; }
      }
    };

    std::size_t res;

    std::vector<bool> nesting(pos_map.size(), false);
    
    for (std::size_t i{}; i < v.size(); i++) {
      std::size_t cl = v[i].class_;
      temp = pos_map[cl];
      if (temp.size() == 1) { continue; }

      // the position of
      // the current element in
      // in temp
      res = find(i);

      if (nesting[cl]) {
        pr_idx = pr_idx == 0 ? 0 : t.get_parent(pr_idx);
        nesting[cl] = false;
      }
      else if (i+1 < v.size() && cl != v[i+1].class_ && res != temp.size() - 1 ) {
        // enter a new region

        if (counter == 0) {
            pr_idx = 0;
        } 

        t.add_vertex(pr_idx,counter, cl);

        nesting[cl] = true;
        pr_idx = counter;        
        ++counter;


      }
      else {
        t.add_vertex(pr_idx,counter, cl);
        ++counter;        
      }
      
    }

    t.print_dot(true);
}


  
spanning_tree::Tree FlowGraph::compute_pst_again() {
  // based on dfs
  spanning_tree::Tree t = spanning_tree::Tree(this->size_internal());

  std::vector<std::size_t> bub_v;
  std::vector<Foo> foo_v;

  std::vector<std::size_t> in_degree(this->size_internal(), 0);
  std::vector<std::size_t> out_degree(this->size_internal(), 0);
  for (std::size_t i{0}; i < this->size_internal(); ++i) {
    for (auto a : this->get_vertex_internal(i).adj_vertices()) {
      if (a < i) { ++in_degree[i]; } else { ++out_degree[i]; }
    }
  }

  std::set<std::size_t> visited;
  std::stack<std::size_t> s;
  std::stack<std::size_t> bi;

  std::size_t current_vertex{this->start_node_id};
  s.push(current_vertex);
  bool explored{true};

  while (!s.empty()) {
    current_vertex = s.top();

    if (in_degree[current_vertex] > 0) {
      while (s.top() != bi.top()) { s.pop(); }
      bi.pop();
      current_vertex = s.top();
    }

    visited.insert(current_vertex);

    explored = true;
    Vertex const& v =  this->get_vertex_internal(current_vertex);

    if (out_degree[current_vertex] > 1) { bi.push(current_vertex); }

    for (auto adj : v.get_adjacent_vertices()) {
      Edge& edge = this->edges.at(adj.e_idx);
      std::size_t a = adj.v_idx;

      if ( visited.find(a) == visited.end()) {
        // we have discovered a new edge current_vertex -- a

        // std::cout << "" << current_vertex << "-- (" << edge.get_weight()  << ") --" << a << std::endl;

        bub_v.push_back(edge.get_weight());
        foo_v.push_back(Foo{.src = current_vertex, .tgt = 1, .class_ = edge.get_weight()});

        s.push(a);
        --in_degree[a];
        explored = false;

        //if (current_vertex == this->start_node_id && a == this->stop_node_internal_idx()) { continue; }
        //handle_edge(edge.get_weight(), current_vertex, a);

        break;
      }
    }

    if (explored) { s.pop(); }
  }

  std::cout << std::endl;

  construct_pst(bub_v);
  construct_pvst(foo_v);

  return t;
}


void FlowGraph::print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TB;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
    "\t{} [label=\"S\", color=\"green\"];\n"
    "\t{} [label=\"E\", color=\"green\"];\n",
    this->start_node_id,
    this->stop_node_internal_idx());

  //std::set<Edge> reported{};
  std::set<std::size_t> reported{};

  for (std::size_t i{}; i < this->size_internal(); i++) {
    for (std::size_t  edge_idx : this->get_edge_indexes_internal(i)) {
      if (reported.count(edge_idx)) { continue; }

      reported.insert(edge_idx);

      //std::cout << "i" << i << " edge idx " << edge_idx << std::endl;

      Edge& edge = this->edges.at(edge_idx);

      std::size_t to = edge.right() == i ? edge.left() : edge.right();

      std::string color_str =
        this->edges.at(edge_idx).get_color() == core::color::black
        ? "black" : "gray";

      if (edge.get_weight() == core::constants::UNDEFINED_INT) {
        std::cout <<
          std::format("\t{} -- {}  [color=\"{}\"];\n", i, to, color_str);
      } else {
        std::cout <<
          std::format("\t{} -- {}  [color=\"{}\", label=\"{}\"];\n",
                      i, to, color_str, edge.get_weight());
      }

      // std::cout <<
      // std::format("\t{} -- {}  [color=\"{}\"];\n", i, to, color_str);
    }
  }

  std::cout << "}" << std::endl;
}

}; // namespace u_graph
