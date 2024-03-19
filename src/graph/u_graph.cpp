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

  //std::cout << "size: " << size << std::endl;

  // initialize the flow graph
  // with a dummy start and stop node and connect them
  std::vector<Vertex> d(size+2, Vertex());
  this->adj_list = std::move(d);
  //std::cout << "size: " << this->adj_list.size() << std::endl;

  std::size_t dummy_start_node_idx = 0;
  std::size_t dummy_stop_node_idx = size+1;
  std::size_t offset = 1;


  std::size_t edge_idx = this->edges.size();

  this->edges.push_back(u_graph::Edge{0, dummy_stop_node_idx});

  // connect start to end
  this->adj_list[0].add_edge_idx(edge_idx, dummy_stop_node_idx);
  this->adj_list[dummy_stop_node_idx].add_edge_idx(edge_idx, 0);

   /*
    handle the start of the graph
   */
  // connect all digraph start nodes to flow graph dummy start node
  for (auto const& start_node : di_graph.starts()) {
    //std::cout << "start node" << start_node << std::endl;
    this->set_start_node(start_node);
    continue;

    //edge_idx = this->edges.size();
    //this->edges.push_back(u_graph::Edge{0, start_node+1});
    //this->adj_list[0].add_edge_idx(edge_idx, start_node+1);
    //this->edges.push_back(u_graph::Edge{0, start_node});
    //this->adj_list[0].add_edge_idx(edge_idx, start_node);
    // no inc because zero index
    //this->adj_list[start_node].add_edge_idx(edge_idx, 0);
  }

  /*
    handle the middle of the graph
   */
  /*
    handle the ends of the graph
   */
  for (std::size_t i{}; i < size; ++i) {
    digraph::Vertex const& v = di_graph.get_vertex(i);
    if (v.get_handle() == "" || (v.out().empty() && v.in().empty())) {
      continue;
    }

    for (auto const& edge : di_graph.get_vertex(i).out()) {
      //std::cout <<"i: " << i << "edge: " << edge.to() << std::endl;
      // add_edge takes two default args: edge weight and inc
      this->add_edge(i, edge.to(), edge.get_color());
    }
  }

  /*
    handle the ends of the graph
   */
  // connect all digraph end nodes to stop node
  for (auto const& stop_node : di_graph.stops()) {
    //std::cout << "stop node" << stop_node << std::endl;
    this->set_stop_node(stop_node);
    continue;

    //edge_idx = this->edges.size();


    //this->edges.push_back(u_graph::Edge{stop_node, this->stop_node_internal_idx()});

    //this->adj_list[this->stop_node_internal_idx()].add_edge_idx(edge_idx, stop_node);
    //this->adj_list[stop_node].add_edge_idx(edge_idx, this->stop_node_internal_idx());


    //this->edges.push_back(u_graph::Edge{stop_node+1, this->stop_node_internal_idx()});

    //this->adj_list[this->stop_node_internal_idx()].add_edge_idx(edge_idx, stop_node+1);
    //this->adj_list[stop_node+1].add_edge_idx(edge_idx, this->stop_node_internal_idx());
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

  for (std::size_t i{}; i < t.size(); ++i) {
    //std::size_t i = t.get_sorted(j);

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

spanning_tree::Tree FlowGraph::compute_spanning_tree_two(digraph::DiGraph const& g) {
  spanning_tree::Tree t = spanning_tree::Tree(this->size_internal());

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

    // std::cout << "current vertex: " << current_vertex << std::endl;

    // TODO: better condition here for speedup
    for (auto adj : v.get_adjacent_vertices()) {
      Edge& edge = this->edges.at(adj.e_idx);
      std::size_t a = adj.v_idx;

      bool x {false};

      if (current_vertex != this->start_node_id && a != this->stop_node_internal_idx()) {
        if  (current_vertex == this->stop_node_internal_idx()) {continue;}
        for (auto e : g.get_vertex(current_vertex - 1 ).out()) {
          if (e.to() == a - 1) {
            x = true;
            break; }
        }
      }
      else if (current_vertex == this->start_node_id && a != this->stop_node_internal_idx()) {
        x = true;
      } else if (current_vertex != this->start_node_id && a == this->stop_node_internal_idx()) {
        x = true;
      } else if (current_vertex == this->start_node_id && a == this->stop_node_internal_idx()) {
        t.add_be(current_vertex, a, edge.get_weight(), false, edge.get_color());
      }

      if (g.starts().count(current_vertex)) {

      }
      else if (g.stops().count(current_vertex)) {

      }
      else {

      }

      if ( x && seen.find(a) == seen.end()) {
        t.add_tree_edge(current_vertex, a, edge.get_weight(), edge.get_color());
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
        t.add_be(current_vertex, a, edge.get_weight(), false, edge.get_color());
      }
    }

    if (!f) { visited.pop(); }
  }

  return t;
}


void construct_pst(std::vector<std::size_t> const& v) {
  tree::Tree t = tree::Tree(v.size(), false);

    std::map<std::size_t, std::vector<std::size_t>> pos_map;
    for (std::size_t i{}; i < v.size(); ++i) {
      pos_map[v[i]].push_back(i);
    }

    std::size_t counter{};
    std::size_t cr_idx{};

    std::vector<std::size_t> temp{};
    auto find = [&](std::size_t q) {
      for (std::size_t j{}; j < temp.size(); ++j) {
        if (q == temp[j]) { return j; }
      }
    };

    std::vector<bool> nesting(pos_map.size(), false);

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
            // t.add_vertex(0,0, v[i]);
            t.get_vertex_mut(0).set_class(v[i]);
            cr_idx = 0;
          } else {
            t.add_vertex(cr_idx,counter, v[i]);
          }

          nesting[v[i]] = true;

          //cr = v[i];
          cr_idx = counter;
          ++counter;
        } else {

        }
      }
    }
    t.print_dot(true);
}

tree::Tree FlowGraph::construct_pst(std::vector<Edge> const& v) const {
  tree::Tree t = tree::Tree(v.size(), false);

    std::map<std::size_t, std::vector<std::size_t>> pos_map;
    for (std::size_t i{}; i < v.size(); ++i) {
      pos_map[v[i].get_weight()].push_back(i);
    }

    std::size_t counter{};
    std::size_t cr_idx{};

    std::vector<std::size_t> temp{};
    auto find = [&](std::size_t q) {
      for (std::size_t j{}; j < temp.size(); ++j) {
        if (q == temp[j]) { return j; }
      }
    };

    std::vector<bool> nesting(pos_map.size(), false);

    for (std::size_t i{}; i < v.size(); ++i) {
      if (pos_map[v[i].get_weight()].size() == 1) { continue; }
      temp = pos_map[v[i].get_weight()];

      std::size_t res;

      // find
      res = find(i);

      if (nesting[v[i].get_weight()]) {
        // the last
        cr_idx = t.get_parent(cr_idx);
        nesting[v[i].get_weight()] = false;

        if ( i+1 < v.size() && v[i] == v[i+1]) {
          std::size_t p = t.get_parent(cr_idx);
          t.add_vertex(p, counter, v[i].get_weight());
          cr_idx = counter;
          ++counter;
        }
      }
      else {
        // why check if i > 0?
        if (i > 0 && v[i] == v[i-1] && res != temp.size() - 1) {
          std::size_t p = t.get_parent(cr_idx);

          t.add_vertex(p, counter, v[i].get_weight());
          cr_idx = counter;
          ++counter;

          if (v[i] != v[i+1]) {
            nesting[v[i].get_weight()] = true;
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

          t.add_vertex(cr_idx,counter, v[i].get_weight());
          cr_idx = counter;
          ++counter;

        }
        else if (v[i].get_weight() == v[i+1].get_weight() && res == temp.size() - 2){
          // first of a possible pair?
          // happens once

          t.add_vertex(cr_idx,counter, v[i].get_weight());
          ++counter;
          //cr_idx = counter;

        }
        else if (v[i].get_weight() != v[i+1].get_weight() && i > 0 && v[i].get_weight() == v[i-1].get_weight()) {

        }
        else if (v[i].get_weight() != v[i+1].get_weight() ) {

          if (counter == 0) {
            // t.add_vertex(0,0, v[i].get_weight());
            t.get_vertex_mut(0).set_class(v[i].get_weight());
            cr_idx = 0;
          } else {
            t.add_vertex(cr_idx,counter, v[i].get_weight());
          }

          nesting[v[i].get_weight()] = true;

          //cr = v[i].get_weight();
          cr_idx = counter;
          ++counter;
        } else {

        }
      }
    }

    //t.print_dot(true);

    return t;
}

tree::Tree FlowGraph::construct_pvst(std::vector<Edge> const& v) const {
  std::size_t last_node_idx = v.back().right() +10;
  tree::Tree t = tree::Tree(last_node_idx, true);

  std::map<std::size_t, std::vector<std::size_t>> pos_map;
  for (std::size_t i{}; i < v.size(); ++i) {
    pos_map[v[i].get_weight()].push_back(i);
  }

  std::size_t res{};
  std::vector<bool> nesting(pos_map.size(), false);

  std::vector<std::size_t> temp{};
  std::size_t counter{};
  std::size_t pr_idx{};

  /*
    the position of the current element in temp
   */
  auto find = [&](std::size_t q) {
    for (std::size_t j{}; j < temp.size(); ++j) {
      if (q == temp[j]) { return j; }
    }
  };

  for (std::size_t i{}; i < v.size(); i++) {
    std::size_t cl = v[i].get_weight();
    temp = pos_map[cl];

    res = find(i);

    if (temp.size() > 1 && nesting[cl]) {
      pr_idx = pr_idx == 0 ? 0 : t.get_parent(pr_idx);
      nesting[cl] = false;
      t.add_vertex(pr_idx,counter, cl);
      t.get_vertex_mut(counter).set_meta(std::to_string(v[i].right()));
      ++counter;
    }
    else if (
      temp.size() > 1 &&
      i+1 < v.size() &&
      cl != v[i+1].get_weight() &&
      res != temp.size() - 1
      ) {
      // enter a new region
      if (counter == 0) {
        // we already have a root
        t.add_vertex(0,++counter, cl);
      } else {
        t.add_vertex(pr_idx,counter, cl);
      }

      t.get_vertex_mut(counter).set_meta(std::to_string(v[i].right()));
      pr_idx = counter;
      ++counter;

      nesting[cl] = true;
    }
    else {
      // matches when temp.size() == 1 and
      // ...

      t.add_vertex(pr_idx,counter, cl);
      t.get_vertex_mut(counter).set_meta(std::to_string(v[i].right()));
      ++counter;
    }
  }

  //t.print_dot(true);

  return t;
}

/**
   DFS with backtracking
 */
std::vector<Edge> FlowGraph::compute_edge_stack() {

  std::vector<std::size_t> in_degree(this->size_internal(), 0);
  std::vector<std::size_t> out_degree(this->size_internal(), 0);

  // populate the vectors in_degree and out_degree
  for (std::size_t i{0}; i < this->size_internal(); ++i) {
    for (auto a : this->get_vertex_internal(i).adj_vertices()) {
      if (a < i) { ++in_degree[i]; } else { ++out_degree[i]; }
    }
  }

  std::set<Edge> added;

  std::set<std::size_t> visited;
  std::stack<std::size_t> s;
  std::stack<std::size_t> bi;

  std::size_t current_vertex{this->start_node_id};
  s.push(current_vertex);
  bool explored{true};

  std::vector<Edge> foo_v;

  std::size_t previous_vertex{};

  //std::size_t furthest{};


  std::map<std::size_t, std::set<std::size_t>> m;

  //std::set<std::size_t> seen;

  while (!s.empty()) {
    current_vertex = s.top();

    if (in_degree[current_vertex] > 0) {
      while (s.top() != bi.top()) {
        m[previous_vertex].insert(s.top());
        s.pop();
      }
      bi.pop();
      current_vertex = s.top();


    }

    previous_vertex = current_vertex;

    visited.insert(current_vertex);

    explored = true;
    Vertex const& v =  this->get_vertex_internal(current_vertex);

    if (out_degree[current_vertex] > 1) { bi.push(current_vertex); }

    for (auto adj : v.get_adjacent_vertices()) {
      Edge& edge = this->edges.at(adj.e_idx);
      std::size_t a = adj.v_idx;

      if ( !m[current_vertex].count(a) && !visited.count(a)

           //&& !seen.count(a)
        ) {
        // we have discovered a new edge current_vertex -- a

        //std::cout << "adding edge " << edge.left() << " -- " << edge.get_weight() << " -- " << edge.right() << std::endl;

        // if (!added.count(edge) ) {}

        //seen.insert(a);

        foo_v.push_back(edge);
        added.insert(edge);

        s.push(a);
        --in_degree[a];
        explored = false;

        break;
      }
    }

    if (explored) { s.pop(); }
  }

  return foo_v;
}

/**
Edge Stack

0-- (0) --1
1-- (7) --2
2-- (8) --4
4-- (8) --6
6-- (10) --9
9-- (10) --11
6-- (9) --10
10-- (9) --11
11-- (8) --12
2-- (11) --12
12-- (7) --16
1-- (2) --3
3-- (2) --5
5-- (1) --7
7-- (1) --13
5-- (3) --8
8-- (4) --13
13-- (6) --14
8-- (5) --14
14-- (2) --15
15-- (2) --16
16-- (0) --17
17-- (0) --18
0-- (0) --18


Edge Stack

0-- (0) --1
1-- (7) --2
2-- (5) --3
3-- (3) --5
4-- (2) --5
4-- (1) --6
6-- (0) --7
0-- (0) --7


*/

void FlowGraph::print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = LR;\n"
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
