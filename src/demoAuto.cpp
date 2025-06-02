#include <vector>
#include <set>
#include <unordered_map>
#include <iostream>
#include <numeric>
#include <cassert>

struct Graphlet {
  std::vector<int> nodes;
  std::unordered_map<int, std::set<int>> edges;
  void addNode(int label) {
    nodes.push_back(label);
  }
  void addEdge(int from, int to) {
    edges[from].insert(to);
    edges[to].insert(from);
  }
  int label(int node) const {
    return nodes[node];
  }
  std::set<int> neighbors(int node) const {
    return edges.at(node);
  }
  int degree(int node) const {
    return neighbors(node).size();
  }
  size_t size() const {
    return nodes.size();
  }
  Graphlet(std::vector<int> labels, std::vector<std::pair<int, int>> edges) {
    for (auto label : labels) {
      addNode(label);
    }
    for (auto edge : edges) {
      addEdge(edge.first, edge.second);
    }
  }
};

void autos(std::vector<std::vector<int>> &automorphisms_vec,
    std::vector<int> &partial, const std::set<int> &remaining,
    const Graphlet &g) {
  int curr = partial.size();
  if (remaining.empty() && partial.size() == g.size()) {
    automorphisms_vec.push_back(partial);
    return;
  }
  assert(!remaining.empty());
  assert(partial.size() < g.size());

  std::set<int> remcop(remaining);

  auto allPreviouslyMapped = [&curr, &partial, &g](int x) {
    for (int y : g.neighbors(curr)) {
      if (y >= curr) continue;
      if (g.neighbors(x).find(partial[y]) == g.neighbors(x).end()) {
        return false;
      }
    }
    return true;
  };

  for (int x : remaining) {
    if (g.degree(x) != g.degree(curr) || g.label(x) != g.label(curr))
      continue;
    if (!allPreviouslyMapped(x)) continue;

    partial.push_back(x);
    remcop.erase(x);
    autos(automorphisms_vec, partial, remcop, g);
    partial.pop_back();
    remcop.insert(x);
  }
}

std::vector<std::vector<int>> generateAutomorphisms(const Graphlet &g) {
  std::set<int> autorem;
  for (int i = 0; i < g.size(); i++) {
    autorem.insert(i);
  }
  std::vector<int> init;
  std::vector<std::vector<int>> automorphisms;
  autos(automorphisms, init, autorem, g);
  return automorphisms;
}

void reduced(std::vector<int> &ord, std::vector<int> &restrictions,
    const Graphlet &g) {
  std::cout << "reduced: ";
  for (auto node : ord) {
    std::cout << node << " ";
  }
  std::cout << std::endl;
  std::cout << "restrictions: ";
  for (auto node : restrictions) {
    std::cout << node << " ";
  }
  std::cout << std::endl;
  throw std::runtime_error("yet to complete");
}

void process_all(std::vector<int> &restrictions,
    std::vector<std::vector<int>> &remaining_auto,
    std::vector<int> &partial, std::set<int> &used,
    std::set<int> &remaining, const Graphlet &g) {
  if (partial.size() == 2) {
    throw std::runtime_error("do something");
  }
  if (remaining.empty() && partial.size() == g.size()) {
    reduced(partial, restrictions, g);
    return;
  }

  assert(!remaining.empty());
  assert(partial.size() < g.size());

  std::vector<bool> remCare(g.size(), false);
  for (int x : remaining) {
    if (remCare[x]) continue;
    std::vector<int> rests(restrictions);
    std::vector<std::vector<int>> remAutos;

    for (auto automorphism : remaining_auto) {
      remCare[automorphism[x]] = true;
      if (automorphism[x] == x) remAutos.push_back(automorphism);
      else rests[automorphism[x]] = x;
    }

    std::set<int> nused(used);
    std::set<int> nrem(remaining);
    std::vector<int> npartial(partial);
    npartial.push_back(x);
    nused.insert(x);
    nrem.erase(x);
    for (int j : g.neighbors(x)) {
      if (used.find(j) == used.end()) 
        nrem.insert(j);
    }

    process_all(rests, remAutos, npartial, nused, nrem, g);
  }
}

void process(const Graphlet &g,
    const std::vector<std::vector<int>> &automorphisms) {
  std::vector<bool> processed(g.size(), false);

  for (int i = 0; i < g.size(); ++i) {
    if (processed[i]) continue;

    std::vector<int> rests(g.size(), -1);
    std::vector<std::vector<int>> remAutos;
    for (auto automorphism : automorphisms) {
      processed[automorphism[i]] = true;
      if (automorphism[i] == i) remAutos.push_back(automorphism);
      else rests[automorphism[i]] = i;
    }

    std::set<int> used;
    std::set<int> rem;
    std::vector<int> partial;
    partial.push_back(i);
    used.insert(i);
    for (int j : g.neighbors(i))
      rem.insert(j);

    process_all(rests, remAutos, partial, used, rem, g);
  }
}

int main() {

  std::vector<int> labels = {1, 2, 2};
  std::vector<std::pair<int, int>> edges = {{0, 1}, {1, 2}, {2, 0}};

  Graphlet g(labels, edges);
  auto automorphisms = generateAutomorphisms(g);
  std::cout << "labels : \n";
  for (auto label : g.nodes) {
    std::cout << label << " ";
  }
  std::cout << std::endl;
  std::cout << "------------------\n";
  for (auto automorphism : automorphisms) {
    for (auto node : automorphism) {
      std::cout << node << " ";
    }
    std::cout << std::endl;
  }

  process(g, automorphisms);


  return 0;

}
