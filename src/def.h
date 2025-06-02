#pragma once

#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <bitset>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <numeric>
#include <utility>

#include <sstream>
#include <cassert>

#include <cmath>

#include <map>

#define ERR_STR "INVALID\n"
#define GLOBAL_AVERAGE_DEGREE 5.0
#define GLOBAL_VERTEX_COUNT 1000.

#ifdef GPU_GENCODE
#define IFGPU(expression) expression
#define NOGPU(expression)
#else
#define IFGPU(expression)
#define NOGPU(expression) expression
#endif

#define LIST_FLAG 1
#define BATCH_FLAG 2
#define COUNT_FLAG 0
#define STREAM_FLAG 4

double expected_size(const std::pair<std::set<int>, std::set<int>> &clause);

class Graphlet;
class DirectedGraphlet;
class ExecutionPlan;
class MultiPlan;
class HyperPlan;
class MultiRed;
class RestPlan;
class MultiRestPlan;
class Graphlet
{
  public:
    int vertices;
    std::vector<std::set<int>> adjacency;
    std::vector<int> labels;
    Graphlet(int n);
    void add_edge(int, int);
    void add_edge(std::pair<int, int>);
    void add_all(const Graphlet &);
    bool connected() const;
    bool isomorphic(const Graphlet &) const;
    int multiplicity() const;
    Graphlet swapped(int, int) const;
    static std::vector<Graphlet> all_connected(int);
    static Graphlet clique(int);
    static Graphlet from_file(std::string);
    std::string toString() const;
    friend class DirectedGraphlet;
};

class DirectedGraphlet
{
  public:
    const int vertices;
    std::vector<std::pair<std::set<int>, std::set<int>>> adjacency;
    DirectedGraphlet();
    DirectedGraphlet(int);
    DirectedGraphlet(const DirectedGraphlet &);
    DirectedGraphlet(const Graphlet &);
    DirectedGraphlet(const ExecutionPlan &);
    bool is_edge_t(int, int) const;
    bool is_edge_f(int, int) const;
    void add_edge(int, int, bool);
    void generate_complement();
    std::vector<std::tuple<int, int, bool>> edgelist() const;
    void reverse_edge(int, int);
    bool isomorphic(const DirectedGraphlet &) const;
    bool cyclic() const;
    bool operator==(const DirectedGraphlet &) const;
    bool undir_automorphism(std::vector<int>) const;
    bool root_symmetric() const;
    int multiplicity(bool) const;
    DirectedGraphlet root_mirror() const;
    DirectedGraphlet canonical();
    std::vector<DirectedGraphlet> acyc_can_vars() const;
    std::string toString() const;
    bool good() const;
};

class ExecutionPlan
{
  public:
    ExecutionPlan();
    ExecutionPlan(
        DirectedGraphlet &, std::vector<int>, std::vector<int>, std::vector<int>);
    std::string toString() const;
    std::string toCode() const;
    double time_complexity() const;
    double data_complexity() const;

    bool same(const std::vector<int> &matchingOrder);

    friend std::ostream &operator<<(std::ostream &out, const ExecutionPlan &data);

    const int vertices;
    bool labelled;
    std::vector<int> mMatchingOrder;
    // restrictions[i] means i's ID should be less than j'd ID. That is, i'd
    // degree should be larger than j's degree.
    std::vector<int> restrictions;//-1 wherever it doesn't exist
                                  // depend indicates the structure of the graphlet
                                  // the first set determines the set intersections to perform and the second
                                  // set the difference operations. For example, ({0,1,3}, {2}) means
                                  // N(v0)&N(v1)&N(v3)-N(v2)
    std::vector<std::pair<std::set<int>, std::set<int>>> depend;
    std::vector<int> labels;
};

class MultiPlan
{
  public:
    std::vector<ExecutionPlan> plans;
    MultiPlan(const std::vector<ExecutionPlan> &);
    void add(ExecutionPlan);
    static int divergence_point(const ExecutionPlan &, const ExecutionPlan &);
    static std::string varname(std::pair<std::set<int>, std::set<int>>);
    static std::vector<std::pair<std::set<int>, std::set<int>>>
      all_priors(std::pair<std::set<int>, std::set<int>>);
    static std::pair<std::string, double> compute_from_priors(
        const std::pair<std::set<int>, std::set<int>> &,
        const std::set<std::pair<std::set<int>, std::set<int>>> &);
    template<typename T> static T &indent(T &, int);
    template<typename stream_type, typename scope_type>
      double gen_codeblock(
          stream_type &, int, const std::vector<bool> &, scope_type &) const;
    template<typename stream_type, typename dp_type, typename scope_type>
      double gen_level(
          stream_type &,
          int,
          const std::vector<bool> &,
          const dp_type &,
          const scope_type &) const;
    std::pair<std::string, double> gen_code() const;
};

class HyperPlan
{

  // std::map<std::tuple<int, int, int, int>, std::string>
  // candidate_implementations; std::vector<std::pair<std::tuple<int, int, int,
  // int>, std::string>> all_implementations;
  std::string prefix;
  int global_plan_idx;

  public:
  std::vector<std::vector<ExecutionPlan>> candidate_plans;
  void add_plan(int, const ExecutionPlan &);
  int sweep(std::string);
  // void generate_pareto(std::string) const;
  std::string greedy() const;

  private:
  void output_plan(std::pair<std::tuple<int, int, int, int>, std::string>);
};

class MultiRed
{
  public:
    MultiRed(Graphlet &);
    MultiRed(Graphlet &, const std::vector<int> &ordering);
    MultiRed(Graphlet &, std::vector<std::vector<ExecutionPlan>> &);

    // void list_automorphisms(std::vector<std::vector<int>>);
    void autos(
        std::vector<std::vector<int>> &, std::vector<int> &, const std::set<int> &);
    void reduced(std::vector<int> &, std::vector<int> &);
    void process_all(
        std::vector<int> &,
        std::vector<std::vector<int>> &,
        std::vector<int> &,
        std::set<int> &,
        std::set<int> &);
    void process_all(
        std::vector<int> &,
        std::vector<std::vector<int>> &,
        std::vector<int> &,
        std::set<int> &,
        const std::vector<int> &);

    [[nodiscard]] std::vector<ExecutionPlan> plans() const;

    // these are functions defined in multip_red_dg
    // void
    // multiplicity_reduced_graph(DirectedGraphlet&,std::vector<ExecutionPlan>&);

  private:
    Graphlet g;
    int n;
    std::vector<ExecutionPlan> mPlans;
    std::vector<std::vector<ExecutionPlan>> *iterplans;
    std::vector<std::vector<int>> automorphisms;
};

class RestSet
{
  public:
    int depth;
    // increasing order
    std::vector<int> ins;
    // increasing order
    std::vector<int> out;
    std::vector<int> restrict;
    std::vector<int> res_chain;
    int label = -1;
    bool tranResChain = false;
    RestSet(
        const std::vector<int> &,
        const std::vector<int> &,
        const std::vector<int> &,
        const int &l = -1);
    bool operator<(const RestSet &) const;

    std::string varname;
    bool valid();
    // get the restriction at the current level
    int restriction() const;
    // computes a parent, if it exists, otherwise returns self;
    RestSet parent() const;
    void append_calc_to_stream(
        std::ostream &,
        int index,
        std::set<RestSet> &,
        const std::string &) const;
    void append_iter_to_stream(std::ostream &, int id) const;
    double data_complexity_ignoring_restrictions() const;
    double time_complexity_ignoring_restrictions() const;
    double time_complexity_ignoring_restrictions(bool, std::set<RestSet> &) const;
    friend std::ostream &operator<<(std::ostream &os, const RestSet &rs);

  private:
    std::string var_name();
};

class RestPlan
{
  public:
    RestPlan(const ExecutionPlan &, int, bool noCM);
    // ExecutionPlan* parent;
    // first element is the set for v1
    std::vector<RestSet> loopons;
    // first element is the list of sets needed from v0
    std::vector<std::set<RestSet>> depends;
    int id() const;
    int vertices() const;
    std::vector<int> labels() const;
    double data_complexity() const;
    double time_complexity() const;
    friend std::ostream &operator<<(std::ostream &os, RestPlan &rp);

  private:
    int mId, mVertices;
    bool labelled;
    std::vector<int> mLabels;
    std::vector<int> rest;
};

class MultiRestPlan
{
  public:
    int depth = 0;
    bool labelled;
    // npls only matters for the root
    //# of plans
    int npls = 0;
    // make a multirestplan of this depth
    MultiRestPlan(int, bool, bool pVerbose = false);
    ~MultiRestPlan();
    // sets to track at this level
    std::set<RestSet> atlev;
    // things to do below
    // std::map<RestSet,MultiRestPlan*> children;
    std::set<int> labels;
    std::map<RestSet, MultiRestPlan *> children;
    std::map<RestSet, int> counters;
    // add an execution plan with a given id;
    void add_rest_plan(RestPlan &);
    void add_ex_plan(ExecutionPlan &, int, bool noCM);
    int maxdepth(const std::map<RestSet, MultiRestPlan *> &pChildren) const;
    std::pair<std::string, std::string> to_code(int, int, int, bool);
    double data_complexity();
    double time_complexity();
    [[nodiscard]] size_t size() const;
    [[nodiscard]] RestPlan getPlan(size_t idx) const;

  private:
    const bool mVerbose;
    std::vector<RestPlan> plans;
    void append_to_stream(std::ostream &, std::ostream &, int, int, int);
    void append_to_stream_nocm(std::ostream &, std::ostream &, std::set<std::string> &);
    void append_to_stream_nocm_rec(std::ostream &, std::ostream &, std::set<std::string> &);
    void instanciateVariables(std::ostream &, std::string &, const RestSet &);
    double rec_data_complexity(double, std::vector<int> &, std::vector<int> &);
    double rec_time_complexity(double, std::vector<int> &, std::vector<int> &);
    void append_initializers(std::ostream &);
    void append_first_loop(std::ostream &, std::string &);
    void append_callback_func(std::ostream &, std::string &, std::string &);
    void append_deintializers(std::ostream &);
};

template<class T> void print_vector(std::vector<T> &v);

template<class T> void print_vector(const std::vector<T> &v);

template<class T> void print_set(std::set<T> &v);


std::ostream &operator<<(std::ostream &os, const Graphlet &g);
std::ostream &operator<<(std::ostream &os, const std::vector<ExecutionPlan> &plans);

std::ostream &operator<<(std::ostream &os, const RestSet &rs);

std::ostream &operator<<(std::ostream &os, const MultiRestPlan &mrp);

std::ostream &operator<<(std::ostream &os, RestPlan &rp);
std::vector<Graphlet> parseGraphList(std::istream &graphs, bool);


Graphlet::Graphlet(int n) : vertices(n), adjacency(n), labels(n, -1) {}

void Graphlet::add_edge(int src, int dst)
{
  assert(src >= 0 && src < vertices);
  assert(dst >= 0 && dst < vertices);
  adjacency.at(src).insert(dst);
  adjacency.at(dst).insert(src);
}

void Graphlet::add_edge(std::pair<int, int> e) { add_edge(e.first, e.second); }

void Graphlet::add_all(const Graphlet &other)
{
  assert(vertices >= other.vertices);
  for (int s = 0; s < other.vertices; s++) {
    for (int d : other.adjacency.at(s)) { add_edge(s, d); }
  }
}

bool Graphlet::connected() const
{
  std::set<int> reachable;
  reachable.insert(0);
  for (int i = 0; i < vertices; i++) {
    std::set<int> to_add;
    for (int s : reachable) {
      to_add.insert(adjacency.at(s).begin(), adjacency.at(s).end());
    }
    reachable.insert(to_add.begin(), to_add.end());
  }
  return reachable.size() == vertices;
}

bool Graphlet::isomorphic(const Graphlet &other) const
{
  if (vertices != other.vertices) return false;
  std::vector<int> relabel(vertices);
  for (int i = 0; i < vertices; i++) relabel.at(i) = i;
  do {
    bool match = true;
    for (int i = 0; match && i < vertices; i++) {
      int local = i;
      int remote = relabel.at(i);
      auto &set_l = adjacency.at(local);
      auto &set_r = other.adjacency.at(remote);
      if (set_l.size() != set_r.size()) {
        match = false;
        continue;
      }
      for (int v : set_l) {
        if (0 == other.adjacency.at(remote).count(relabel.at(v))) {
          match = false;
          break;
        }
      }
    }
    if (match) return true;
  } while (std::next_permutation(relabel.begin(), relabel.end()));
  return false;
}

long factorial(int x)
{
  if (x <= 1)
    return 1;
  else
    return x * factorial(x - 1);
}

int Graphlet::multiplicity() const
{
  // indistinguishable sets
  std::vector<std::set<int>> indis(vertices);
  for (int i = 0; i < vertices; i++) { indis.at(i).insert(i); }
  for (int i = 0; i < vertices; i++) {
    if (indis.at(i).size() == 0) continue;
    for (int j = i + 1; j < vertices; j++) {
      if (indis.at(j).size() == 0) continue;
      int u = *indis.at(i).begin();
      int v = *indis.at(j).begin();
      Graphlet other = swapped(u, v);
      if (adjacency == other.adjacency) {
        indis.at(i).insert(indis.at(j).begin(), indis.at(j).end());
        indis.at(j).clear();
      }
    }
  }
  std::vector<std::set<int>> temp_indis;
  std::copy_if(
      indis.begin(),
      indis.end(),
      std::inserter(temp_indis, temp_indis.end()),
      [](const auto &s) { return s.size() > 0; });
  indis = temp_indis;
  // interchangeable sets of indistinguishable vertices
  std::vector<std::vector<std::set<int>>> inter(indis.size());
  for (int i = 0; i < indis.size(); i++) { inter.at(i).push_back(indis.at(i)); }
  for (int i = 0; i < inter.size(); i++) {
    if (inter.at(i).size() == 0) continue;
    for (int j = i + 1; j < inter.size(); j++) {
      if (inter.at(j).size() == 0) continue;
      const std::set<int> &su = *inter.at(i).begin();
      const std::set<int> &sv = *inter.at(j).begin();
      if (su.size() != sv.size()) continue;
      std::vector<int> vec_u(su.begin(), su.end());
      std::vector<int> vec_v(sv.begin(), sv.end());
      int nverts = su.size();
      Graphlet other = *this;
      for (int k = 0; k < nverts; k++) {
        other = other.swapped(vec_u.at(k), vec_v.at(k));
      }
      if (adjacency == other.adjacency) {
        inter.at(i).insert(
            inter.at(i).end(), inter.at(j).begin(), inter.at(j).end());
        inter.at(j).clear();
      }
    }
  }
  std::vector<std::vector<std::set<int>>> temp_inter;
  std::copy_if(
      inter.begin(),
      inter.end(),
      std::inserter(temp_inter, temp_inter.end()),
      [](const auto &v) { return v.size() > 0; });
  inter = temp_inter;
  int multip = 1;
  for (const auto &s : indis) {// multiplicity from indistinguishable
    multip *= factorial(s.size());
  }
  for (const auto &v : inter) {// multiplicity from interchangeable
    multip *= factorial(v.size());
  }
  return multip;
}

Graphlet Graphlet::swapped(int u, int v) const
{
  std::vector<int> remap(vertices);
  for (int k = 0; k < vertices; k++) {
    remap.at(k) = k;
    if (k == u) remap.at(k) = v;
    if (k == v) remap.at(k) = u;
  }
  Graphlet out(vertices);
  for (int i = 0; i < vertices; i++) {
    for (int j : adjacency.at(i)) { out.add_edge(remap.at(i), remap.at(j)); }
  }
  return out;
}

std::vector<Graphlet> Graphlet::all_connected(int n)
{
  if (n <= 0)
    return {};
  else if (n == 1)
    return { Graphlet(1) };
  else {
    std::vector<Graphlet> graphlets;
    for (Graphlet gs : all_connected(n - 1)) {
      Graphlet gb(n);
      gb.add_all(gs);
      for (uint32_t mask = 1; mask < (1 << (n - 1)); mask++) {
        std::bitset<32> target(mask);
        Graphlet g(gb);
        for (int i = 0; i < n - 1; i++) {
          if (target[i]) g.add_edge(i, n - 1);
        }
        volatile bool is_unique = true;
        //#pragma omp parallel for
        for (int i = 0; i < graphlets.size(); i++) {
          if (is_unique && g.isomorphic(graphlets.at(i))) { is_unique = false; }
        }
        if (is_unique) {
          graphlets.push_back(g);
          // std::cerr << "found " << graphlets.size() << " of size " << n <<
          // "\n";
        }
      }
    }
    return graphlets;
  }
}

Graphlet Graphlet::clique(int N)
{
  Graphlet out(N);
  for (int i = 0; i < N; i++) {
    for (int j = i + 1; j < N; j++) { out.add_edge(i, j); }
  }
  return out;
}

Graphlet Graphlet::from_file(std::string fname)
{
  std::ifstream infile(fname.c_str());
  int N, src, dst;
  infile >> N;
  Graphlet out(N);
  while (infile >> src >> dst) { out.add_edge(src, dst); }
  return out;
}

std::string Graphlet::toString() const
{
  std::ostringstream oss;
  oss << "Graphlet(" << vertices << " [";
  for (int s = 0; s < vertices; s++) {
    for (int d : adjacency.at(s)) {
      if (s < d) oss << "(" << s << "," << d << ")";
    }
  }
  oss << "])";
  return oss.str();
}

std::vector<Graphlet> parseGraphList(std::istream &graphs, bool labels)
{
  int ngraphs;
  graphs >> ngraphs;
  std::vector<Graphlet> res;
  while (ngraphs-- > 0) {
    int vert, edge;
    graphs >> vert >> edge;
    Graphlet g(vert);
    if (labels) {
      for (int i = 0; i < vert; ++i) {
        assert(!graphs.eof());
        graphs >> g.labels[i];
      }
    }
    for (int i = 0; i < edge; ++i) {
      int a, b;
      assert(!graphs.eof());
      graphs >> a >> b;
      g.add_edge(a, b);
    }
    res.push_back(g);
  }
  return res;
}


DirectedGraphlet::DirectedGraphlet() : vertices(0) {}

DirectedGraphlet::DirectedGraphlet(int n_vertices) : vertices(n_vertices), adjacency(n_vertices) {}

DirectedGraphlet::DirectedGraphlet(const DirectedGraphlet &other) :
  vertices(other.vertices), adjacency(other.vertices) {
    for(int src = 0; src < vertices; src++) {
      for(int dst = 0; dst < vertices; dst++) {
        if(other.is_edge_t(src, dst)) add_edge(src, dst, true);
        if(other.is_edge_f(src, dst)) add_edge(src, dst, false);
      }
    }
  }

DirectedGraphlet::DirectedGraphlet(const Graphlet &other) :
  vertices(other.vertices), adjacency(other.vertices) {
    for(int src = 0; src < vertices; src++) {
      auto &d = other.adjacency.at(src);
      for(int dst : d) {
        if(src < dst) add_edge(src, dst, true);
      }
    }
    generate_complement();
  }

DirectedGraphlet::DirectedGraphlet(const ExecutionPlan &other) :
  vertices(other.vertices), adjacency(other.vertices) {
    for(int src = 0; src < vertices; src++) {
      auto &d = other.depend.at(src);
      for(int dst : d.first) add_edge(src, dst, true);
      for(int dst : d.second) add_edge(src, dst, false);
    }
  }

bool DirectedGraphlet::is_edge_t(int src, int dst) const {
  if(src >= vertices || dst >= vertices) return false;
  return adjacency.at(src).first.find(dst) != adjacency.at(src).first.end();
}

bool DirectedGraphlet::is_edge_f(int src, int dst) const {
  if(src >= vertices || dst >= vertices) return false;
  return adjacency.at(src).second.find(dst) != adjacency.at(src).second.end();
}

void DirectedGraphlet::add_edge(int src, int dst, bool exists=true) {
  assert(src != dst);
  assert(!is_edge_t(src, dst));
  assert(!is_edge_f(src, dst));
  assert(!is_edge_t(dst, src));
  assert(!is_edge_f(dst, src));
  if(exists) {
    adjacency.at(src).first.insert(dst);
  } else {
    adjacency.at(src).second.insert(dst);
  }
}

void DirectedGraphlet::generate_complement() {
  for(int src = 0; src < vertices; src++) {
    for(int dst = src+1; dst < vertices; dst++) {
      if(!is_edge_t(src, dst) && !is_edge_t(dst, src)) {
        if(!is_edge_f(src, dst) && !is_edge_f(dst, src)) {
          add_edge(src, dst, false);
        }
      }
    }
  }
}

std::vector<std::tuple<int, int, bool>> DirectedGraphlet::edgelist() const {
  std::vector<std::tuple<int, int, bool>> el;
  el.reserve(vertices * (vertices - 1) / 2);
  for(int src = 0; src < vertices; src++) {
    for(int dst = src + 1; dst < vertices; dst++) {
      if(is_edge_t(src, dst))
        el.push_back(std::make_tuple(dst, src, true));
      if(is_edge_t(dst, src))
        el.push_back(std::make_tuple(src, dst, true));
      if(is_edge_f(src, dst))
        el.push_back(std::make_tuple(dst, src, false));
      if(is_edge_f(dst, src))
        el.push_back(std::make_tuple(src, dst, false));
    }
  }
  assert(el.size() == vertices * (vertices - 1) / 2);
  return el;
}

void DirectedGraphlet::reverse_edge(int src, int dst) {
  generate_complement();
  if(is_edge_t(src, dst)) {
    adjacency.at(src).first.erase(adjacency.at(src).first.find(dst));
    adjacency.at(dst).first.insert(src);
    return;
  }
  if(is_edge_t(dst, src)) {
    adjacency.at(dst).first.erase(adjacency.at(dst).first.find(src));
    adjacency.at(src).first.insert(dst);
    return;
  }
  if(is_edge_f(src, dst)) {
    adjacency.at(src).second.erase(adjacency.at(src).second.find(dst));
    adjacency.at(dst).second.insert(src);
    return;
  }
  if(is_edge_f(dst, src)) {
    adjacency.at(dst).second.erase(adjacency.at(dst).second.find(src));
    adjacency.at(src).second.insert(dst);
    return;
  }
  assert(false);
}

bool DirectedGraphlet::isomorphic(const DirectedGraphlet &other) const {
  if(vertices != other.vertices) return false;
  std::vector<int> relabel(vertices);
  for(int i = 0; i < vertices; i++) relabel.at(i) = i;
  while(std::next_permutation(relabel.begin(), relabel.end())) {
    bool match = true;
    for(int i = 0; match && i < vertices; i++) {
      int local = i;
      int remote = relabel.at(i);
      auto &sets_l = adjacency.at(local);
      auto &sets_r = other.adjacency.at(remote);
      if(sets_l.first.size() != sets_r.first.size()) {
        match = false;
        continue;
      }
      if(sets_l.second.size() != sets_r.second.size()) {
        match = false;
        continue;
      }
      for(int v : sets_l.first) {
        if(!other.is_edge_t(remote, relabel.at(v))) {
          match = false;
          continue;
        }
      }
      for(int v : sets_l.second) {
        if(!other.is_edge_f(remote, relabel.at(v))) {
          match = false;
          continue;
        }
      }
    }
    if(match) return true;
  }
  return false;
}

bool DirectedGraphlet::cyclic() const {
  std::vector<std::set<int>> successors(vertices);
  for(int i = 0; i < vertices; i++) {
    successors.at(i).insert(adjacency.at(i).first.begin(), adjacency.at(i).first.end());
    successors.at(i).insert(adjacency.at(i).second.begin(), adjacency.at(i).second.end());
  }
  for(int round = 0; round < vertices; round++) {
    for(int i = 0; i < vertices; i++) {
      std::set<int> to_add;
      for(int v : successors.at(i)) {
        to_add.insert(adjacency.at(v).first.begin(), adjacency.at(v).first.end());
        to_add.insert(adjacency.at(v).second.begin(), adjacency.at(v).second.end());
      }
      if(to_add.find(i) != to_add.end()) return true;
      successors.at(i).insert(to_add.begin(), to_add.end());
    }
  }
  return false;
}

bool DirectedGraphlet::operator ==(const DirectedGraphlet &other) const {
  if(other.adjacency.size() != vertices) return false;
  for(int i = 0; i < vertices; i++) {
    if(adjacency.at(i) != other.adjacency.at(i)) return false;
  }
  return true;
}

bool DirectedGraphlet::undir_automorphism(std::vector<int> mapping) const {
  assert(mapping.size() >= vertices);
  std::vector<std::set<int>> original(vertices), mapped(vertices);
  for(int src = 0; src < vertices; src++) {
    for(int dst : adjacency.at(src).first) {
      original.at(src).insert(dst);
      original.at(dst).insert(src);
      mapped.at(mapping.at(src)).insert(mapping.at(dst));
      mapped.at(mapping.at(dst)).insert(mapping.at(src));
    }
  }
  for(int src = 0; src < vertices; src++) {
    if(original.at(src) != mapped.at(src)) {
      //std::cout << "mismatch at " << src << "\n";
      return false;
    }
  }
  return true;
}

bool DirectedGraphlet::root_symmetric() const {
  return multiplicity(true) >= 2;
}

int DirectedGraphlet::multiplicity(bool limit=false) const {
  assert(vertices >= 2);
  std::vector<int> mapping(vertices);
  int mul = 0;
  // for each edge bidirectionally
  int src, dst;
  bool exists;
  auto el = edgelist();
  if(limit) el = { std::make_tuple(0, 1, true) };
  for(auto edge : el) {
    std::tie(src, dst, exists) = edge;
    if(!exists) continue;
    std::vector<std::pair<int, int>> es = {
      std::make_pair(src, dst),
      std::make_pair(dst, src)
    };
    for(auto e : es) {
      mapping.at(0) = e.first;
      mapping.at(1) = e.second;
      int pos = 2;
      for(int i = 0; i < vertices; i++) {
        if(i != e.first && i != e.second) mapping.at(pos++) = i;
        assert(pos <= vertices);
      }
      assert(pos == vertices);
      bool has_one = false;
      do {
        if(undir_automorphism(mapping)) {
          mul++;
          has_one = true;
        }
      } while((!limit || !has_one) && std::next_permutation(mapping.begin()+2, mapping.end()));
    }
  }
  return mul;
}

DirectedGraphlet DirectedGraphlet::root_mirror() const {
  DirectedGraphlet temp(*this);
  temp.reverse_edge(0, 1);
  return temp.canonical();
}

DirectedGraphlet DirectedGraphlet::canonical() {
  std::vector<std::set<int>> successors(vertices);
  for(int i = 0; i < vertices; i++) {
    successors.at(i).insert(adjacency.at(i).first.begin(), adjacency.at(i).first.end());
    successors.at(i).insert(adjacency.at(i).second.begin(), adjacency.at(i).second.end());
  }
  for(int round = 0; round < vertices; round++) {
    for(int i = 0; i < vertices; i++) {
      std::set<int> to_add;
      for(int v : successors.at(i)) {
        to_add.insert(adjacency.at(v).first.begin(), adjacency.at(v).first.end());
        to_add.insert(adjacency.at(v).second.begin(), adjacency.at(v).second.end());
      }
      if(to_add.find(i) != to_add.end()) return DirectedGraphlet(0);
      successors.at(i).insert(to_add.begin(), to_add.end());
    }
  }
  DirectedGraphlet out(vertices);
  std::vector<int> map1(vertices, -1), map2(vertices, -1);
  for(int i = 0; i < vertices; i++) {
    int deg = successors.at(i).size();
    assert(deg < vertices);
    if(map1.at(deg) >= 0) return out;
    map1.at(deg) = i;
    map2.at(i) = deg;
  }
  for(int i = 0; i < vertices; i++) {
    int src = map2.at(i);
    for(int v : adjacency.at(i).first) out.add_edge(src, map2.at(v), true);
    for(int v : adjacency.at(i).second) out.add_edge(src, map2.at(v), false);
  }
  return out;
}

std::vector<DirectedGraphlet> DirectedGraphlet::acyc_can_vars() const {
  std::vector<DirectedGraphlet> vars;
  std::vector<int> order(vertices);
  for(int i = 0; i < vertices; i++) order.at(i) = i;
  do {
    DirectedGraphlet dg(vertices);
    for(int i = 0; i < vertices; i++) {
      for(int v : adjacency.at(i).first) {
        if(order.at(i) > order.at(v))
          dg.add_edge(i, v, true);
        else
          dg.add_edge(v, i, true);
      }
      for(int v : adjacency.at(i).second) {
        if(order.at(i) > order.at(v))
          dg.add_edge(i, v, false);
        else
          dg.add_edge(v, i, false);
      }
    }
    volatile bool is_unique = true;
    //    #pragma omp parallel for
    for(int i = 0; i < vars.size(); i++) {
      if(is_unique && dg.isomorphic(vars.at(i))) {
        is_unique = false;
      }
    }
    if(is_unique) {
      vars.push_back(dg.canonical());
    }
  } while(std::next_permutation(order.begin(), order.end()));
  return vars;
}

std::string DirectedGraphlet::toString() const {
  std::ostringstream oss;
  oss << "DirectedGraphlet(";
  auto el = edgelist();
  std::sort(el.begin(), el.end());
  for(auto t : el) {
    int src, dst;
    bool exists;
    std::tie(src, dst, exists) = t;
    oss << "(" << src << "->" << dst << (exists ? " true)" : " false)");
  }
  oss << " " << (cyclic() ? "Cyclic" : "Acyclic");
  bool rs = root_symmetric();
  oss << " " << (rs ? "RS" : "RA");
  int mul = multiplicity();
  //if(rs) mul /= 2;
  oss << " " << mul << ")\n";
  return oss.str();
}

bool DirectedGraphlet::good() const {
  for(auto &adj : adjacency) {
    if(adj.first.size() == 0 && adj.second.size() > 0) return false;
  }
  return true;
}


ExecutionPlan::ExecutionPlan() : vertices(0), labelled(false) {}

ExecutionPlan::ExecutionPlan(
    DirectedGraphlet &g,
    std::vector<int> rest,
    std::vector<int> lab,
    std::vector<int> pMatchingOrder)
  : vertices(g.vertices),
  labelled(true),
  mMatchingOrder(std::move(pMatchingOrder)),
  restrictions(std::move(rest)),
  depend(g.vertices),
  labels(std::move(lab))
{
  for (auto t : g.edgelist()) {
    int src, dst;
    bool exists;
    std::tie(src, dst, exists) = t;
    auto &d = depend.at(dst);
    if (exists)
      d.first.insert(src);
    else
      d.second.insert(src);
  }
  /*
     for(int i=0;i<vertices;++i){
     std::cout<<i<<":";
     for(int x : depend.at(i).first)std::cout<<"+"<<x;
     for(int x : depend.at(i).second)std::cout<<"-"<<x;
     std::cout<<std::endl;
     }
     */
}

std::string ExecutionPlan::toString() const
{
  std::ostringstream oss;
  oss << "s0 = vertices\n";
  for (int i = 0; i < vertices - 1; i++) {
    for (int j = 0; j < i; j++) oss << "  ";
    oss << "for(v" << i << " : s" << i << "){\n";
    // indentation
    for (int j = 0; j < i + 1; j++) oss << "  ";
    const auto &d = depend.at(i + 1);
    oss << "s" << i + 1 << " = ";
    if (d.first.size() == 0) return ERR_STR;
    for (auto it = d.first.begin(); it != d.first.end(); ++it) {
      if (it != d.first.begin()) oss << " & ";
      oss << "N(v" << *it << ")";
    }
    for (auto it = d.second.begin(); it != d.second.end(); ++it) {
      oss << " - N(v" << *it << ")";
    }
    oss << "\n";
  }
  for (int j = 0; j < vertices - 1; j++) oss << "  ";
  oss << "counter += cardinality(s" << vertices - 1 << ")\n";
  for (int j = vertices - 1; j >= 0; --j) {
    for (int i = j; i; --i) oss << "  ";
    oss << "}\n";
  }
  return oss.str();
}

std::string ExecutionPlan::toCode() const
{
  std::ostringstream oss;
  oss << "uint64_t count = 0;" << std::endl;
  oss << "const vidType vertices = g.V();\n{\n";
  oss << "#pragma omp parallel for reduction(+:count) schedule(dynamic,64)\n";
  for (int i = 0; i < vertices - 1; i++) {
    for (int j = 0; j < i; j++) oss << "  ";
    if (i > 0)
      oss << "for(vidType v" << i << " : s" << i << "){\n";
    else
      oss << "for(vidType v0 = 0; v0 < vertices; ++v0) {\n";
    // indentation
    for (int j = 0; j < i + 1; j++) oss << "  ";
    const auto &d = depend.at(i + 1);
    if (i == vertices - 2)
      oss << "count += ";
    else
      oss << "VertexSet s" << i + 1 << " = ";
    bool count = (i == vertices - 2);
    for (auto it = d.second.begin(); it != d.second.end(); ++it) {
      oss << "difference_" << (count ? "num(" : "set(");
      count = false;
    }
    if (d.first.size() == 0) return ERR_STR;
    for (auto it = d.first.begin(); it != d.first.end(); ++it) {
      if (it != d.first.begin()) {
        oss << "intersection_" << (count ? "num(" : "set(");
        count = false;
      }
      // oss << "N(v" << *it << ")"
    }
    // if(d.first.size() == 0) return ERR_STR;
    int rest = restrictions.at(i + 1);
    bool delayedrest = false;
    for (auto it = d.first.begin(); it != d.first.end(); ++it) {
      if (it != d.first.begin()) oss << ", ";
      if (i == 0 && *it == rest) {
        // only one thing
        oss << "bounded(g.N(v" << *it << "),"
          << "v0)";
        if (i == vertices - 2) oss << ".size()";
        continue;
      }
      oss << "g.N(v" << *it << ")";
      if (it != d.first.begin()) {
        // i+1 is restricted by this one
        if (*it == rest) { oss << ",v" << rest; }
        // i+1 is restricted by the original
        if (delayedrest) {
          oss << ",v" << rest;
          delayedrest = false;
        }
        oss << ")";
      } else if (*it == rest) {
        delayedrest = true;
      }
    }
    for (auto it = d.second.begin(); it != d.second.end(); ++it) {
      oss << ", g.N(v" << *it << ")";
      // i+1 is restricted by this one
      if (*it == rest) { oss << ",v" << rest; }
      // i+1 is restricted by the original
      if (delayedrest) {
        oss << ",v" << rest;
        delayedrest = false;
      }
      oss << ")";
    }
    oss << ";\n";
  }
  for (int j = vertices - 2; j >= 0; --j) {
    for (int i = j; i; --i) oss << "  ";
    oss << "}\n";
  }
  oss << "}std::cout<<\"Counted \"<< count <<\" instances.\";";
  return oss.str();
}
double ExecutionPlan::time_complexity() const
{
  double c = 1.0;
  for (int i = 0; i < vertices; i++) { c *= expected_size(depend.at(i)); }
  return c;
}

  std::set<std::pair<std::set<int>, std::set<int>>>
prefixes(const std::pair<std::set<int>, std::set<int>> &clause)
{
  std::set<std::pair<std::set<int>, std::set<int>>> ps;
  int largest = clause.first.size() + clause.second.size() - 1;
  auto copy = clause;
  for (int i = largest; i >= 0; i--) {
    copy.first.erase(i);
    copy.second.erase(i);
    int t_size = copy.first.size();
    int f_size = copy.second.size();
    int both = t_size + f_size;
    if (both <= 1 || t_size <= 0) break;
    ps.insert(copy);
  }
  return ps;
}

double data_cost(
    const std::pair<std::set<int>, std::set<int>> &expr,
    const std::set<std::pair<std::set<int>, std::set<int>>> &priors)
{
  // find first available prior and compute with it
  std::pair<std::set<int>, std::set<int>> selected_prior, current_expr = expr,
    remaining_expr = expr;
  // bool using_prior = false;
  int largest = expr.first.size() + expr.second.size() - 1;
  while (current_expr.first.size() > 0) {
    current_expr.first.erase(largest);
    current_expr.second.erase(largest);
    if (priors.find(current_expr) != priors.end()) {
      selected_prior = current_expr;
      // using_prior = true;
      break;
    }
    largest--;
  }
  for (int x : selected_prior.first) remaining_expr.first.erase(x);
  for (int x : selected_prior.second) remaining_expr.second.erase(x);
  return expected_size(selected_prior)
    + GLOBAL_AVERAGE_DEGREE * remaining_expr.first.size()
    + GLOBAL_AVERAGE_DEGREE * remaining_expr.second.size();
}

double ExecutionPlan::data_complexity() const
{
  std::set<std::pair<std::set<int>, std::set<int>>> stored_sets;
  for (int i = 2; i < vertices; i++) {
    stored_sets.insert(depend.at(i));
    auto ps = prefixes(depend.at(i));
    stored_sets.insert(ps.begin(), ps.end());
  }
  double complexity = data_cost(depend.at(vertices - 1), stored_sets);
  for (int i = vertices - 2; i >= 0; i--) {
    complexity += data_cost(depend.at(i), stored_sets);
    if (i > 0) complexity *= expected_size(depend.at(i - 1));
  }
  return complexity;
}

bool ExecutionPlan::same(const std::vector<int> &matchingOrder)
{
  if (matchingOrder.size() != mMatchingOrder.size())
    return false;
  if (std::equal(mMatchingOrder.begin(), mMatchingOrder.end(), matchingOrder.begin()))
    return true;
  return false;
}

std::ostream &operator<<(std::ostream &out, const ExecutionPlan &data)
{
  out << "vertices: " << data.vertices;
  out << " labelled: " << data.labelled << "\n";
  out << "labels: ";
  std::copy(data.labels.begin(), data.labels.end(),
      std::ostream_iterator<int>(out, " "));
  out << "restrictions: ";
  std::copy(data.restrictions.begin(), data.restrictions.end(),
      std::ostream_iterator<int>(out, " "));
  out << "\n";
  out << "matching order: ";
  std::copy(data.mMatchingOrder.begin(), data.mMatchingOrder.end(),
      std::ostream_iterator<int>(out, " "));
  out << "\n";
  for (const auto &[first, second] : data.depend) {
    out << "f: ";
    std::copy(first.begin(), first.end(),
        std::ostream_iterator<int>(out, " "));
    out << "s: ";
    std::copy(second.begin(), second.end(),
        std::ostream_iterator<int>(out, " "));
    out << "\n";
  }
  out << "\n";

  return out;
}


double expected_size(const std::pair<std::set<int>, std::set<int>> &clause) {
  double p = (double)GLOBAL_AVERAGE_DEGREE / (double)GLOBAL_VERTEX_COUNT;
  return GLOBAL_VERTEX_COUNT * pow(p, clause.first.size()) * pow(1-p, clause.second.size());
}

MultiPlan::MultiPlan(const std::vector<ExecutionPlan> &in) : plans(in) {}

void MultiPlan::add(ExecutionPlan p) {
  plans.push_back(p);
}

int MultiPlan::divergence_point(const ExecutionPlan &one, const ExecutionPlan &two) {
  int max_v = std::min(one.vertices, two.vertices);
  for(int v = 0; v < max_v; v++) {
    if(one.depend.at(v).first != two.depend.at(v).first) return v;
    if(one.depend.at(v).second != two.depend.at(v).second) return v;
  }
  return max_v;
}

std::string MultiPlan::varname(std::pair<std::set<int>, std::set<int>> expr) {
  if(expr.first.size() + expr.second.size() == 0) return "sV";
  std::ostringstream oss;
  int largest = expr.first.size() + expr.second.size();
  for(int i = 0; i < largest; i++) {
    if(expr.first.find(i) != expr.first.end()) oss << "y" << i;
    else if(expr.second.find(i) != expr.second.end()) oss << "n" << i;
    else assert(false);
  }
  return oss.str();
}

std::vector<std::pair<std::set<int>, std::set<int>>> MultiPlan::all_priors(
    std::pair<std::set<int>, std::set<int>> expr) {
  //std::cerr << varname(expr) << " has priors: ";
  std::vector<std::pair<std::set<int>, std::set<int>>> priors;
  int largest = expr.first.size() + expr.second.size() - 1;
  for(int i = largest; i > 0; i--) {
    expr.first.erase(i);
    if(expr.first.size() == 0) break;
    expr.second.erase(i);
    priors.push_back(expr);
  }
  //for(auto e : priors) std::cerr << varname(e) << ", ";
  //std::cerr << "\n";
  return priors;
}

size_t max_sets_in_scope = 0;

std::pair<std::string, double> MultiPlan::compute_from_priors(
    const std::pair<std::set<int>, std::set<int>> &expr,
    const std::set<std::pair<std::set<int>, std::set<int>>> &priors) {
  max_sets_in_scope = std::max(max_sets_in_scope, priors.size());
  // find first available prior and compute with it
  std::pair<std::set<int>, std::set<int>> selected_prior,
    current_expr = expr, remaining_expr = expr;
  bool using_prior = false;
  int largest = expr.first.size() + expr.second.size() - 1;
  while(current_expr.first.size() > 0) {
    current_expr.first.erase(largest);
    current_expr.second.erase(largest);
    if(priors.find(current_expr) != priors.end()) {
      selected_prior = current_expr;
      using_prior = true;
      break;
    }
    largest--;
  }
  for(int x : selected_prior.first) remaining_expr.first.erase(x);
  for(int x : selected_prior.second) remaining_expr.second.erase(x);
  std::ostringstream oss;
  double complexity = 0.0;
  if(using_prior) oss << varname(selected_prior);
  complexity += expected_size(selected_prior);
  for(auto it = remaining_expr.first.begin(); it != remaining_expr.first.end(); ++it) {
    if(using_prior || it != remaining_expr.first.begin()) oss << " & ";
    oss << "g.N(v" << *it << ")";
    complexity += GLOBAL_AVERAGE_DEGREE;
  }
  for(auto it = remaining_expr.second.begin(); it != remaining_expr.second.end(); ++it) {
    oss << " - g.N(v" << *it << ")";
    complexity += GLOBAL_AVERAGE_DEGREE;
  }
  auto s = oss.str();
  if(s.length() == 0) {
    return { "vertices", GLOBAL_VERTEX_COUNT };
  } else {
    return { s, complexity };
  }
}

template<typename T>
T& MultiPlan::indent(T& stream, int level) {
  for(int i = 0; i < level+1; i++) stream << "  ";
  return stream;
}

template<typename stream_type, typename scope_type>
double MultiPlan::gen_codeblock(stream_type &stream, int level,
    const std::vector<bool> &active,
    scope_type &scope) const {
  if(level == 0) return 0;
  scope_type target_expressions;
  for(int plan_id = 0; plan_id < plans.size(); plan_id++) {
    if(!active.at(plan_id)) continue; // find priors for active plans
    assert(level < plans.at(plan_id).vertices);
    for(int lvl_id = 0; lvl_id < plans.at(plan_id).vertices; lvl_id++) {
      auto expr = plans.at(plan_id).depend.at(lvl_id);
      auto expr_priors = all_priors(expr);
      for(auto &e : expr_priors) {
        if(e.first.size() + e.second.size() > 1) target_expressions.insert(e);
      }
    }
    target_expressions.insert(plans.at(plan_id).depend.at(level));
  }
  double complexity = 0.0;
  for(auto &expr : target_expressions) {
    if(scope.count(expr)) continue;
    if(expr.first.size() + expr.second.size() > level) continue;
    auto p = compute_from_priors(expr, scope);
    complexity += p.second;
    indent(stream, level) << "VertexSet " << varname(expr) << " = "
      << p.first << ";\n";
    if(expr.first.size() + expr.second.size() > 1) scope.insert(expr);
  }
  return complexity;
}

template<typename stream_type, typename dp_type, typename scope_type>
double MultiPlan::gen_level(stream_type &stream, int level,
    const std::vector<bool> &active,
    const dp_type &div_pts,
    const scope_type &scope) const {
  std::vector<bool> processed(plans.size(), false);
  scope_type current_scope(scope);
  double complexity = gen_codeblock(stream, level, active, current_scope);
  while(processed != active) {
    int first=-1;
    std::vector<bool> mask(plans.size(), false), term_mask(plans.size(), false);
    for(int i = 0; i < plans.size(); i++) {
      if(active.at(i) && !processed.at(i)) {
        if(plans.at(i).vertices <= level+1) {
          //indent(std::cerr, level) << "adding " << i << " to term_mask\n";
          term_mask.at(i) = true;
          continue;
        }
        if(first < 0) {
          first = i;
        }
        if(div_pts.at(first).at(i) > level) {
          mask.at(i) = true;
        }
      }
    }
    std::vector<bool> all_mask(mask);
    for(int i = 0; i < plans.size(); i++) if(term_mask.at(i)) all_mask.at(i) = true;
    //gen_codeblock(stream, level, all_mask, current_scope);
    for(int i = 0; i < plans.size(); i++) {
      if(term_mask.at(i)) {
        indent(stream, level) << "counter[" << i << "] += "
          << varname(plans.at(i).depend.at(level))
          << ".size();\n";
      }
    }
    if(first >= 0) {
      if(level == 0) {
        indent(stream, level) << "for(vidType v0 = 0; v0 < vertices; v0++) {\n";
      } else {
        indent(stream, level) << "for(vidType v" << level << " : "
          << varname(plans.at(first).depend.at(level)) << ") {\n";
      }
      if(level == 1) indent(stream, level+1) << "if(v0 > v1) continue;\n";
      scope_type inner_scope(current_scope);
      double c = gen_level(stream, level+1, mask, div_pts, inner_scope);
      complexity += c * expected_size(plans.at(first).depend.at(level));
      indent(stream, level) << "}\n";
    }
    for(int i = 0; i < plans.size(); i++) {
      processed.at(i) = processed.at(i) || mask.at(i) || term_mask.at(i);
    }
  }
  return complexity;
}

std::pair<std::string, double> MultiPlan::gen_code() const {
  std::ostringstream oss;
  std::vector<bool> full_mask(plans.size(), true);
  std::vector<std::vector<int>> div_pts(plans.size(), std::vector<int>(plans.size()));
  for(int i = 0; i < plans.size(); i++) {
    for(int j = 0; j < plans.size(); j++) {
      div_pts.at(i).at(j) = divergence_point(plans.at(i), plans.at(j));
    }
  }
  std::set<std::pair<std::set<int>, std::set<int>>> empty_scope;
  double complexity = gen_level(oss, 0, full_mask, div_pts, empty_scope);
  printf("number of sets active was %ld\n", max_sets_in_scope);
  return { oss.str(), complexity };
}


int substring_count(std::string large, std::string small) {
  int count = 0;
  std::string::size_type pos = 0;
  while((pos = large.find(small, pos)) != std::string::npos) {
    count++;
    pos += small.length();
  }
  return count;
}

void HyperPlan::add_plan(int id, const ExecutionPlan &p) {
  while(candidate_plans.size() <= id) {
    candidate_plans.push_back(std::vector<ExecutionPlan>(0));
  }
  candidate_plans.at(id).push_back(p);
}
/*void generate_pareto(std::string sweep_prefix="") const {
  std::cout << "\nPareto-optimal configs:\n";
  std::set<std::tuple<int, int, int, int>> opt_keys;
  for(const auto p : candidate_implementations) {
  auto my_k = p.first;
  std::set<std::tuple<int, int, int, int>> to_remove;
  bool is_opt = true;
  for(const auto co_k : opt_keys) {
  int my_0, my_1, my_2, my_3;
  std::tie(my_0, my_1, my_2, my_3) = my_k;
  int co_0, co_1, co_2, co_3;
  std::tie(co_0, co_1, co_2, co_3) = co_k;
  if(my_0 >= co_0 && my_1 >= co_1 &&
  my_2 >= co_2 && my_3 >= co_3) {
  is_opt = false;
  }
  if(my_0 <= co_0 && my_1 <= co_1 &&
  my_2 <= co_2 && my_3 <= co_3) {
  to_remove.insert(co_k);
  }
  }
  for(const auto co_k : to_remove) opt_keys.erase(co_k);
  if(is_opt) opt_keys.insert(my_k);
  }
  int plan_idx = 0;
//for(const auto my_k : opt_keys) { // TODO here
for(const auto my_p : all_implementations) {
auto my_k = my_p.first; //
auto my_v = my_p.second; //
int n_for, n_set, n_int, n_sub;
std::tie(n_for, n_set, n_int, n_sub) = my_k;
std::ostringstream oss;
oss << "  // Plan using:\t" << n_for << " for-loops\t" << n_set << " VertexSet's\t"
<< n_int << " intersections\t" << n_sub << " differences\n";
//oss << candidate_implementations.at(my_k) << "\n";
oss << my_v << "\n";
if(sweep_prefix.length() > 0) {
std::ostringstream fname;
fname << sweep_prefix << plan_idx << ".hpp";
std::ofstream outfile(fname.str().c_str());
assert(outfile);
outfile << oss.str();
outfile.close();
} else {
std::cout << oss.str();
}
plan_idx++;
}
}*/
void HyperPlan::output_plan(std::pair<std::tuple<int, int, int, int>, std::string> plan) {
  int local_idx;
  //  #pragma omp critical
  {
    local_idx = global_plan_idx++;
  }
  std::ostringstream fname;
  fname << prefix << local_idx << ".hpp";
  std::ofstream outfile(fname.str().c_str());
  assert(outfile);
  int n_for, n_set, n_int, n_sub;
  std::tie(n_for, n_set, n_int, n_sub) = plan.first;
  outfile << "  // Plan using:\t" << n_for << " for-loops\t" << n_set << " VertexSet's\t"
    << n_int << " intersections\t" << n_sub << " differences\n";
  outfile << plan.second;
  outfile.close();
}

/**
 * Inserts the reduced execution plan into the vector
 * Should be a way to have the automorphism already integrated at this point,
 */
void MultiRed::reduced(std::vector<int> &ord, std::vector<int> &restrictions) {
  /*DirectedGraphlet d = dg.reorder(d);
    for(int i=0;i<n;++i){
  //reorder dg
  }*/
  /*
     for(int x : ord)
     std::cout<<x<<",";
     std::cout<<std::endl;

     std::cout<<std::endl;
     */
  // apply ord onto restrictions
  DirectedGraphlet dg(n);
  // add a bunch of edges
  std::vector<int> labels(n);
  std::vector<int> applied_restrictions(n, -1);
  for (int i = 0; i < n; ++i) {
    labels[i] = g.labels.at(ord[i]);
    for (int j = 0; j < n; ++j) {
      if (i == j)
        continue;
      // restriction application
      if (ord[i] == restrictions[ord[j]]) {
        applied_restrictions[j] = i;
        assert(i < j);
        // std::cout<<"added restriction "<< j << " "<<i<<std::endl;
      }
      // edges
      if (i < j) {
        if (g.adjacency[ord[i]].find(ord[j]) != g.adjacency[ord[i]].end()) {
          dg.add_edge(j, i, true);
          // std::cout<<"added edge "<<i<<" "<<j<<std::endl;
        } else {
          dg.add_edge(j, i, false);
        }
      }
    }
  }
  // dg.generate_complement();
  ExecutionPlan ex(dg, applied_restrictions, labels, ord);
  // print ordering and restrictions for now
  mPlans.push_back(ex);
  if (iterplans) {
    for (int i = 0; i < n; ++i) {
      // remove unset
      if (applied_restrictions[i] < 2)
        applied_restrictions[i] = -1;
    }
    ExecutionPlan exi(dg, applied_restrictions, labels, ord);
    iterplans->at(iterplans->size() - 1).push_back(exi);
  }
}

void MultiRed::process_all(std::vector<int> & restrictions,
    std::vector<std::vector<int>> &remaining_auto,
    std::vector<int> &partial, std::set<int> &used,
    const std::vector<int> &remaining) {

  if (partial.size() == 2) {
    if (iterplans) {
      std::vector<ExecutionPlan> plas;
      iterplans->push_back(plas);
    }
  }
  if (remaining.empty() && partial.size() == n) {
    reduced(partial, restrictions);
    return;
  }
  assert(!remaining.empty());
  assert(partial.size() < n);
  std::vector<int> rests(restrictions);
  std::set<int> new_used(used);
  std::vector<int> part(partial);
  std::vector<int> rem(remaining.begin() + 1, remaining.end());
  int currentVertex = remaining.front();
  part.push_back(currentVertex);
  new_used.insert(currentVertex);
  std::vector<std::vector<int>> rem_auts;
  for (std::vector<int> aut : remaining_auto) {
    if (aut[currentVertex] == currentVertex)
      rem_auts.push_back(aut);
    else
      rests[aut[currentVertex]] = currentVertex;
  }

  process_all(rests, rem_auts, part, new_used, rem);
}

/**
 * restrictions is based on our g, used for space reduction
 * can be translated into restrictions for output
 *
 */
void MultiRed::process_all(std::vector<int> &restrictions,
    std::vector<std::vector<int>> &remaining_auto,
    std::vector<int> &partial, std::set<int> &used,
    std::set<int> &remaining) {
  if (partial.size() == 2) {
    if (iterplans) {
      std::vector<ExecutionPlan> plas;
      iterplans->push_back(plas);
    }
  }

  if (remaining.empty() && partial.size() == n) {
    // we've processed, add to the plans
    reduced(partial, restrictions);
    return;
  }
  assert(!remaining.empty());
  assert(partial.size() < n);
  // which ones we have already processed;
  std::vector<bool> remcare(n, false);
  // go in increasing order
  for (int x : remaining) {
    // already processed
    if (remcare[x])
      continue;
    std::vector<int> rests(restrictions);
    std::vector<std::vector<int>> rem_auts;
    // generate equivalence class
    for (std::vector<int> aut : remaining_auto) {
      // member of the class, mark it as taken care of
      remcare[aut[x]] = true;
      // stabilizer group
      if (aut[x] == x)
        rem_auts.push_back(aut);
      // restriction
      else
        rests[aut[x]] = x;
    }
    std::set<int> nused(used);
    std::set<int> rem(remaining);
    std::vector<int> part(partial);
    part.push_back(x);
    nused.insert(x);
    rem.erase(x);
    // add neighbors to remaining
    for (int j : g.adjacency.at(x)) {
      // if it hasn't already been processed, add it to the list
      if (used.find(j) == used.end()) {
        rem.insert(j);
      }
    }
    // recurse
    process_all(rests, rem_auts, part, nused, rem);
  }
}

/**
*/
void MultiRed::autos(
    std::vector<std::vector<int>> &automorphisms_vec,
    std::vector<int> &partial, const std::set<int> &remaining) {
  int curr = partial.size();
  // printf("E1 %d %d %d\n",partial.size(),remaining.size(),*remaining.begin());
  if (remaining.empty() && partial.size() == n) {
    automorphisms_vec.push_back(partial);
    return;
  }
  assert(!remaining.empty());
  assert(partial.size() < n);
  // go in increasing order
  std::set<int> remcop(remaining);
  // for all candidates for the next
  for (int x : remaining) {
    // check if mapping curr to x is a valid plan
    // equal degree and label
    if (g.adjacency.at(x).size() != g.adjacency.at(curr).size() ||
        g.labels[x] != g.labels[curr])
      continue;
    // only works for undirected right now
    bool bad = false;
    // must have the correct edges to ancestors
    for (int y : g.adjacency.at(curr)) {
      if (y < curr) { // y has been mapped already
                      // the permuted y and x must be adjacent
                      // or x would not be a valid next candidate
        if (g.adjacency.at(x).find(partial[y]) == g.adjacency.at(x).end()) {
          bad = true;
          break;
        }
      }
    }
    if (bad)
      continue;
    // if x is valid as a next
    partial.push_back(x);
    remcop.erase(x);
    autos(automorphisms_vec, partial, remcop);
    partial.pop_back();
    remcop.insert(x);
  }
}

// This generates all distinct restricted plans
MultiRed::MultiRed(Graphlet &d)
  : g(d), n(d.vertices), mPlans(), iterplans(nullptr) {
    // make the automorphisms
    std::set<int> autorem;
    for (int i = 0; i < n; ++i)
      autorem.insert(i);
    std::vector<int> init;
    autos(automorphisms, init, autorem);

    // use the automorphisms
    std::vector<bool> cared(n, false);
    for (int i = 0; i < cared.size(); ++i) {
      // already processed under another thing
      if (cared[i])
        continue;

      std::vector<int> rests(cared.size(), -1);
      std::vector<std::vector<int>> rem_auts;
      for (std::vector<int> aut : automorphisms) {
        cared[aut[i]] = true;
        if (aut[i] == i)
          rem_auts.push_back(aut);
        else
          rests[aut[i]] = i;
      }
      std::set<int> used;
      std::set<int> rem;
      std::vector<int> partial;
      partial.push_back(i);
      used.insert(i);
      // have to fill rem with connected to rem, in dg.
      for (int j : g.adjacency.at(i)) {
        rem.insert(j);
      }
      // did not have time to process without recursion,
      // is probably inefficient, who knows.
      process_all(rests, rem_auts, partial, used, rem);
    }
  }
MultiRed::MultiRed(Graphlet &d, const std::vector<int> &ordering)
  : g(d), n(d.vertices), mPlans(), iterplans(nullptr) {
    // make the automorphisms
    std::set<int> autorem;
    for (int i = 0; i < n; ++i)
      autorem.insert(i);
    std::vector<int> init;
    autos(automorphisms, init, autorem);

    int currentVertex = ordering.front();
    std::vector<int> rests(n, -1);
    std::vector<std::vector<int>> rem_auts;
    for (std::vector<int> aut : automorphisms) {
      if (aut[currentVertex] == currentVertex)
        rem_auts.push_back(aut);
      else
        rests[aut[currentVertex]] = currentVertex;
    }
    std::set<int> used;
    std::set<int> rem;
    std::vector<int> partial;
    partial.push_back(currentVertex);
    used.insert(currentVertex);
    // have to fill rem with connected to rem, in dg.
    for (int j : g.adjacency.at(currentVertex)) {
      rem.insert(j);
    }
    std::vector<int> new_ordering{ordering.begin()+1, ordering.end()};
    // did not have time to process without recursion,
    // is probably inefficient, who knows.
    process_all(rests, rem_auts, partial, used, new_ordering);
  }
MultiRed::MultiRed(Graphlet &d,
    std::vector<std::vector<ExecutionPlan>> &ip)
  : g(d), n(d.vertices), mPlans(), iterplans(&ip) {
    // make the automorphisms
    std::set<int> autorem;
    for (int i = 0; i < n; ++i)
      autorem.insert(i);
    std::vector<int> init;
    autos(automorphisms, init, autorem);
    /*
       for(std::vector<int> aut: automorphisms){
       for(int i : aut){
       std::cout<<i;
       }
       std::cout<<std::endl;
       }
       */
    // printf("entered \n");
    // use the automorphisms
    std::vector<bool> cared(n, false);
    for (int i = 0; i < n; ++i) {
      // already processed under another thing
      if (cared[i])
        continue;

      std::vector<int> rests(n, -1);
      std::vector<std::vector<int>> rem_auts;
      for (std::vector<int> aut : automorphisms) {
        cared[aut[i]] = true;
        if (aut[i] == i)
          rem_auts.push_back(aut);
        else
          rests[aut[i]] = i;
      }
      std::set<int> used;
      std::set<int> rem;
      std::vector<int> partial;
      partial.push_back(i);
      used.insert(i);
      // have to fill rem with connected to rem, in dg.
      for (int j : g.adjacency.at(i)) {
        rem.insert(j);
      }
      // did not have time to process without recursion,
      // is probably inefficient, who knows.
      process_all(rests, rem_auts, partial, used, rem);
    }
  }

std::vector<ExecutionPlan> MultiRed::plans() const {
  return mPlans;
}

/*
   void multiplicity_reduced_graph(DirectedGraphlet& dg,std::vector<ExecutionPlan>&
   planlist){
//find automorphisms
int n = dg.vertices;
std::vector<std::vector<int>> automorphisms;
//start with an initial partial permutation
list_automorphisms(automorphisms);
//see if the partial is ok
//add stuff
//schedule reduction



//apply all permutations - already permutation reduced


}
*/

#include "def.h"
#include <iostream>

// restrictions should be the restrction vector of an execution plan
RestSet::RestSet(
    const std::vector<int> &i,
    const std::vector<int> &o,
    const std::vector<int> &restrictions,
    const int &l)
  : ins(i), out(o), label(l)
{
  depth = 0;
  if (!i.empty()) depth = std::max(depth, i[i.size() - 1]);
  if (!o.empty()) depth = std::max(depth, o[o.size() - 1]);
  // std::copy(restrictions.begin(), restrictions.begin()+depth+2,
  // back_inserter(restrict));
  restrict = restrictions;
  res_chain = std::vector<int>(depth + 1, -1);
  int t = depth + 1;
  while (restrictions[t] != -1) {
    //     std::cerr << "t: " << t << std::endl;
    res_chain[t - 1] = restrictions[t];
    //     std::cerr << "res_chain[t-1]: " << res_chain[t-1] << std::endl;
    t = res_chain[t - 1];
    if (t <= 0) break;
  }
  //   std::cout << "res_chain: (";
  //   std::copy(res_chain.begin(), res_chain.end(),
  //   std::ostream_iterator<int>(std::cout, " ")); std::cout << ")\n";
  varname = var_name();
}

bool RestSet::operator<(const RestSet &other) const
{
  if (label != other.label) return label < other.label;
  // ins size, outs size, ins values, outs values, restriction.
  if (other.ins.size() > ins.size()) return true;
  if (other.ins.size() < ins.size()) return false;
  if (other.out.size() > out.size()) return true;
  if (other.out.size() < out.size()) return false;
  if (other.res_chain.size() < res_chain.size()) return true;
  if (other.res_chain.size() > res_chain.size()) return false;

  for (int i = 0; i < ins.size(); ++i) {
    if (other.ins.at(i) < ins.at(i)) return true;
    if (other.ins.at(i) > ins.at(i)) return false;
  }
  for (int i = 0; i < out.size(); ++i) {
    if (other.out.at(i) < out.at(i)) return true;
    if (other.out.at(i) > out.at(i)) return false;
  }
  for (int i = 0; i < res_chain.size(); ++i) {
    if (res_chain.at(i) < other.res_chain.at(i)) return true;
    if (res_chain.at(i) > other.res_chain.at(i)) return false;
  }
  return false;
}

std::string RestSet::var_name()
{

  std::ostringstream oss;
  int ini = 0;
  int oui = 0;
  int i = 0;
  for (i = 0; i < ins.size() + out.size(); ++i) {
    if (oui < out.size() && out[oui] == i) {
      oss << "n" << out[oui];
      ++oui;
    } else {
      oss << "y" << ins[ini];
      ++ini;
    }
    if (res_chain[i] >= 0) { oss << "f" << res_chain[i]; }
  }
  if (label != -1) oss << "l" << label;
  return oss.str();
}

int RestSet::restriction() const { return res_chain[depth]; }

// adjusting this function can change how performance ends up going.
RestSet RestSet::parent() const
{
  RestSet parent = *this;
  if (!parent.tranResChain) {
    parent.tranResChain = true;
    std::vector<int> t(parent.res_chain.size(), -1);
    for (auto r : parent.res_chain)
      if (r != -1) t[r] = r;
    parent.res_chain = t;
  }
  assert(!ins.empty());
  // std::cout<<"CALLING PARENT FOR "<<varname<<std::endl;
  // highest variable value remaining
  int lastvar = ins[ins.size() - 1];//
  if (out.size() > 0) lastvar = std::max(lastvar, out[out.size() - 1]);
  if (ins[ins.size() - 1] == lastvar) {
    parent.ins.pop_back();
    parent.res_chain.pop_back();
    parent.depth--;
    parent.varname = parent.var_name();
    return parent;
  } else {
    parent.out.pop_back();
    parent.res_chain.pop_back();
    parent.depth--;
    parent.varname = parent.var_name();
    return parent;
  }
}

bool RestSet::valid() { return ins.size() > 0; }

// available is what is available at that level already
void RestSet::append_iter_to_stream(std::ostream &oss, int id) const
{
  if (label == -1 && id == 0) {
    oss << "for(vidType v0 = 0; v0 < vertices;++v0)";
    return;
  }
  oss << "for(vidType v" << id << ":" << varname << ")";
}
void RestSet::append_calc_to_stream(
    std::ostream &oss,
    int index,
    std::set<RestSet> &available,
    const std::string &tabdep) const
{
  // static int transient_id = 100;
  // starting set, l0
  const std::string labelStr = (label != -1) ? "l" + std::to_string(label) : "";
  if (ins.empty() && out.empty()) {
    if (label != -1)
      oss << "auto " << labelStr << " = g.L(" << label << ");\n";
    // if this is empty it's just iterating over everything anyway
    return;
  }
  // just a starting set, y0/y6l0
  if (ins.size() == 1 && out.empty()) {
    if (res_chain[depth] != -1) {
      oss << "auto &" << varname;
      oss << " = y" << ins[0] << labelStr;
      oss << "; //bounded\n";
    } else {
      oss << "auto " << varname;
      oss << " = g.N(v" << ins[0];
      if (label != -1) oss << "," << label;
      oss << ");\n";
    }
    return;
  }
  RestSet par = parent();
  // not sure whether this will be executed in the new version
  if (-1 == index && res_chain[depth] != -1) {
    std::vector<int> rs(restrict);
    RestSet testpar(ins, out, rs, label);
    testpar.res_chain[depth] = -1;
    testpar.varname = testpar.var_name();
    if (available.count(testpar)) {
      oss << "auto " << varname << " = bounded(" << testpar.varname << ",v"
        << res_chain[depth] << ");\n";
      return;
    }
  }
  // std::cout<<varname<<" is a child of "<<par.varname<<std::endl;
  if (par.ins.size() == 0) {// has transient unnamed intermediate sets
                            // bounded for 0
    if (par.out.size() == 0) {
      // oss << "VertexSet " << varname << " = ";
      // y0l0
      if (res_chain[depth] == 0) {
        oss << "auto y0 = g.N(v0";
        if (label != -1) oss << "," << label;
        oss << "); ";
        oss << "auto " << varname << " = bounded(y0" << labelStr;
        oss << ",v0);\n";
      } else {
        oss << "auto y0 = g.N(v0);\n";
      }
      return;
    }
    // we can't use the parent, it doesn't exist
    // just compute ourselves
    /*if(varname == "n0y1") {
      std::cerr << "strange res_chain: ";
      std::copy(res_chain.begin(), res_chain.end(),
      std::ostream_iterator<int>(std::cerr, " ")); std::cerr << std::endl;
      }*/
    if (!out.empty()) {
      std::string tmpLastvarname = "y" + std::to_string(ins[0]) + "n" + std::to_string(out[0]);
      if (out.size() == 1) {
        oss << "auto " << varname;
      } else {
        oss << "auto " << tmpLastvarname << labelStr << "_t";
      }
      oss << " = SetOp::difference_set MATCHING ("
        << "y" << ins[0] 	<< labelStr << ", "
        << "y" << out[0] << labelStr << ");\n";
      for (int i = 1; i < out.size(); ++i) {
        std::string newtmpLastvarname = tmpLastvarname + "n" + std::to_string(out[i]);
        if (i+1 < out.size()) {
          oss << tabdep << "auto " << newtmpLastvarname << labelStr << "_t";
        } else {
          oss << tabdep << "auto " << varname;
        }
        oss  << " = MATCHING_DIFF MATCHING ("
          << tmpLastvarname << labelStr << "_t" << ", "
          << "y" << out[i] << labelStr << ");\n";
        std::swap(tmpLastvarname, newtmpLastvarname);
      }
    }
    /*
       if (out.size() != 0) { oss << "auto " << varname << " = "; }
       for (int i = 0; i < out.size(); ++i) {
       oss << "SetOp::difference_set"
       << " MATCHING "
       << "(";
       }
       oss << "y" << ins[0] << labelStr;
       for (int i = 0; i < out.size(); ++i) {
       oss << ", y" << out[i] << labelStr;
       if (i == 0 && res_chain[depth] != -1) {
       oss << ", v" << res_chain[depth];
       }
       IFGPU(oss << ", local_indices, local_search, &manager, work_group");
       oss << ")";
       }
       if (index != -1 && out.size() == 0) oss << ".size()";
       oss << ";\n";
       */
    return;
  }

  // whether or not we are just counting here
  // if(index== -1) {
  //   oss << "VertexSet " << varname << " = ";
  // } else {
  //   oss<<"counter["<<index<<"] += ";
  // }
  oss << "auto " << varname << " = ";
  // this set is a difference
  if (out.size() != par.out.size()) {
    // oss<< "difference_"<< (index!=-1 ? "num(" : "set(")
    oss <<"SetOp::difference_set"
      << " MATCHING "
      << "("
      << par.varname << ", y" << out[out.size() - 1] << labelStr;
    if (res_chain[depth] != -1) { oss << ", v" << res_chain[depth]; }
    IFGPU(oss << ", local_indices, local_search, &manager, work_group");
    oss << ");\n";
    return;
  }
  // this set is an intersection
  if (ins.size() != par.ins.size()) {

    // oss<< "intersection_" << (index!=-1?"num(":"set(")
    oss << "SetOp::intersection_set"
      << " MATCHING "
      << "("
      << par.varname << ", y" << ins[ins.size() - 1] << labelStr;

    if (res_chain[depth] != -1) { oss << ", v" << res_chain[depth]; }
    IFGPU(oss << ", local_indices, local_search, &manager, work_group");
    oss << ");\n";
    if (index != -1) {
      // oss << "for(vidType v" << "lastone" << ":" << varname << ") {";
      // oss << "add here \n";
    }
    return;
  }
}

double RestSet::data_complexity_ignoring_restrictions() const
{
  std::set<int> inset(ins.begin(), ins.end());
  std::set<int> outset(out.begin(), out.end());
  std::pair<std::set<int>, std::set<int>> dar(inset, outset);
  double xp = expected_size(dar);
  if (std::isnan(xp)) {
    std::cout << "Nan from size" << ins.size() << " " << out.size() << " "
      << res_chain[depth] << std::endl;
  }
  return expected_size(dar);
}

double RestSet::time_complexity_ignoring_restrictions() const
{
  std::set<RestSet> emp;
  return time_complexity_ignoring_restrictions(false, emp);
}
double RestSet::time_complexity_ignoring_restrictions(
    bool numb, std::set<RestSet> &available) const
{
  if (ins.size() <= 1 && out.size() == 0) return 1;
  /*std::set<int> inset(ins.begin(),ins.end());
    std::set<int> outset(out.begin(),out.end());
    std::pair<std::set<int>,std::set<int>> dar(inset,outset);
    */
  if (!numb && res_chain[depth] != -1) {
    std::vector<int> rs(restrict);
    rs[depth] = -1;
    RestSet testpar(ins, out, rs);
    if (available.count(testpar) == 1) {
      // basically nothing, we just do it bounded
      return 0;// std::min(64,restset.data_complexity_ignoring_restrictions());
    }
  }
  RestSet par = parent();
  if (par.ins.size() == 0) {
    // we compute here.
    // global average degree time 1 + expected size one intersect + expected
    // size 2+ ... expected size num outs -1
    // = (global average degree - expected size numouts intersect)/(1-(1-p))
    double p = (double)GLOBAL_AVERAGE_DEGREE / (double)GLOBAL_VERTEX_COUNT;
    return GLOBAL_AVERAGE_DEGREE * (1 - pow(1 - p, out.size())) / p;
  }

  //  std::cout<<data_complexity_ignoring_restrictions();
  return parent().data_complexity_ignoring_restrictions();
}


RestPlan::RestPlan(const ExecutionPlan &ep, int _id, bool noCM)
  : mId(_id),
  mVertices(ep.vertices),
  labelled(ep.labelled),
  mLabels(ep.labels),
  rest(ep.restrictions)
{
  // construct the plans at each level for what we loop on
  // first element on loopon is the set for v1
  // first construct what they depend on, put them in
  // This loop implements single-plan set operation motion
  // std::cout<<"whoop"<<std::endl;
  if (noCM && false) { // not working for now
    for (int i = 1; i < ep.vertices; ++i) {
      std::vector<int> ins(ep.depend[i].first.begin(), ep.depend[i].first.end());
      std::vector<int> out(ep.depend[i].second.begin(), ep.depend[i].second.end());
      loopons.emplace_back(ins, out, ep.restrictions, ep.labels[i]);
      depends.emplace_back();
      if(i != ep.vertices-1) depends[i-1].insert(loopons[i-1]);

      RestSet curr = loopons[i-1].parent();
      auto depth = curr.depth+1;
      while(!curr.ins.empty()) {
        depends[depth].insert(curr);
        curr = curr.parent();
      }
    }
  } else {
    for (int i = 0; i < ep.vertices; ++i) {
      std::vector<int> ins(ep.depend[i].first.begin(), ep.depend[i].first.end());
      std::vector<int> out(ep.depend[i].second.begin(), ep.depend[i].second.end());
      //     std::vector<int> rs(ep.restrictions.begin()+1,
      //     ep.restrictions.end()); loopons.emplace_back(ins,out,rs);
      loopons.emplace_back(ins, out, ep.restrictions, ep.labels[i]);
      depends.emplace_back();
      // don't include the last one
      // std::cout<<i<<" "<<ep.vertices<<std::endl;
      for (int j = 0; j < i; ++j) {
        std::vector<int> jin(1, j);
        std::vector<int> jout;
        std::vector<int> jres(mVertices, -1);
        RestSet temp(jin, jout, jres, mLabels[i]);
        depends[j + 1].insert(temp);
      }
      depends[i].insert(loopons[i]);
      if (i > 0) {
        RestSet curr = loopons[i].parent();
        while (!curr.ins.empty()) {
          // depends is offset to have
          depends[curr.depth + 1].insert(curr);
          curr = curr.parent();
        }
      }
      // std::cout<<depends[i-1].size()<<std::endl;
    }
  }
}

int RestPlan::vertices() const
{
  return mVertices;
}

int RestPlan::id() const {
  return mId;
}

std::vector<int> RestPlan::labels() const
{
  return mLabels;
}

// expected memory allocated to vertex sets
double RestPlan::data_complexity() const
{
  // go through depends and loopons, iterating down
  std::vector<int> multcounts(mVertices);
  // for(int i=0;i<vertices;++i)multcounts[0] = 0;
  double res = 0;

  // expected times a loop happens is
  // (normal expected times) / product (multcounts)
  double expectedloopness = GLOBAL_VERTEX_COUNT;
  for (int i = 0; i < mVertices - 1; ++i) {
    // go down on loopon
    int mc = i;
    while (mc != -1) {
      ++multcounts[mc];
      mc = rest[mc];
    }
    double inexpect = expectedloopness;
    for (int j = 0; j <= i; ++j) inexpect /= multcounts[j];
    for (RestSet s : depends[i]) {
      res += s.data_complexity_ignoring_restrictions() * inexpect;
    }
    expectedloopness *= loopons[i].data_complexity_ignoring_restrictions();
  }
  // delete[] multcounts;
  // std::cout<<res<<std::endl;
  return res;
}

// expected time to compute
double RestPlan::time_complexity() const
{
  // go through depends and loopons, iterating down
  std::vector<int> multcounts(mVertices);
  // int* multcounts = new int[vertices];
  for (int i = 0; i < mVertices; ++i) multcounts[i] = 0;
  double res = 0;

  // expected times a loop happens is
  // (normal expected times) / product (multcounts)
  double expectedloopness = GLOBAL_VERTEX_COUNT;
  for (int i = 0; i < mVertices - 1; ++i) {
    // go down on loopon
    int mc = i;
    while (mc != -1) {
      ++multcounts[mc];
      mc = rest[mc];
    }
    double inexpect = expectedloopness;
    for (int j = 0; j <= i; ++j) {
      // std::cout<<j<<"has mult "<<multcounts[j]<<" at step "<<i<<std::endl;
      inexpect /= multcounts[j];
    }
    // double oldres = res;
    for (const RestSet& s : depends[i]) {
      res += s.time_complexity_ignoring_restrictions() * inexpect;
    }
    if (i == mVertices - 2) {
      // std::cout<<" Adding complexity of the loop "<<std::endl;
      res += loopons[i].time_complexity_ignoring_restrictions() * inexpect;
    }
    //    std::cout<<i<<" makes " <<oldres<<" res "<<res<<"
    //    "<<res-oldres<<std::endl;
    expectedloopness *= loopons[i].data_complexity_ignoring_restrictions();
  }

  return res;
}

MultiRestPlan::MultiRestPlan(int dep, bool l, bool pVerbose)
  : depth(dep), labelled(l), labels(), mVerbose(pVerbose)
{
  // everything else should be initiated automatically
}
MultiRestPlan::~MultiRestPlan()
{
  // delete children
  for (auto smsm : children) { delete smsm.second; }
}
void MultiRestPlan::add_ex_plan(ExecutionPlan &ep, int id, bool noCM)
{
  RestPlan rp(ep, id, noCM);
  add_rest_plan(rp);
}

void MultiRestPlan::add_rest_plan(RestPlan &rp)
{
  plans.push_back(rp);
  MultiRestPlan *curr = this;
  int dep = 0;
  for (; dep < rp.loopons.size() - 1; ++dep) {
    if (mVerbose) {
      std::cout << "putting into " << dep << "\n";
      std::copy(
          rp.depends[dep].begin(),
          rp.depends[dep].end(),
          std::ostream_iterator<RestSet>(std::cout, " "));
      std::cout << "\n";
    }

    // add all the dependencies, atlev is a set so it doesn't matter if there
    // are repeats
    curr->atlev.insert(rp.depends[dep].begin(), rp.depends[dep].end());
    RestSet lopon = rp.loopons[dep];
    // make the child if it doesn't exist

    if (curr->children.find(lopon) == curr->children.end()) {
      curr->children[lopon] = new MultiRestPlan(dep + 1, labelled);
    }
    // progress to child
    curr = curr->children[lopon];
  }
  // add counter, and atlevs
  curr->atlev.insert(rp.depends[dep].begin(), rp.depends[dep].end());
  if (curr->counters.find(rp.loopons[dep]) == curr->counters.end()) {
    curr->counters[rp.loopons[dep]] = npls;
  } else {
    std::cout << "DUPLICATE PLANS ENCOUNTERED" << std::endl;
  }
  npls = std::max(npls, 1 + rp.id());
}

std::pair<std::string, std::string> MultiRestPlan::to_code(
    int mpdepth, int stdepth, int mode, bool noCM)
{
  std::ostringstream oss1, oss2;
  if (noCM && false) {
    std::set<std::string> variablesInstanced;
    append_to_stream_nocm(oss1, oss2, variablesInstanced);
  } else {
    append_to_stream(oss1, oss2, mpdepth, stdepth, mode);
  }
  return { oss1.str(), oss2.str() };
}

int MultiRestPlan::maxdepth(const std::map<RestSet, MultiRestPlan *> &pChildren) const
{
  if (pChildren.empty()) return 0;

  return maxdepth(pChildren.begin()->second->children) + 1;
}

void MultiRestPlan::append_initializers(std::ostream &oss) {
  oss << "const int n_counters = " << npls << ";\n";
  oss << "std::vector<std::vector<uint64_t>> "
    "global_counters(omp_get_max_threads());\n";
  oss << "std::vector<std::vector<vidType>> "
    << "mMatchingVec(omp_get_max_threads(), "
    << "std::vector<vidType>(" << maxdepth(children) + 1 << ",0)"
    << "); //ignore\n";
  oss << "std::vector<Callback> callbackVec"
    << "(omp_get_max_threads(), Callback(print));"
    << " //ignore\n";
  oss << "#pragma omp parallel\n{\n";
  oss << "  const int tid = omp_get_thread_num(); //ignore\n";
  oss << "  std::vector<uint64_t> &counter = "
    "global_counters.at(omp_get_thread_num());\n";
  oss << "  counter.resize(n_counters, 0);\n";
}

void MultiRestPlan::append_first_loop(std::ostream &oss, std::string &tabdep) {
  int dynsize = 64;
  oss << tabdep <<  "#pragma omp for schedule(dynamic," << dynsize << ") nowait\n";
  oss << tabdep;
  children.begin()->first.append_iter_to_stream(oss, 0);
  oss << "{\n";
  oss << tabdep << "  auto &mMatching = mMatchingVec[tid]; // ignore\n";
}

void MultiRestPlan::instanciateVariables(std::ostream &oss, std::string &tabdep, const RestSet &restSet) {

}

void MultiRestPlan::append_to_stream_nocm_rec(
    std::ostream &oss_head,
    std::ostream &oss,
    std::set<std::string> &instanciated) {

  std::cout << "depth " << depth << std::endl;
  std::ostringstream tabdepthoss;
  tabdepthoss << "  ";
  for (int i = 0; i < depth; ++i) tabdepthoss << "  ";
  std::string tabdep = tabdepthoss.str();
  if (children.empty()) {
    std::string lastVarname = atlev.rbegin()->varname;
    append_callback_func(oss, tabdep, lastVarname);
    return;
  }
  std::cout << "calling child with depth " << depth << std::endl;
  const auto &restSet = children.begin()->first;
  instanciateVariables(oss, tabdep, restSet);
  oss << tabdep;
  children.begin()->first.append_iter_to_stream(oss, 0);
  oss << "{\n";
  children.begin()->second->append_to_stream_nocm_rec(oss_head, oss, instanciated);
  oss << tabdep << "}\n";

}

void MultiRestPlan::append_to_stream_nocm(
    std::ostream &oss_head,
    std::ostream &oss,
    std::set<std::string> &instanciated) {
  append_initializers(oss);
  std::cout << "depth " << depth << std::endl;
  std::ostringstream tabdepthoss;
  tabdepthoss << "  ";
  for (int i = 0; i < depth; ++i) tabdepthoss << "  ";
  std::string tabdep = tabdepthoss.str();
  std::cout << "calling child with depth " << depth << std::endl;
  append_first_loop(oss, tabdep);
  oss << tabdep << "  mMatching[" << depth << "] = v" << depth << "; //ignore \n";
  auto nextChildIter = children.begin();
  nextChildIter->second->append_to_stream_nocm_rec(oss_head, oss, instanciated);
  oss << tabdep << "}\n";
  append_deintializers(oss);
}

void MultiRestPlan::append_callback_func(std::ostream &oss, std::string &tabdep,
    std::string &lastVarname) {
  oss << tabdep << "for(vidType v" << depth << ":" << lastVarname << ") {\n";
  oss << tabdep << "  mMatching.at(" << depth << ")"
    << " = v" << depth << "; //ignore\n";
  oss << tabdep << "  callbackVec[tid](mMatching); //counter\n";
  oss << tabdep << "}\n";
}

void MultiRestPlan::append_deintializers(std::ostream &oss) {
  oss << "}//217\n";// pragma omp parallel
  oss << "const uint64_t data_complexity = " << data_complexity() << ";\n";
  oss << "const uint64_t time_complexity = " << time_complexity() << ";\n";
  oss << "Callback result = std::accumulate(callbackVec.begin(),"
    << " callbackVec.end(), Callback(false));\n";
  oss << "return result;\n";
}

void MultiRestPlan::append_to_stream(
    std::ostream &oss_head,
    std::ostream &oss,
    int fordepth,
    int start,
    int mode)
{
  if (depth == 0) {
    append_initializers(oss);
  }
  // print all the sets we use
  std::ostringstream tabdepthoss;
  tabdepthoss << "  ";
  for (int i = 0; i < depth; ++i) tabdepthoss << "  ";
  std::string tabdep = tabdepthoss.str();
  std::string lastVarname;
  for (const auto &restset : atlev) {
    oss << tabdep;
    restset.append_calc_to_stream(oss, -1, atlev, tabdep);
    lastVarname = restset.varname;
  }
  for (auto child : children) {
    if (depth == 0) {
      int dynsize = 64;
      oss << tabdep << "#pragma omp for schedule(dynamic," << dynsize
        << ") nowait\n";
    }
    oss << tabdep;
    child.first.append_iter_to_stream(oss, depth);
    oss << " {\n";
    if (depth == 0)
      oss << tabdep << "  auto &mMatching = mMatchingVec.at(tid); //ignore\n";
    oss << tabdep << "  mMatching.at(" << depth << ")"
      << " = v" << depth << "; //ignore\n";
    child.second->append_to_stream(
        oss_head, oss, fordepth, start, mode);
    oss << tabdep << "}\n";
  }
  if (children.empty()) {
    append_callback_func(oss, tabdep, lastVarname);
  }
  if (depth == 0) {
    append_deintializers(oss);
  }
}

double MultiRestPlan::data_complexity()
{
  std::vector<int> temp;
  std::vector<int> rest;
  rest.push_back(-1);
  temp.push_back(0);
  return rec_data_complexity(GLOBAL_VERTEX_COUNT, rest, temp);
}

double MultiRestPlan::rec_data_complexity(
    double loopedness,
    std::vector<int> &restrictions,
    std::vector<int> &multcounts)
{
  // std::cout<<"loopy"<<loopedness<<std::endl;
  double res = 0;
  double comploop = loopedness;
  {
    int temp = depth;
    while (temp != -1) {
      ++multcounts.at(temp);
      temp = restrictions[temp];
    }
  }
  for (int i : multcounts) { comploop /= i; }
  for (RestSet s : atlev) {
    // these would be time complexity if we are calculating time complexity
    res += comploop * s.data_complexity_ignoring_restrictions();
  }
  // check for counts not calculated at the level
  for (auto m : counters) {
    if (atlev.count(m.first) == 0) {
      res += comploop * m.first.data_complexity_ignoring_restrictions();
    }
  }
  // pair
  for (auto m : children) {
    if (atlev.count(m.first) == 0) {
      res += comploop * m.first.data_complexity_ignoring_restrictions();
    }
    restrictions.push_back(m.first.restriction());
    multcounts.push_back(0);

    res += m.second->rec_data_complexity(
        loopedness * m.first.data_complexity_ignoring_restrictions(),
        restrictions,
        multcounts);

    multcounts.pop_back();
    restrictions.pop_back();
  }

  {
    int temp = restrictions.back();
    while (temp != -1) {
      --multcounts.at(temp);
      temp = restrictions[temp];
    }
  }
  return res;
}

double MultiRestPlan::time_complexity()
{
  std::vector<int> temp;
  std::vector<int> rest;
  rest.push_back(-1);
  temp.push_back(0);
  return rec_time_complexity(GLOBAL_VERTEX_COUNT, rest, temp);
}

size_t MultiRestPlan::size() const { return plans.size(); }

RestPlan MultiRestPlan::getPlan(size_t idx) const { return plans[idx]; }

double MultiRestPlan::rec_time_complexity(
    double loopedness,
    std::vector<int> &restrictions,
    std::vector<int> &multcounts)
{
  // std::cout<<"loopy"<<loopedness<<std::endl;
  double res = 0;
  double comploop = loopedness;
  // add the multiplicity calculation.
  {
    int temp = depth;
    while (temp != -1) {
      ++multcounts.at(temp);
      temp = restrictions[temp];
    }
  }
  // this is how many would be achieved
  for (int i : multcounts) { comploop /= i; }
  for (RestSet s : atlev) {
    // these would be time complexity if we are calculating time complexity
    res += comploop * s.time_complexity_ignoring_restrictions(false, atlev);
  }

  for (auto d : counters) {
    // for(auto m:d.second){
    if (atlev.count(d.first) == 0) {
      res +=
        comploop * d.first.time_complexity_ignoring_restrictions(true, atlev);
    }
    //}
  }
  // pair
  for (auto d : children) {

    if (atlev.count(d.first) == 0 && counters.count(d.first) == 0) {
      res +=
        comploop * d.first.time_complexity_ignoring_restrictions(true, atlev);
    }
    // update restrictions to be processed by the child
    restrictions.push_back(d.first.restriction());
    multcounts.push_back(0);

    res += d.second->rec_time_complexity(
        loopedness * d.first.data_complexity_ignoring_restrictions(),
        restrictions,
        multcounts);

    multcounts.pop_back();
    restrictions.pop_back();
    // undo the restrictions
  }

  {
    int temp = restrictions.back();
    while (temp != -1) {
      --multcounts.at(temp);
      temp = restrictions[temp];
    }
  }
  // return restrictions to how it was originally
  return res;
}

// template<class T>
// void print_vector(const std::vector<T> &v)
// {
//   std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
//   //std::cout << std::endl;
// }
//
// void print_vector(std::vector<T> &v)
// {
//   std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
//   //std::cout << std::endl;
// }

template<class T>
void print_set(std::set<T> &v)
{
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(std::cout, " "));
  std::cout << std::endl;
}

std::ostream & operator<<(std::ostream & os, const RestSet& rs)
{
  os << "RestSet: varname:" << rs.varname << " ins: (";
  std::copy(rs.ins.begin(), rs.ins.end(), std::ostream_iterator<int>(os, " "));
  os << ") out: (";
  std::copy(rs.out.begin(), rs.out.end(), std::ostream_iterator<int>(os, " "));
  os << ") restriction: (";
  for(auto i : rs.res_chain)
    os << i << ",";
  os << ") label: " << rs.label;
  return os;
}

std::ostream & operator<<(std::ostream & os, RestPlan & rp)
{
  std::cout << "rest: (";
  //print_vector<int>(rp.rest);
  std::cout << ") \n loopns: \n";
  for(const auto &rs : rp.loopons)
    std::cout << rs << "\n";
  //std::cout << "depends:\n";
  int i=0;
  for(const auto &d : rp.depends) {
    std::cout << "depends[" << i++ << "]: ";
    for(const auto &rs : d) {
      std::cout << "[" << rs << "]";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  return os;
}

std::ostream & operator<<(std::ostream & os, const Graphlet& g) {
  os << "Graphlet(" << g.vertices << ") ";
  os << "[ ";
  for (int i = 0; i < g.adjacency.size(); i++) {
    os << i << " -> { ";
    for (auto neigh : g.adjacency[i])
      os << neigh << " ";
    if (i+1 == g.adjacency.size())
      os << "}";
    else
      os << "}, ";
  }
  os << " ]" << std::endl;
  return os;
}

std::ostream & operator<<(std::ostream & os, const std::vector<ExecutionPlan>& plans) {
  os << "std::vector<ExecutionPlan>(" << plans.size() << ") \n";
  for (const auto &plan : plans) {
    os << "vertices : " << plan.vertices << "\n";
    os << "restrictions : {";
    for (const auto res : plan.restrictions)
      os << res << " ";
    os << "}\n";
    os << "depends (" << plan.depend.size() << ") : ";
    for (const auto& [f, s] : plan.depend) {
      os << "{ ";
      for (auto e : f) os << e << " ";
      os << "}, ";
      os << "{ ";
      for (auto e : s) os << e << " ";
      os << "} \n";
    }
    os << "\n";
  }
  os << std::endl;
  return os;
}

std::ostream & operator<<(std::ostream & os, const MultiRestPlan &mrp) {
  for (const auto &rs : mrp.atlev) {
    if (rs.ins.empty() && rs.out.empty()) {
      if (rs.label != -1)
        os << "g.L(" << rs.label << ")\n";
      continue;
    }
    if (rs.ins.size() == 1 && rs.out.empty()) {
      if (rs.res_chain[rs.depth] == -1) {
        os << rs.varname << " = g.N(v" << rs.ins[0] << "," << rs.label << ");\n";
      }
      continue;
    }
    RestSet par = rs.parent();
    if (par.ins.empty()) {
      if (par.out.empty()) {
        assert(false);
        continue;
      }
      // recursive difference set
    } else if (rs.out.size() != par.out.size()) {
      os << rs.varname << " = DS(" << par.varname << ", y" << rs.out[rs.out.size()-1]
                                            << "l" << rs.label;
      if (rs.res_chain[rs.depth] != -1) os << ", v" << rs.res_chain[rs.depth];
      os << ")\n";
    } else if (rs.ins.size() != par.ins.size()) {
      os << rs.varname << " = IS(" << par.varname << ", y" << rs.ins[rs.ins.size()-1]
      << "l" << rs.label;
      if (rs.res_chain[rs.depth] != -1) os << ", v" << rs.res_chain[rs.depth];
      os << ")\n";
    }
  }
  for (const auto& child : mrp.children) {
    os << "for (vidType v<num> : " << child.first.varname << ")\n";
    os << *(child.second);

  }
  return os;
}
