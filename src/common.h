#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/lookup_edge.hpp>

#include <bliss/digraph.hh>

#include <fmt/core.h>
#include <stdint.h>

#include "def.h"

using Permutation = std::vector<uint32_t>;
using VertexProperty = boost::property<
  boost::vertex_name_t, uint32_t,
  boost::property<boost::vertex_index_t, uint>>;
using EdgeProperty = boost::property<boost::edge_name_t, uint32_t>;
using BoostGraph = boost::adjacency_list<
  boost::listS, boost::vecS, boost::bidirectionalS, VertexProperty, EdgeProperty>;
  //boost::property<boost::vertex_index_t, uint>,
  //boost::property<boost::edge_index_t, uint>>;

using VertexDescriptor = typename BoostGraph::vertex_descriptor;
using EdgeDescriptor = typename BoostGraph::edge_descriptor;

using NeighborInfo = std::pair<VertexDescriptor, int32_t>;
using MaybeEdge = std::pair<bool, int32_t>;
using CanonicalShorthand = std::vector<std::pair<int32_t, std::vector<std::pair<MaybeEdge, MaybeEdge>>>>;

inline constexpr uint32_t kMarkedLabel = std::numeric_limits<uint32_t>::max();
inline constexpr uint32_t kDefaultEdgeLabel = 0;

namespace Type {
  using Labels = uint32_t;
  using LabelsVector = std::vector<Labels>;
  using Edges = std::unordered_map<uint32_t, std::set<uint32_t>>;
  using LabelledEdges = std::unordered_map<uint32_t, std::map<uint32_t, uint32_t>>;
  using RemapEdges = std::unordered_map<uint32_t, uint32_t>;
  using RemapEdgesReverse = std::unordered_map<uint32_t, uint32_t>;
};

struct MatchingSchedule {
  std::vector<uint32_t> mVertexToMatchingOrder{};
  std::vector<uint32_t> mMatchingOrderToVertex{};
  MatchingSchedule() = default;
  MatchingSchedule(uint32_t pSize);
};

namespace GraphHelpers {

  MatchingSchedule generateSchedule(size_t pNumVertices, Type::LabelledEdges &pEdges,
      const Type::LabelsVector &pLabels, Type::LabelledEdges &pIncomingEdges);
  std::tuple<std::shared_ptr<bliss::Digraph>, Permutation, Permutation> canonicalForm(const BoostGraph &pBG);
  bool isClique(const BoostGraph &graph);
  bool isClique(const Type::LabelledEdges &pEdges, const uint32_t pNumNodes);
  Type::LabelsVector getLabels(const BoostGraph &pBG);
  Type::LabelledEdges getEdges(const BoostGraph &pBG);
  Type::LabelledEdges getIncomingEdges(const BoostGraph &pBG);
  CanonicalShorthand getCanonicalShorthand(const BoostGraph &pBG, const Permutation &pPermutation);
  BoostGraph boostFromShorthand(const CanonicalShorthand& pShorthand);
  Graphlet generateGraphlet(const std::vector<uint32_t> &pLabels, const Type::LabelledEdges &pEdges);
  Graphlet generateGraphlet(const std::vector<uint32_t> &pLabels);
  MultiRestPlan generatePlan(Graphlet &pGraphlet);
  template <typename Instruction>
  std::vector<Instruction> generateInstructions(const std::vector<uint32_t> &pLabels, const Type::LabelledEdges &pEdges);
  template<typename GraphClass>
  void dot(const GraphClass &pGraph, const std::string &pFilename);
  template<typename GraphClass>
  void dotStream(const GraphClass &pGraph, std::ostream &pStream);
};

template <typename... T>
void flushedPrint(fmt::format_string<T...> fmt, T &&... args) {
  fmt::print(fmt, std::forward<T>(args)...);
  fflush(stdout);
}
