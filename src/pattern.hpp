#pragma once

#include "common.h"
#include "def.h"

#include <algorithm>
#include <boost/range/iterator_range_core.hpp>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <fmt/core.h>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <set>
#include <type_traits>
#include <queue>
#include <memory>
#include <vector>

#include <iterator>
#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/lookup_edge.hpp>
#include <boost/graph/graphviz.hpp>

#include <bliss/digraph.hh>

#include <vf3lib/ARGraph.hpp>
#include <vf3lib/VFLib.h>

#include "generateInstructions.hpp"

MatchingSchedule::MatchingSchedule(uint32_t pSize) { // NOLINT(misc-definitions-in-headers)
  mVertexToMatchingOrder.resize(pSize, std::numeric_limits<uint32_t>::max());
  mMatchingOrderToVertex.resize(pSize, std::numeric_limits<uint32_t>::max());
}

class Core;
class Vf3Pattern : public vflib::ARGLoader<uint32_t, uint32_t> {
  private:

    using Labels = Type::Labels;
    using LabelsVector = Type::LabelsVector;
    using Edges = Type::Edges;
    using LabelledEdges = Type::LabelledEdges;

    LabelsVector mLabels;
    LabelledEdges mEdges;
    LabelledEdges mIncomingEdges;

    MatchingSchedule mMatchingSchedule;

    friend class Core;

    bool mIsClique{false};
    Permutation mCliqueCanonPermutation{};
    CanonicalShorthand mCliqueShorthand{};

    BoostGraph mBoostGraph;

    std::shared_ptr<bliss::Digraph> mBlissGraph; // canonical_form_
    Permutation mCanonicalPermutation;
    Permutation mCanonicalPermutationInverse;

#ifdef D_CUSTOM_FREQUENCY
    instruction::Instructions mInstructions;
#endif

    template <typename T>
    void edgeArrow(T &ss, uint i, uint j) const {
      auto outI = mEdges.find(i);
      auto outJ = mEdges.find(j);
      bool IJ = outI->second.find(j) != outI->second.end();
      bool JI = outJ->second.find(i) != outJ->second.end();
      if (IJ && JI) {
        auto IJLabel = outI->second.find(j)->second;
        auto JILabel = outJ->second.find(i)->second;
        ss << " <"<<JILabel<<"-"<<IJLabel<<"> ";
      } else if (IJ) {
        auto IJLabel = outI->second.find(j)->second;
        ss << " ---"<<IJLabel<<"> ";
      } else if (JI) {
        auto JILabel = outJ->second.find(i)->second;
        ss << " <"<<JILabel<<"--- ";
      } else {
        ss << "     ";
      }
    }

    std::string stringGraph(uint i, uint j) const {
      std::stringstream ss;
      ss << mLabels[i];
      edgeArrow(ss, i, j);
      ss << mLabels[j];
      return ss.str();
    }

    std::ostream &printGraph(std::ostream &os, uint i, uint j) const {
      os << stringGraph(i, j);
      return os;
    }

    std::tuple<uint, uint, uint> reorder(uint i, uint j, uint k) const {
      // switch i, j and k if an edge is absent between two vertices
      // then make them i and k
      auto outI = mEdges.find(i);
      auto outJ = mEdges.find(j);
      auto outK = mEdges.find(k);
      bool IJ = outI->second.find(j) != outI->second.end() || outJ->second.find(i) != outJ->second.end();
      bool JK = outJ->second.find(k) != outJ->second.end() || outK->second.find(j) != outK->second.end();
      bool KI = outK->second.find(i) != outK->second.end() || outI->second.find(k) != outI->second.end();
      if (!IJ && JK && KI) {
        std::swap(j, k);
        std::swap(outJ, outK);
      } else if (!JK && KI && IJ) {
        std::swap(i, j);
        std::swap(outI, outJ);
      }
      return std::make_tuple(i, j, k);
    }

    std::string stringGraph(uint i, uint j, uint k) const {
      std::stringstream ss;
      std::tie(i, j, k) = reorder(i, j, k);
      edgeArrow(ss, k, i);
      ss << mLabels[i];
      edgeArrow(ss, i, j);
      ss << mLabels[j];
      edgeArrow(ss, j, k);
      ss << mLabels[k];
      edgeArrow(ss, k, i);
      return ss.str();
    }

    std::ostream &printGraph(std::ostream &os, uint i, uint j, uint k) const {
      os << stringGraph(i, j, k);
      return os;
    }

    void minimumPermutationRecurse(
      std::vector<unsigned int>& current_minimum, std::set<unsigned int>& unused,
      std::vector<unsigned int>& partial, unsigned int n, const BoostGraph &graph) {
    unsigned int curr = partial.size();
    if (n == curr) {
      return;
    }
    std::set<unsigned int> unused_copy(unused.begin(), unused.end());
    for (unsigned int next : unused) {
      // test if leq
      // recurse deeper
      bool isless = false;
      bool ismore = false;
      // check label of current vertex, make it smaller
      unsigned int old = current_minimum[curr];
      int32_t old_label = boost::get(boost::vertex_name, graph, old);
      int32_t next_label = boost::get(boost::vertex_name, graph, next);
      if (old_label == next_label) {
        // compare edges
        for (unsigned int neighbor = 0; neighbor < curr; ++neighbor) {
          // lookup edge is very slow, would probably benefit from having a
          // clique-specific method for having these benched
          auto to_new_edge =
              boost::lookup_edge(partial[neighbor], next, graph);
          auto to_old_edge =
              boost::lookup_edge(current_minimum[neighbor], old, graph);
          //edge not existing is smaller
          if (to_old_edge.second != to_new_edge.second) {
            isless = to_old_edge.second;
            ismore = to_new_edge.second;
            break;
          }
          if (to_old_edge.second && to_new_edge.second) {
            /* TODO: see if creates bugs
            int32_t to_old =
                boost::get(boost::edge_name_t(), graph, to_old_edge.first);
            int32_t to_new =
                boost::get(boost::edge_name_t(), graph, to_new_edge.first);

            if (to_new != to_old) {
              isless = (to_new < to_old);
              ismore = !isless;
              break;
            }
            */
          }
          auto from_new_edge =
              boost::lookup_edge(next, partial[neighbor], graph);
          auto from_old_edge =
              boost::lookup_edge(old, current_minimum[neighbor], graph);
          if (from_old_edge.second != from_new_edge.second) {
            // existence is larger than non-existence
            isless = from_old_edge.second;
            ismore = from_new_edge.second;
            break;
          }
          if (from_old_edge.second && from_new_edge.second) {
            /* TODO: see if creates bugs
            int32_t from_new =
                boost::get(boost::edge_name_t(), graph, from_new_edge.first);
            int32_t from_old =
                boost::get(boost::edge_name_t(), graph, from_old_edge.first);
            if (from_new != from_old) {
              isless = (from_new < from_old);
              ismore = !isless;
              break;
            }
            */
          }
        }
      } else {
        isless = next_label < old_label;
        ismore = !isless;
      }
      // check labels/existence of edges from 0 to curr-1 (inclusive)
      // if less, set current_minimum to partial+unused
      if (isless) {
        size_t i;
        for (i = 0; i < curr; ++i) {
          current_minimum[i] = partial[i];
        }
        current_minimum[curr] = next;
        auto it = unused.begin();
        for (i = curr + 1; i < n; ++i) {
          if (*it == next)
            ++it;
          current_minimum[i] = *it;
          ++it;
        }
      }
      if (!ismore) {
        partial.push_back(next);
        unused_copy.erase(next);
        minimumPermutationRecurse(current_minimum, unused_copy, partial, n, graph);
        unused_copy.insert(next);
        partial.pop_back();
      }
    }
  }

  std::vector<uint32_t> MinimumPermutation(const BoostGraph &graph) {
    auto num = boost::num_vertices(graph);
    std::vector<uint> current_minimum(num);
    std::iota(current_minimum.begin(), current_minimum.end(), 0);
    std::set<uint> unused(current_minimum.begin(), current_minimum.end());
    std::vector<uint> partial;
    minimumPermutationRecurse(current_minimum, unused, partial, num, graph);
    return current_minimum;
  }

  CanonicalShorthand getCanonicalShorthand(const BoostGraph &graph,
      const std::vector<uint32_t> &permutation) {
    CanonicalShorthand result;
    for (uint i = 0; i < permutation.size(); ++i) {
      unsigned int curr = permutation[i];
      int32_t label = boost::get(boost::vertex_name, graph, curr);
      std::vector<std::pair<MaybeEdge, MaybeEdge>> edges;
      for (unsigned int neighbor = 0; neighbor < i; ++neighbor) {
        unsigned int other = permutation[neighbor];
        auto to_edge = boost::lookup_edge(other, curr, graph);
        auto from_edge = boost::lookup_edge(curr, other, graph);
        auto toMaybeEdge = std::make_pair(to_edge.second, 0);
        auto fromMaybeEdge = std::make_pair(from_edge.second, 0);
        //auto toMaybeEdge = std::make_pair(
        //    to_edge.second, to_edge.second ? boost::get(boost::edge_name_t(), graph, to_edge.first) : 0);
        //auto fromMaybeEdge = std::make_pair(
        //    from_edge.second, from_edge.second ? boost::get(boost::edge_name_t(), graph, from_edge.first) : 0);
        edges.emplace_back(toMaybeEdge, fromMaybeEdge);
      }
      result.emplace_back(label, edges);
    }
    return result;
  }

  public:
    Vf3Pattern(LabelsVector pLabels, std::vector<std::pair<uint32_t, uint32_t>> pEdges)
        : mLabels(pLabels)
    {
#ifdef D_CUSTOM_FREQUENCY
      assert(mInstructions.size() == 0);
#endif
      for (auto [src, dst] : pEdges) {
        mEdges[src].insert(std::make_pair(dst, kDefaultEdgeLabel));
        mIncomingEdges[dst].insert(std::make_pair(src, kDefaultEdgeLabel));
      }
      mMatchingSchedule = GraphHelpers::generateSchedule(mLabels.size(), mEdges, mLabels, mIncomingEdges);

#ifdef D_CUSTOM_FREQUENCY
      mInstructions = GraphHelpers::generateInstructions<instruction::Instruction>(mLabels, mEdges);
#endif

      for (auto labels : mLabels) {
        boost::add_vertex(labels, mBoostGraph);
      }
      for (auto [src, dst] : pEdges) {
        boost::add_edge(src, dst, kDefaultEdgeLabel, mBoostGraph);
      }
      mIsClique = GraphHelpers::isClique(mEdges, mLabels.size());
      if (mIsClique) {
        mCliqueCanonPermutation = MinimumPermutation(mBoostGraph);
        mCliqueShorthand = getCanonicalShorthand(mBoostGraph, mCliqueCanonPermutation);
      }
      std::tie(mBlissGraph, mCanonicalPermutation, mCanonicalPermutationInverse) = GraphHelpers::canonicalForm(mBoostGraph);
    }

    Vf3Pattern(LabelsVector pLabels, std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> pLabelledEdges)
        : mLabels(pLabels)
    {
#ifdef D_CUSTOM_FREQUENCY
      assert(mInstructions.size() == 0);
#endif
      for (auto [src, dst, label] : pLabelledEdges) {
        mEdges[src].insert(std::make_pair(dst, label));
        mIncomingEdges[dst].insert(std::make_pair(src, label));
      }
      mMatchingSchedule = GraphHelpers::generateSchedule(mLabels.size(), mEdges, mLabels, mIncomingEdges);

#ifdef D_CUSTOM_FREQUENCY
      mInstructions = GraphHelpers::generateInstructions<instruction::Instruction>(mLabels, mEdges);
#endif

      for (auto labels : mLabels) {
        boost::add_vertex(labels, mBoostGraph);
      }
      for (auto [src, dst, label] : pLabelledEdges) {
        boost::add_edge(src, dst, label, mBoostGraph);
      }
      mIsClique = GraphHelpers::isClique(mEdges, mLabels.size());
      if (mIsClique) {
        mCliqueCanonPermutation = MinimumPermutation(mBoostGraph);
        mCliqueShorthand = getCanonicalShorthand(mBoostGraph, mCliqueCanonPermutation);
      }
      std::tie(mBlissGraph, mCanonicalPermutation, mCanonicalPermutationInverse) = GraphHelpers::canonicalForm(mBoostGraph);
    }

    explicit Vf3Pattern(const CanonicalShorthand &pShorthand) : mBoostGraph(GraphHelpers::boostFromShorthand(pShorthand)) {
#ifdef D_CUSTOM_FREQUENCY
      assert(mInstructions.size() == 0);
#endif
      mIsClique = GraphHelpers::isClique(mBoostGraph);
      if (mIsClique) {
        mCliqueCanonPermutation = MinimumPermutation(mBoostGraph);
        mCliqueShorthand = getCanonicalShorthand(mBoostGraph, mCliqueCanonPermutation);
      }

      mLabels = GraphHelpers::getLabels(mBoostGraph);
      mEdges = GraphHelpers::getEdges(mBoostGraph);
      mIncomingEdges = GraphHelpers::getIncomingEdges(mBoostGraph);

      mMatchingSchedule = GraphHelpers::generateSchedule(mLabels.size(), mEdges, mLabels, mIncomingEdges);

#ifdef D_CUSTOM_FREQUENCY
      mInstructions = GraphHelpers::generateInstructions<instruction::Instruction>(mLabels, mEdges);
#endif

      std::tie(mBlissGraph, mCanonicalPermutation, mCanonicalPermutationInverse) = GraphHelpers::canonicalForm(mBoostGraph);
    }

    Vf3Pattern(const BoostGraph &pBG) : mBoostGraph(pBG), mIsClique(GraphHelpers::isClique(pBG)) {
#ifdef D_CUSTOM_FREQUENCY
      assert(mInstructions.size() == 0);
#endif
      if (mIsClique) {
        mCliqueCanonPermutation = MinimumPermutation(mBoostGraph);
        mCliqueShorthand = getCanonicalShorthand(mBoostGraph, mCliqueCanonPermutation);
      }

      mLabels = GraphHelpers::getLabels(mBoostGraph);
      mEdges = GraphHelpers::getEdges(mBoostGraph);
      mIncomingEdges = GraphHelpers::getIncomingEdges(mBoostGraph);

      mMatchingSchedule = GraphHelpers::generateSchedule(mLabels.size(), mEdges, mLabels, mIncomingEdges);

#ifdef D_CUSTOM_FREQUENCY
      mInstructions = GraphHelpers::generateInstructions<instruction::Instruction>(mLabels, mEdges);
#endif

      std::tie(mBlissGraph, mCanonicalPermutation, mCanonicalPermutationInverse) = GraphHelpers::canonicalForm(mBoostGraph);
    }

    virtual ~Vf3Pattern() = default;

    std::string name() const {
      if (mLabels.size() == 2) {
        assert(mEdges.size() <= 2);
        return stringGraph(0, 1);
      } if (mLabels.size() == 3) {
        return stringGraph(0, 1, 2);
      }
      throw std::runtime_error("Not implemented");
    }

#ifdef D_CUSTOM_FREQUENCY
    instruction::Instructions instructions() const { return mInstructions; }
#endif

    virtual uint32_t NodeCount() const override { return mLabels.size(); }
    virtual uint32_t GetNodeAttr(vflib::nodeID_t node) const override {
      return mLabels.at(node);
    }
    virtual uint32_t OutEdgeCount(vflib::nodeID_t node) const override {
      return (mEdges.find(node) == mEdges.end()) ? 0 : mEdges.at(node).size();
    }
    virtual vflib::nodeID_t GetOutEdge(vflib::nodeID_t node, uint32_t index, uint32_t *pattr) const override {
      auto [neigh, edgeLabel] = GetOutEdgeWeight(node, index);
      pattr = nullptr;
      return neigh;
    }
    virtual std::pair<vflib::nodeID_t, uint32_t> GetOutEdgeWeight(vflib::nodeID_t node, uint32_t i) const override {
      assert(node >= 0);
      assert(node < mLabels.size());
      assert(i >= 0);
      assert(i < mEdges.at(node).size());
      uint neigh = std::next(mEdges.at(node).begin(), i)->first;
      uint edgeLabel = std::next(mEdges.at(node).begin(), i)->second;
      return std::make_pair(neigh, edgeLabel);
    }

    BoostGraph boostGraph() const {
      return mBoostGraph;
    }
    bool isClique() const {
      return mIsClique;
    }
    CanonicalShorthand canonicalShorthand() const {
      return mCliqueShorthand;
    }
    Permutation cliqueCanonPermutation() const {
      return mCliqueCanonPermutation;
    }
    uint32_t nodes() const { return mLabels.size(); }
    uint32_t label(uint32_t pVertex) const { return mLabels[pVertex]; }
    std::vector<uint32_t> labels() const { return mLabels; }
    uint32_t matchOrderToVertex(uint32_t pOrder) const { return mMatchingSchedule.mMatchingOrderToVertex[pOrder]; }
    uint32_t vertexToMatchOrder(uint32_t pVertex) const { return mMatchingSchedule.mVertexToMatchingOrder[pVertex]; }

    bool operator<(const Vf3Pattern &pOther) const {
      if (boost::num_vertices(mBoostGraph) != boost::num_vertices(pOther.mBoostGraph)) {
        return boost::num_vertices(mBoostGraph) < boost::num_vertices(pOther.mBoostGraph);
      }
      auto aHash = mBlissGraph->get_hash();
      auto bHash = pOther.mBlissGraph->get_hash();
      if (aHash != bHash) {
        return aHash < bHash;
      }
      auto blisscmp = mBlissGraph->cmp(*(pOther.mBlissGraph));
      return blisscmp < 0;
    }

    bool operator==(const Vf3Pattern &pOther) const {
      if (boost::num_vertices(mBoostGraph) != boost::num_vertices(pOther.mBoostGraph)) {
        return false;
      }
      auto aHash = mBlissGraph->get_hash();
      auto bHash = pOther.mBlissGraph->get_hash();
      if (aHash != bHash) {
        return false;
      }
      auto blisscmp = mBlissGraph->cmp(*(pOther.mBlissGraph));
      return blisscmp == 0;
    }

    void visualize(std::ostream& os) const {
      std::string fileName = fmt::format("figures/patterns/p{}.dot", mBlissGraph->get_hash());
      GraphHelpers::dot(mBoostGraph, fileName);
      if (mLabels.size() == 2) {
        os << stringGraph(0, 1);
        return;
      }
      if (mLabels.size() == 3) {
        os << stringGraph(0, 1, 2);
        return;
      }
      else if (isClique()) {
        os << "Pattern(" << mLabels.size() << " nodes, " << mEdges.size() << " edges, clique ";
        std::copy(mLabels.begin(), mLabels.end(), std::ostream_iterator<uint32_t>(os, " "));
        os << fmt::format(" ) [{}] ", fileName);
      }
      else {
        os << fmt::format("Pattern( {} nodes, {} edges) [{}] \n", mLabels.size(), mEdges.size(), fileName);
        for (uint32_t src = 0; src < mLabels.size(); src++) {
          auto dstsIter = mEdges.find(src);
          if (dstsIter == mEdges.end()) {
            os << src << " (" << mLabels[src] << ") -> " << std::endl;
            continue;
          }
          os << src << " (" << mLabels[src] << ") -> ";
          for (auto [dst, label] : mEdges.at(src)) {
            os << dst << " [" << label << "] ";
          }
          os << std::endl;
        }
      }
    }

    friend std::ostream& operator<<(std::ostream& os, const Vf3Pattern &patt) {
      // print v id label
      // print e src dst label
      size_t numEdges = 0;
      for (auto [src, dstsTuples] : patt.mEdges) {
        numEdges += dstsTuples.size();
      } 
      os << "p " << patt.mLabels.size() << " " << numEdges << "\n";
      for (size_t i = 0; i < patt.mLabels.size(); i++) {
        os << "v " << i << " " << patt.mLabels[i] << "\n";
      }
      for (auto [src, dstsTuples] : patt.mEdges) {
        for (auto [dst, label] : dstsTuples) {
          os << "e " << src << " " << dst << " " << label << "\n";
        }
      }
      return os;
    }

  public:
    Permutation getCanonicalPermutation() const {
      return mCanonicalPermutation;
    }
    Permutation getCanonicalPermutationInverse() const {
      return mCanonicalPermutationInverse;
    }

};

template <typename Graph>
std::ostream& operator<<(std::ostream& os, const std::vector<Graph> &graph) {
  os << graph.size() << "\n";
  for (auto &g : graph) {
    os << g;
  }
  return os;
}

template <typename Graph>
std::ostream& operator<<(std::ostream& os, const std::set<Graph> &graph) {
  os << graph.size() << "\n";
  for (auto &g : graph) {
    os << g;
  }
  return os;
}

namespace GraphHelpers {

  using Labels = Type::Labels;
  using LabelsVector = Type::LabelsVector;
  using Edges = Type::Edges;
  using LabelledEdges = Type::LabelledEdges;

   std::unordered_map<uint32_t, uint32_t> getCandidatesLabelFreq( // NOLINT(misc-definitions-in-headers)
      std::vector<uint32_t> pCandidates, const LabelsVector &pLabels) {
    std::unordered_map<uint32_t, uint32_t> labelFreq;
    for (auto vid : pCandidates) {
      labelFreq[pLabels[vid]]++;
    }
    return labelFreq;
  }

  std::vector<std::pair<std::uint32_t, std::pair<uint32_t, uint32_t>>> getDegrees( // NOLINT(misc-definitions-in-headers)
      std::vector<uint32_t> pCandidates, Edges &pEdges,
      const std::unordered_map<uint32_t, uint32_t> &pLabelFreq,
      Edges &pIncomingEdges,
      const LabelsVector &pLabels) {
    std::vector<std::pair<std::uint32_t, std::pair<uint32_t, uint32_t>>> degrees;
    degrees.reserve(pCandidates.size());
    for (auto vid : pCandidates) {
      auto label = pLabels[vid];
      uint32_t freq = pLabelFreq.at(label);
      auto degreeFreq = std::make_pair(pEdges[vid].size() + pIncomingEdges[vid].size(), freq);
      degrees.emplace_back(std::make_pair(vid, degreeFreq));
    }
    return degrees;
  }

  uint32_t pickCandidateAndExtend(std::vector<uint32_t> &pCandidates, const LabelsVector &pLabels, Edges &pEdges, Edges &pIncomingEdges, MatchingSchedule &pMatchingSchedule) { // NOLINT(misc-definitions-in-headers)
    std::unordered_map<uint32_t, uint32_t> labelFreq = getCandidatesLabelFreq(pCandidates, pLabels);

    std::vector<std::pair<uint32_t, std::pair<uint32_t, uint32_t>>> degrees = getDegrees(pCandidates, pEdges, labelFreq, pIncomingEdges, pLabels);

    std::sort(degrees.begin(), degrees.end(), [](const auto &lhs, const auto &rhs) {
        auto [lhsDegree, lhsFreq] = lhs.second;
        auto [rhsDegree, rhsFreq] = rhs.second;
        return (lhsDegree != rhsDegree) ? lhsDegree < rhsDegree : lhsFreq > rhsFreq;
        });

    uint32_t retId = degrees.back().first;


    std::copy_if(pEdges.at(retId).begin(), pEdges.at(retId).end(), std::back_inserter(pCandidates), [&pCandidates, vertVec = pMatchingSchedule.mMatchingOrderToVertex] (uint32_t pId) {
        return std::find(pCandidates.begin(), pCandidates.end(), pId) == pCandidates.end() &&
              std::find(vertVec.begin(), vertVec.end(), pId) == vertVec.end();
        });

    std::copy_if(pIncomingEdges.at(retId).begin(), pIncomingEdges.at(retId).end(), std::back_inserter(pCandidates), [&pCandidates, vertVec = pMatchingSchedule.mMatchingOrderToVertex] (uint32_t pId) {
        return std::find(pCandidates.begin(), pCandidates.end(), pId) == pCandidates.end() &&
              std::find(vertVec.begin(), vertVec.end(), pId) == vertVec.end();
        });

    pCandidates.erase(std::remove(pCandidates.begin(), pCandidates.end(), retId), pCandidates.end());

    return retId;
  }

  MatchingSchedule generateSchedule(size_t pNumVertices, LabelledEdges &pEdges, // NOLINT(misc-definitions-in-headers)
      const LabelsVector &pLabels, LabelledEdges &pIncomingEdges) {
    MatchingSchedule schedule(pNumVertices);
    std::vector<std::pair<uint32_t, uint32_t>> degrees(pNumVertices);
    for (uint32_t vid = 0; vid < pNumVertices; ++vid) {
      degrees[vid] = std::make_pair(vid, pEdges[vid].size());
    }
    std::sort(degrees.begin(), degrees.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.second > rhs.second;
        });
    assert(degrees[0].second != 0);
    //assert(std::all_of(degrees.begin(), degrees.end(), [](const auto &p) {
    //      return p.second != 0;
    //      }));
    Edges edges;
    Edges incomingEdges;
    for (auto [src, dstsLabels] : pEdges) {
      for (auto [dst, edgeLabel] : dstsLabels) {
        edges[src].insert(dst);
      }
    }
    for (auto [src, dstsLabels] : pIncomingEdges) {
      for (auto [dst, edgeLabel] : dstsLabels) {
        incomingEdges[src].insert(dst);
      }
    }
    std::vector<uint32_t> candidates = { degrees[0].first };
    for (uint32_t i = 0; i < pNumVertices; ++i) {
      uint32_t vid = pickCandidateAndExtend(candidates, pLabels, edges, incomingEdges, schedule);
      schedule.mMatchingOrderToVertex[i] = vid;
      schedule.mVertexToMatchingOrder[vid] = i;
    }
    return schedule;
  }

  std::pair<BoostGraph, BoostGraph> convertToBoost(const BoostGraph &originalGraph, Permutation canonicalForm) {
    BoostGraph boostGraph, shortGraph;

    std::vector<std::pair<uint32_t, uint32_t>> vertexLabels;
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> edges;
    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> shortEdges;
    uint max_vertex_color = 0;

    for (size_t i = 0; i < boost::num_vertices(originalGraph); ++i) {
      uint32_t label = boost::get(boost::vertex_name, originalGraph, i);
      if (label != kMarkedLabel) {
        max_vertex_color = std::max(max_vertex_color, canonicalForm[i]);
      }
      vertexLabels.emplace_back(i, label);
    }
    uint32_t edge_num = boost::num_vertices(originalGraph);
    for (auto edge : boost::make_iterator_range(boost::edges(originalGraph))) {
      uint32_t edge_label = boost::get(boost::edge_name, originalGraph, edge);
      uint32_t blissEdgeLabel = edge_label + 1 + max_vertex_color;
      uint32_t source = boost::source(edge, originalGraph);
      uint32_t target = boost::target(edge, originalGraph);
      vertexLabels.emplace_back(edge_num, blissEdgeLabel);
      shortEdges[canonicalForm[source]].insert(canonicalForm[target]);
      edges[canonicalForm[source]].insert(canonicalForm[edge_num]);
      edges[canonicalForm[edge_num]].insert(canonicalForm[target]);
      edge_num++;
    }

    // permute the vertex labels with the canonical form
    std::sort(vertexLabels.begin(), vertexLabels.end(), [&canonicalForm](const auto &lhs, const auto &rhs) {
        return canonicalForm[lhs.first] < canonicalForm[rhs.first];
        });

    for (int i = 0; i < vertexLabels.size(); i++) {
      boost::add_vertex(vertexLabels[i].second, boostGraph);
      if (i < boost::num_vertices(originalGraph)) {
        boost::add_vertex(vertexLabels[i].second, shortGraph);
      }
    }

    for (auto [src, dsts] : edges) {
      for (auto dst : dsts) {
        boost::add_edge(src, dst, boostGraph);
      }
    }

    for (auto [src, dsts] : shortEdges) {
      for (auto dst : dsts) {
      boost::add_edge(src, dst, shortGraph);
      }
    }

    return std::make_pair(boostGraph, shortGraph);
  }

  std::tuple<std::shared_ptr<bliss::Digraph>, Permutation, Permutation>
    canonicalForm(const BoostGraph& pBG) { //NOLINT(misc-definitions-in-headers)
    bliss::Digraph graphv0;
    uint max_vertex_color = 0;

    for (auto vert : boost::make_iterator_range(boost::vertices(pBG))) {
      uint32_t label = boost::get(boost::vertex_name, pBG, vert);
      uint32_t vertex = graphv0.add_vertex(label);
      if (label != kMarkedLabel) {
        max_vertex_color = std::max(max_vertex_color, label);
      }
    }
    bliss::Stats stats;
#if 0
    const unsigned int* canon = graphv0.canonical_form(stats);
    Permutation canonical_permutation = Permutation(canon,
        canon + boost::num_vertices(pBG));
#else  // NOTE: needed for using double edges instead of uni edges
    for (auto edge : boost::make_iterator_range(boost::edges(pBG))) {
      uint32_t edge_label = boost::get(boost::edge_name, pBG, edge);
      uint32_t blissEdgeLabel = edge_label + 1 + max_vertex_color;
      uint32_t edge_num = graphv0.add_vertex(blissEdgeLabel);
      uint32_t source = boost::source(edge, pBG);
      uint32_t target = boost::target(edge, pBG);
      graphv0.add_edge(source, edge_num);
      graphv0.add_edge(edge_num, target);
    }
    const unsigned int* canon = graphv0.canonical_form(stats);
    std::vector<unsigned int> canonEdges(canon, canon + boost::num_vertices(pBG) + boost::num_edges(pBG));
    Permutation canonical_permutation = Permutation(canon,
        canon + boost::num_vertices(pBG) + boost::num_edges(pBG));
#endif
    Permutation canonical_permutation_inverse = Permutation(canonical_permutation.size());
    std::shared_ptr<bliss::Digraph> canonical_form = std::unique_ptr<bliss::Digraph>(graphv0.permute(canon));
    for (unsigned int j = 0; j < canonical_permutation.size(); ++j) {
      canonical_permutation_inverse[canonical_permutation[j]] = j;
    }
    return std::make_tuple(canonical_form, canonical_permutation, canonical_permutation_inverse);
  }

  std::tuple<std::shared_ptr<bliss::Digraph>, Permutation, Permutation>
    canonicalFormShort(const BoostGraph& pBG) { //NOLINT(misc-definitions-in-headers)
      bliss::Digraph graphv0;
      uint max_vertex_color = 0;

      for (auto vert : boost::make_iterator_range(boost::vertices(pBG))) {
        uint32_t label = boost::get(boost::vertex_name, pBG, vert);
        uint32_t vertex = graphv0.add_vertex(label);
        if (label != kMarkedLabel) {
          max_vertex_color = std::max(max_vertex_color, label);
        }
      }
      bliss::Stats stats;
      for (auto edge : boost::make_iterator_range(boost::edges(pBG))) {
        uint32_t source = boost::source(edge, pBG);
        uint32_t target = boost::target(edge, pBG);
        graphv0.add_edge(source, target);
      }
      const unsigned int* canon = graphv0.canonical_form(stats);
      Permutation canonical_permutation = Permutation(canon,
          canon + boost::num_vertices(pBG));

      Permutation canonical_permutation_inverse = Permutation(canonical_permutation.size());
      std::shared_ptr<bliss::Digraph> canonical_form = std::unique_ptr<bliss::Digraph>(graphv0.permute(canon));
      for (unsigned int j = 0; j < canonical_permutation.size(); ++j) {
        canonical_permutation_inverse[canonical_permutation[j]] = j;
      }
      return std::make_tuple(canonical_form, canonical_permutation, canonical_permutation_inverse);
    }

  bool isClique(const LabelledEdges &pLabelledEdges, const uint32_t pNumNodes) { // NOLINT(misc-definitions-in-headers)
    for (size_t i = 0; i < pNumNodes; ++i) {
      for (size_t j = i + 1; j < pNumNodes; ++j) {
        if (pLabelledEdges.at(i).find(j) == pLabelledEdges.at(i).end() &&
            pLabelledEdges.at(j).find(i) == pLabelledEdges.at(j).end()) {
          return false;
        }
      }
    }
    return true;
  }

  bool isClique(const BoostGraph &graph) { // NOLINT(misc-definitions-in-headers)
    for (VertexDescriptor i = 0; i < boost::num_vertices(graph); ++i) {
      for (VertexDescriptor j = i + 1; j < boost::num_vertices(graph); ++j) {
        if (!boost::edge(i, j, graph).second && !boost::edge(j, i, graph).second) {
          return false;
        }
      }
    }
    return true;
  }

  LabelsVector getLabels(const BoostGraph &pBG) { // NOLINT(misc-definitions-in-headers)
    LabelsVector labels(boost::num_vertices(pBG));
    for (VertexDescriptor i = 0; i < boost::num_vertices(pBG); ++i) {
      labels[i] = boost::get(boost::vertex_name, pBG, i);
    }
    return labels;
  }

  LabelledEdges getEdges(const BoostGraph &pBG) { // NOLINT(misc-definitions-in-headers)
    LabelledEdges edges;

    for (auto edge : boost::make_iterator_range(boost::edges(pBG))) {
      VertexDescriptor src = boost::source(edge, pBG);
      VertexDescriptor dst = boost::target(edge, pBG);
      auto weight = boost::get(boost::edge_name, pBG, edge);
      edges[src].insert(std::make_pair(dst, weight));
    }
    return edges;
  }

  LabelledEdges getIncomingEdges(const BoostGraph &pBG) { // NOLINT(misc-definitions-in-headers)
    LabelledEdges edges;
    for (auto edge : boost::make_iterator_range(boost::edges(pBG))) {
      VertexDescriptor src = boost::source(edge, pBG);
      VertexDescriptor dst = boost::target(edge, pBG);
      auto weight = boost::get(boost::edge_name, pBG, edge);
      edges[dst].insert(std::make_pair(src, weight));
    }
    return edges;
  }

  CanonicalShorthand getCanonicalShorthand(const BoostGraph &graph, // NOLINT(misc-definitions-in-headers)
      const Permutation &permutation) {
    CanonicalShorthand result;
    for (uint i = 0; i < permutation.size(); ++i) {
      unsigned int curr = permutation[i];
      int32_t label = boost::get(boost::vertex_name, graph, curr);
      std::vector<std::pair<MaybeEdge, MaybeEdge>> edges;
      for (unsigned int neighbor = 0; neighbor < i; ++neighbor) {
        unsigned int other = permutation[neighbor];
        auto to_edge = boost::lookup_edge(other, curr, graph);
        auto from_edge = boost::lookup_edge(curr, other, graph);
        auto toMaybeEdge = std::make_pair(to_edge.second,  0);
        auto fromMaybeEdge = std::make_pair(from_edge.second, 0);
        edges.emplace_back(toMaybeEdge, fromMaybeEdge);
      }
      result.emplace_back(label, edges);
    }
    return result;
  }

  Graphlet generateGraphlet(const std::vector<uint32_t> &pLabels, const LabelledEdges &pEdges) { // NOLINT(misc-definitions-in-headers)
    Graphlet graphlet(pLabels.size());
    for (size_t i = 0; i < pLabels.size(); ++i) {
      graphlet.labels[i] = pLabels[i];
    }
    for (auto [src, dstsLabels] : pEdges) {
      for (auto [dst, edgeLabel] : dstsLabels) {
        graphlet.add_edge(src, dst);
      }
    }
    return graphlet;
  }

  Graphlet generateGraphlet(const std::vector<uint32_t> &pLabels) { // NOLINT(misc-definitions-in-headers)
    Graphlet graphlet(pLabels.size());
    for (size_t i = 0; i < pLabels.size(); ++i) {
      graphlet.labels[i] = pLabels[i];
    }
    return graphlet;
  }

  MultiRestPlan generatePlan(Graphlet &pGraphlet) { // NOLINT(misc-definitions-in-headers)
    MultiRestPlan mrp(0, true);
    int total = 0;
    std::vector<std::vector<ExecutionPlan>> planss;
    MultiRed mr(pGraphlet);
    planss.push_back(mr.plans());
    total += mr.plans().size();
    for (auto &plans : planss) {
      double bestcomplex = std::numeric_limits<double>::infinity();
      int bindex = 0;
      for (int j = 0; j < plans.size(); ++j) {
        RestPlan test(plans[j], 0, false);
        double comple = test.time_complexity();
        if (comple < bestcomplex) {
          bestcomplex = comple;
          bindex = j;
        }
      }
      mrp.add_ex_plan(plans[bindex], 0, false);
    }
    return mrp;
  }

  template<typename Instruction>
  std::vector<Instruction> generateInstructions(const std::vector<uint32_t> &pLabels,
      const LabelledEdges &pEdges) {
      auto graphlet = GraphHelpers::generateGraphlet(pLabels, pEdges);
      auto mrp = GraphHelpers::generatePlan(graphlet);
    return generatorInstructions::generateInstructions<RestSet>(mrp);
  }

  BoostGraph boostFromShorthand(const CanonicalShorthand& pShorthand) { // NOLINT(misc-definitions-in-headers)
    BoostGraph bg;
    std::vector<size_t> vertices;
    size_t n = pShorthand.size();
    for (size_t i = 0; i <n; ++i) {
      vertices.emplace_back(boost::add_vertex(
            VertexProperty(pShorthand[i].first), bg));
    }
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < i; ++j) {
        if (pShorthand[i].second[j].first.first)
          boost::add_edge(vertices[j], vertices[i],
              pShorthand[i].second[j].first.second, bg);
        if (pShorthand[i].second[j].second.first)
          boost::add_edge(vertices[i], vertices[j],
              pShorthand[i].second[j].second.second, bg);
      }
    }
    return bg;
  }

  template<typename GraphClass>
  void dot(const GraphClass &graph, const std::string &filename) {
    std::ofstream out(filename);
    boost::write_graphviz(out, graph, boost::make_label_writer(boost::get(boost::vertex_name, graph)),
        boost::make_label_writer(boost::get(boost::edge_name, graph)));
  }

  template<typename GraphClass>
  void dotstream(const GraphClass &graph, std::ostream &out) {
    boost::write_graphviz(out, graph, boost::make_label_writer(boost::get(boost::vertex_name, graph)),
        boost::make_label_writer(boost::get(boost::edge_name, graph)));
  }

};
