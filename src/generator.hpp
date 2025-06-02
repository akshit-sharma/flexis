#pragma once

#include "common.h"

#include <cassert>
#include <chrono>
#include <cstdint>
#include <fmt/core.h>
#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <numeric>
#include <functional>
#include <queue>
#include <set>
#include <tuple>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/lookup_edge.hpp>

#include <bliss/digraph.hh>

namespace PatternGenerator {

  std::vector<std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>> edgeSetGenerator(uint x, uint y, uint elabel1, uint elabel2) { // NOLINT(misc-definitions-in-headers)
    return std::vector<std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>>{
      {{x, y, elabel1}}, // a -> b
      {{y, x, elabel2}}, // a <- b
      {{x, y, elabel1}, {y, x, elabel2}}, // a <-> b
    };
  }

  template <class PatternGraph>
  std::set<PatternGraph> dEdges(int pLabels) {
    std::set<PatternGraph> patterns;
    // sum n ... 1 = n * (n + 1) / 2
    int size = (pLabels * (pLabels + 1)) / 2;
    for (uint32_t i = 0; i < pLabels; i++) {
      for (uint32_t j = i; j < pLabels; j++) {
        patterns.insert(PatternGraph({i, j}, {{0, 1}, {1, 0}}));
      }
    }
    return patterns;
  }

  template <class PatternGraph>
  std::set<PatternGraph> edges(uint pLabels, uint pLabelEdges) {
    std::set<PatternGraph> patterns;
    for (uint e1 = 0; e1 < pLabelEdges; e1++) {
      for (uint e2 = e1; e2 < pLabelEdges; e2++) {
        auto edgeSet = edgeSetGenerator(0, 1, e1, e2);
        for (uint32_t i = 0; i < pLabels; i++) {
          for (uint32_t j = i; j < pLabels; j++) {
            for (auto edge : edgeSet) {
              patterns.insert(PatternGraph({i, j}, {edge}));
            }
          }
        }
      }
    }
    return patterns;
  }

//  using VertexProperty = boost::property<
//    boost::vertex_name_t, int32_t,
//    boost::property<boost::vertex_index_t, int>>;
//  using BoostGraph = boost::adjacency_list<
//    boost::listS, boost::vecS, boost::bidirectionalS,
//    VertexProperty>;

  using VertexDescriptor = typename BoostGraph::vertex_descriptor;
//  using Permutation = typename std::vector<unsigned int>;

  using VertexNameMap = typename boost::property_map<BoostGraph, boost::vertex_name_t>::type;

  std::tuple<VertexDescriptor, uint32_t,
   std::vector<std::pair<VertexDescriptor, int32_t>>,
   std::vector<std::pair<VertexDescriptor, int32_t>>> markRemovedVertex(BoostGraph &pBG, uint32_t pIndex) { //NOLINT(misc-definitions-in-headers)
      uint32_t index = pIndex;
      auto it = boost::vertices(pBG).first;
      while(pIndex--) {
        ++it;
      }
      VertexDescriptor marked_vertex = *it;
      uint32_t label = boost::get(boost::vertex_name, pBG, marked_vertex);
      std::vector<std::pair<VertexDescriptor, int32_t>> in_neighbours;
      for (auto edge : boost::make_iterator_range(boost::in_edges(marked_vertex, pBG))) {
        in_neighbours.push_back({boost::source(edge, pBG), boost::get(boost::edge_name, pBG, edge)});
      }
//      for (auto [i, i_end] = boost::in_edges(marked_vertex, pBG); i != i_end; ++i) {
//        in_neighbours.push_back({boost::source(*i, pBG), boost::get(boost::edge_name, pBG, *i)});
//      }
      std::vector<std::pair<VertexDescriptor, int32_t>> out_neighbours;
      for (auto edge : boost::make_iterator_range(boost::out_edges(marked_vertex, pBG))) {
        out_neighbours.push_back({boost::target(edge, pBG), boost::get(boost::edge_name, pBG, edge)});
      }
//      for (auto [i, i_end] = boost::out_edges(marked_vertex, pBG); i != i_end; ++i) {
//        out_neighbours.push_back({boost::target(*i, pBG), boost::get(boost::edge_name, pBG, *i)});
//      }

      boost::clear_vertex(marked_vertex, pBG);
      VertexNameMap name_map = boost::get(boost::vertex_name, pBG);
      boost::put(name_map, marked_vertex, kMarkedLabel);

      return std::tuple(marked_vertex, label, in_neighbours, out_neighbours);
  }

class BlissTimer {
  static std::chrono::nanoseconds time_in_ns;
  std::chrono::high_resolution_clock::time_point mStart;
public:
  std::chrono::milliseconds timeMS() { return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds(time_in_ns)); }
  
  void start() { mStart = std::chrono::high_resolution_clock::now(); }

  void end() {
    auto end = std::chrono::high_resolution_clock::now();
    time_in_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(end - mStart);
  }

};

std::chrono::nanoseconds BlissTimer::time_in_ns = std::chrono::nanoseconds(0);

  template<class PatternGraph>
  class Core {
    public:
      Core(PatternGraph *pattern, uint32_t markedVertex) : mParent(pattern), mMarkedVertexIdx(markedVertex) {
        mBoostGraph = pattern->boostGraph();
        /*
        for (const auto &label : pattern->labels()) {
          boost::add_vertex(VertexProperty(label), mBoostGraph);
        }
        for (const auto &[src, dst] : pattern->edges()) {
          boost::add_edge(src, dst, mBoostGraph);
        }
        */
        std::tie(mMarkedVertex, mLabel, mInNeighbors, mOutNeighbors) = markRemovedVertex(mBoostGraph, mMarkedVertexIdx);
        BlissTimer bt;
        bt.start();
        std::tie(mBlissGraph, mCanonicalPermutation, mCanonicalPermutationInverse) = GraphHelpers::canonicalForm(mBoostGraph);
        bt.end();
        //canonicalForm(*this);
      }

      std::shared_ptr<bliss::Digraph> blissGraph() const {
        return mBlissGraph;
      }

      PatternGraph *parent() const {
        return mParent;
      }

      uint32_t markedVertexIdx() const {
        return mMarkedVertexIdx;
      }

      uint32_t numVertices() const {
        // NOTE: assert will be false with edges being converted into intermediate vertices
        //assert(boost::num_vertices(mBoostGraph) == mBlissGraph->get_nof_vertices());
        return boost::num_vertices(mBoostGraph);
      }

      uint32_t label() const {
        return mLabel;
      }

      bool operator<(const Core &other) const {
        if (boost::num_vertices(mBoostGraph) != boost::num_vertices(other.mBoostGraph)) {
          return boost::num_vertices(mBoostGraph) < boost::num_vertices(other.mBoostGraph);
        }
        auto aHash = mBlissGraph->get_hash();
        auto bHash = other.mBlissGraph->get_hash();
        if (aHash != bHash) {
          return aHash < bHash;
        }
        auto blisscmp = mBlissGraph->cmp(*other.mBlissGraph);
        return blisscmp < 0;
      }

      std::vector<Permutation> findAllAutomorphisms(const std::vector<std::shared_ptr<Core>> &pToIgnore) {
        if (!mGeneratorsGenerated) {
          mGeneratorsGenerated = true;
          bliss::Stats stats;
          mBlissGraph->find_automorphisms(stats, [&](uint n, const uint *aut) {
                this->mGenerators.emplace_back(aut, aut + n);
              });
        }

        uint n = mBlissGraph->get_nof_vertices();

        return generateAutomorphismsFromGenerators(n, mGenerators);
      }

      std::vector<Permutation> findAutomorphisms(const std::vector<std::shared_ptr<Core>> &pToIgnore) {
        if (!mGeneratorsGenerated) {
          mGeneratorsGenerated = true;
          bliss::Stats stats;
          mBlissGraph->find_automorphisms(stats, [&](uint n, const uint *aut) {
                this->mGenerators.emplace_back(aut, aut + n);
              });
        }

        uint n = mBlissGraph->get_nof_vertices();
        bool allSame = false;

        for (auto &cg : pToIgnore) {
          bool isConstant = true;
          for (auto &gen : mGenerators) {
            bool isDifferent = false;
            std::vector<NeighborInfo*> repositionedOut(n);
            for (auto &outNeigh : cg->mOutNeighbors) {
              repositionedOut[outNeigh.first] = &outNeigh;
            }
            for (auto &outNeigh : cg->mOutNeighbors) {
              auto remapped = cg->mCanonicalPermutationInverse[gen[cg->mCanonicalPermutation[outNeigh.first]]];
              if (!repositionedOut[remapped]) {
                isDifferent = true;
                break;
              }
              if (repositionedOut[remapped]->second != outNeigh.second) {
                isDifferent = true;
                break;
              }
            }
            std::vector<NeighborInfo*> repositionedIn(n);
            for (auto &inNeigh : cg->mInNeighbors) {
              repositionedIn[inNeigh.first] = &inNeigh;
            }
            for (auto &inNeigh : cg->mInNeighbors) {
              auto remapped = cg->mCanonicalPermutationInverse[gen[cg->mCanonicalPermutation[inNeigh.first]]];
              if (!repositionedIn[remapped]) {
                isDifferent = true;
                break;
              }
              if (repositionedIn[remapped]->second != inNeigh.second) {
                isDifferent = true;
                break;
              }
            }
            if (isDifferent) {
              isConstant = false;
              break;
            }
          }
          if (isConstant) {
            allSame = true;
            break;
          }
        }
        if (!allSame) {
          return generateAutomorphismsFromGenerators(n, mGenerators);
        } else {
          Permutation result(n);
          std::iota(result.begin(), result.end(), 0);
          return {result};
        }
        return mGenerators;
      }

      std::vector<Permutation> generateAutomorphismsFromGenerators(unsigned int n,
          std::vector<Permutation> generators) {
        Permutation identity(n);
        std::iota(identity.begin(), identity.end(), 0);
        std::set<Permutation> result_set;
        std::queue<Permutation> to_process;
        to_process.push(identity);
        result_set.insert(identity);
        while(!to_process.empty()) {
          Permutation x = std::move(to_process.front());
          to_process.pop();
          for (Permutation gen : generators) {
            Permutation genx(n);
            for (unsigned int i = 0; i < n; ++i) {
              genx[i] = gen[x[i]];
            }
            if (result_set.insert(genx).second) {
              to_process.push(genx);
            }
          }
        }
        std::vector<Permutation> result(result_set.begin(), result_set.end());
        return result;
      }

      std::unique_ptr<PatternGraph> merge(Core& other, const Permutation &permutation) {
        VertexNameMap nameMap = boost::get(boost::vertex_name, mBoostGraph);
        boost::put(nameMap, mMarkedVertex, mLabel);
        for (const auto &inNeigh : mInNeighbors) {
          //boost::add_edge(inNeigh.first, mMarkedVertex, mBoostGraph);
          boost::add_edge(inNeigh.first, mMarkedVertex, inNeigh.second, mBoostGraph);
        }
        for (const auto &outNeigh : mOutNeighbors) {
          //boost::add_edge(mMarkedVertex, outNeigh.first, mBoostGraph);
          boost::add_edge(mMarkedVertex, outNeigh.first, outNeigh.second, mBoostGraph);
        }
        VertexDescriptor newVertex = boost::add_vertex(other.mLabel, mBoostGraph);
        for (const auto &inNeighOther : other.mInNeighbors) {
          NeighborInfo inNeigh{mCanonicalPermutationInverse[permutation[other.mCanonicalPermutation[inNeighOther.first]]], inNeighOther.second};
          //boost::add_edge(inNeigh.first, newVertex, mBoostGraph);
          boost::add_edge(inNeigh.first, newVertex, inNeigh.second, mBoostGraph);
        }
        for (const auto &outNeighOther : other.mOutNeighbors) {
          NeighborInfo outNeigh{mCanonicalPermutationInverse[permutation[other.mCanonicalPermutation[outNeighOther.first]]], outNeighOther.second};
          //boost::add_edge(newVertex, outNeigh.first, mBoostGraph);
          boost::add_edge(newVertex, outNeigh.first, outNeigh.second, mBoostGraph);
        }

        auto result = std::make_unique<PatternGraph>(mBoostGraph);
        boost::clear_vertex(newVertex, mBoostGraph);
        boost::remove_vertex(newVertex, mBoostGraph);
        std::tie(mMarkedVertex, mLabel, mInNeighbors, mOutNeighbors) = markRemovedVertex(mBoostGraph, mMarkedVertexIdx);
        return result;
      }

      std::unique_ptr<PatternGraph> mergeClique(Core& other, const Permutation &permutation, uint matchedVertex, Core& matchingBase, Core& third, const Permutation &permutation2) {
        assert(parent()->isClique());
        assert(other.parent()->isClique());
        assert(third.parent()->isClique());

        if (other.label() != third.label()) {
          return nullptr;
        }

        addMarkedVertex(mBoostGraph);

        std::map<uint, int32_t> newIns;
        std::map<uint, int32_t> newOuts;
        std::map<uint, int32_t> newIns2;
        std::map<uint, int32_t> newOuts2;

        for (const auto& inNeighOther : other.mInNeighbors) {
          uint vertex = mCanonicalPermutationInverse[permutation[other.mCanonicalPermutation[inNeighOther.first]]];
          newIns[vertex] = inNeighOther.second;
          if (vertex == matchedVertex) {
            newIns2[vertex] = inNeighOther.second;
          }
        }
        for (const auto& inNeighThird : third.mInNeighbors) {
          uint vertex = matchingBase.canonicalPermutationInverse()[permutation2[third.mCanonicalPermutation[inNeighThird.first]]];
          newIns2[vertex] = inNeighThird.second;
          if (vertex == mMarkedVertexIdx) {
            newIns[vertex] = inNeighThird.second;
          }
        }

        for (const auto& outNeighOther : other.mOutNeighbors) {
          uint vertex = mCanonicalPermutationInverse[permutation[other.mCanonicalPermutation[outNeighOther.first]]];
          newOuts[vertex] = outNeighOther.second;
          if (vertex == matchedVertex) {
            newOuts2[vertex] = outNeighOther.second;
          }
        }
        for (const auto& outNeighThird : third.mOutNeighbors) {
          uint vertex = matchingBase.canonicalPermutationInverse()[permutation2[third.mCanonicalPermutation[outNeighThird.first]]];
          newOuts2[vertex] = outNeighThird.second;
          if (vertex == mMarkedVertexIdx) {
            newOuts[vertex] = outNeighThird.second;
          }
        }

        VertexDescriptor newVertex = boost::add_vertex(other.mLabel, mBoostGraph);
        if (newIns != newIns2 || newOuts != newOuts2) {
          removeMarkedVertex(mBoostGraph, newVertex);
          return nullptr;
        }
        for (const auto& inNeigh: newIns) {
          //boost::add_edge(inNeigh.first, newVertex, mBoostGraph);
          boost::add_edge(inNeigh.first, newVertex, inNeigh.second, mBoostGraph);
        }
        for (const auto& outNeigh: newOuts) {
          //boost::add_edge(newVertex, outNeigh.first, mBoostGraph);
          boost::add_edge(newVertex, outNeigh.first, outNeigh.second, mBoostGraph);
        }
        auto result = std::make_unique<PatternGraph>(mBoostGraph);
        removeMarkedVertex(mBoostGraph, newVertex);
        return result;

      }

      void addMarkedVertex(BoostGraph &pBG) {
        VertexNameMap nameMap = boost::get(boost::vertex_name, pBG);
        boost::put(nameMap, mMarkedVertex, mLabel);
        for (const auto &inNeigh : mInNeighbors) {
          //boost::add_edge(inNeigh.first, mMarkedVertex, pBG);
          boost::add_edge(inNeigh.first, mMarkedVertex, inNeigh.second, pBG);
        }
        for (const auto &outNeigh : mOutNeighbors) {
          //boost::add_edge(mMarkedVertex, outNeigh.first, pBG);
          boost::add_edge(mMarkedVertex, outNeigh.first, outNeigh.second, pBG);
        }
      }

      void removeMarkedVertex(BoostGraph &pBG, VertexDescriptor markedVertex) {
        boost::clear_vertex(markedVertex, pBG);
        boost::remove_vertex(markedVertex, pBG);
        std::tie(mMarkedVertex, mLabel, mInNeighbors, mOutNeighbors) = markRemovedVertex(pBG, mMarkedVertexIdx);
      }

      Permutation canonicalPermutation() const {
        return mCanonicalPermutation;
      }

      Permutation canonicalPermutationInverse() const {
        return mCanonicalPermutationInverse;
      }

    private:
      PatternGraph* mParent;
      BoostGraph mBoostGraph;

      uint32_t mMarkedVertexIdx;
      VertexDescriptor mMarkedVertex;
      uint32_t mLabel;
      std::vector<std::pair<VertexDescriptor, int32_t>> mInNeighbors;
      std::vector<std::pair<VertexDescriptor, int32_t>> mOutNeighbors;

      std::shared_ptr<bliss::Digraph> mBlissGraph; // canonical_form_
      Permutation mCanonicalPermutation;
      Permutation mCanonicalPermutationInverse;

      bool mGeneratorsGenerated{false};
      std::vector<Permutation> mGenerators;

    public:
      static uint64_t bliss_timing_ns;
  };

  template <template<class> class CoreGraph, class PatternGraph>
  std::vector<CoreGraph<PatternGraph>> generateCoreGraphs(PatternGraph &pattern) {
    std::vector<CoreGraph<PatternGraph>> coreGraphs;
    for (uint32_t i = 0; i < boost::num_vertices(pattern.boostGraph()); ++i) {
      coreGraphs.emplace_back(&pattern, i);
    }
    return coreGraphs;
  }

  template <template<class> class CoreGraph, class PatternGraph>
  class Merger {
    public:
      using AllCores = std::vector<std::shared_ptr<CoreGraph<PatternGraph>>>;
      using CoreMap = std::map<PatternGraph*, std::vector<std::shared_ptr<CoreGraph<PatternGraph>>>>; // patternCores
      using CoreGroups = std::map<std::shared_ptr<CoreGraph<PatternGraph>>, std::vector<std::shared_ptr<CoreGraph<PatternGraph>>>,
            std::function<bool(const std::shared_ptr<CoreGraph<PatternGraph>>&, const std::shared_ptr<CoreGraph<PatternGraph>>&)>>;
      using InputCliques = std::set<CanonicalShorthand>;

    public:

    AllCores generateCores(std::vector<PatternGraph> &patterns) const {
      AllCores allCores;
      for ( auto &pattern : patterns ) {
        auto cores = generateCoreGraphs<CoreGraph>(pattern);
          for (auto &core : cores) {
            auto shared = std::make_shared<CoreGraph<PatternGraph>>(std::move(core));
            allCores.push_back(shared);
          }
      }
      return allCores;
    }

    InputCliques generateInputCliques(std::vector<PatternGraph> &patterns) const {
      InputCliques inputCliques;
      for ( auto &pattern : patterns) {
        if (pattern.isClique()) {
          inputCliques.insert(pattern.canonicalShorthand());
        }
      }
      return inputCliques;
    }

    CoreMap generateCoreMap(AllCores allCores) const {
      CoreMap coreMap;
      for (auto core_shared : allCores) {
        coreMap[core_shared->parent()].push_back(core_shared);
      }
      return coreMap;
    }

    CoreGroups generateCoreGroups(AllCores allCores) const {
        CoreGroups coreGroups([](const std::shared_ptr<CoreGraph<PatternGraph>> &a, const std::shared_ptr<CoreGraph<PatternGraph>> &b) {
          return *a < *b;
        });
        for (auto& core : allCores) {
          coreGroups[core].push_back(core);
        }
        return coreGroups;
    }

    std::tuple<AllCores, CoreMap, CoreGroups, InputCliques> initialize(std::vector<PatternGraph> &patterns) {

        AllCores allCores;
        CoreMap coreMap; // patternCores
        InputCliques inputCliques;

        for (auto &pattern : patterns) {
          if (pattern.isClique()) {
            inputCliques.insert(pattern.canonicalShorthand());
          }
          auto cores = generateCoreGraphs<CoreGraph>(pattern);
          for (auto &core : cores) {
            auto shared = std::make_shared<CoreGraph<PatternGraph>>(std::move(core));
            coreMap[&pattern].push_back(shared);
            allCores.push_back(shared);
          }
        }

        CoreGroups coreGroups([](const std::shared_ptr<CoreGraph<PatternGraph>> &a, const std::shared_ptr<CoreGraph<PatternGraph>> &b) {
          return *a < *b;
        });
        for (auto& core : allCores) {
          coreGroups[core].push_back(core);
        }
        return std::make_tuple(allCores, coreMap, coreGroups, inputCliques);
      }

    std::tuple<bool, typename CoreGroups::iterator> checkCliques(const CoreGraph<PatternGraph> &core_i,
        CoreGroups &coreGroups,
        CoreMap &coreMap) const {
      bool core_i_check_cliques = core_i.parent()->isClique()
        && core_i.markedVertexIdx() == core_i.parent()->cliqueCanonPermutation().back();
      auto cliqueGroup = coreGroups.end();
      if (core_i_check_cliques) {
        cliqueGroup = coreGroups.find(coreMap[core_i.parent()]
            [core_i.parent()->cliqueCanonPermutation()[core_i.numVertices() - 2]]);
        core_i_check_cliques = cliqueGroup != coreGroups.end();
      }
      return std::make_tuple(core_i_check_cliques, cliqueGroup);
    }

  };

  template <class PatternGraph>
  std::pair<std::set<PatternGraph>, uint> extendDuplicates(std::vector<PatternGraph> &from) {
    uint duplicatesIncluded = 0;
    Merger<Core, PatternGraph> merger;
    auto [allCores, coreMap, coreGroups, inputCliques] = merger.initialize(from);

    std::set<std::unique_ptr<PatternGraph>, std::function<bool(const std::unique_ptr<PatternGraph>&, const std::unique_ptr<PatternGraph>&)>> results(
        [](const std::unique_ptr<PatternGraph> &a, const std::unique_ptr<PatternGraph> &b) {
          return *a < *b;
        });
    for (auto &[coreKey, coreGroup] : coreGroups) {
      for (auto iiter = coreGroup.begin(); iiter < coreGroup.end(); ++iiter) {
        Core<PatternGraph> &core_i = **iiter;
        auto [core_i_check_cliques, cliqueGroup] = merger.checkCliques(core_i, coreGroups, coreMap);
        for (auto jiter = iiter; jiter < coreGroup.end(); ++jiter) {
          Core<PatternGraph> &core_j = **jiter;
          std::vector<std::shared_ptr<Core<PatternGraph>>> fixers;
          fixers.push_back(*iiter);
          fixers.push_back(*jiter);
          std::vector<Permutation> automorphisms = coreKey->findAllAutomorphisms(fixers);
          for (auto &automorphism : automorphisms) {
            std::unique_ptr<PatternGraph> newPattern = core_i.merge(*jiter->get(), automorphism);
            auto res = results.insert(std::move(newPattern));
            duplicatesIncluded++;
            if (core_i_check_cliques && core_j.parent()->isClique()) {
              uint matchedVertex = core_i.parent()->cliqueCanonPermutation()[core_i.numVertices() - 2];
              Core<PatternGraph> &cliqueBase = *(coreMap[core_i.parent()][matchedVertex]);
              for (auto& clique3 : cliqueGroup->second) {
                if (!clique3->parent()->isClique()) {
                  continue;
                }
                if (clique3->label() != core_j.label()) {
                  continue;
                }
                std::vector<Permutation> automorphisms2 = clique3->findAutomorphisms(
                    std::vector<std::shared_ptr<Core<PatternGraph>>>(0));
                for (auto& automorphism2 : automorphisms2) {
                  std::unique_ptr<PatternGraph> p = core_i.mergeClique(core_j, automorphism, matchedVertex, cliqueBase, *clique3, automorphism2);
                  if (p) {
                    bool allChildrenGood = true;
                    for (uint exclude = 0; exclude < p->nodes(); ++exclude) {
                      std::vector<unsigned int> allBut(p->nodes() - 1);
                      std::generate(allBut.begin(), allBut.end(),
                          [exclude, n=-1]() mutable { return (++n == exclude) ? ++n : n; } );
                      const CanonicalShorthand cs = GraphHelpers::getCanonicalShorthand(p->boostGraph(), allBut);
                      PatternGraph p2(cs);
                      if (inputCliques.find(p2.canonicalShorthand()) == inputCliques.end()) {
                        allChildrenGood = false;
                        break;
                      }
                    }
                    if (allChildrenGood) {
                      duplicatesIncluded++;
                      results.insert(std::move(p));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

#if 0
    std::vector<PatternGraph> newPatterns;
    newPatterns.reserve(results.size());
    for (auto &result : results) {
      newPatterns.push_back(std::move(*result));
    }
#else
    std::set<PatternGraph> newPatterns;
    for (auto &result : results) {
      newPatterns.insert(std::move(*result));
    }
#endif
    return std::make_pair(newPatterns, duplicatesIncluded);
  }
  
  template <class PatternGraph>
  std::pair<std::set<PatternGraph>, uint> extend(std::vector<PatternGraph> &from) {
    uint prunedPatternsCount = 0;
    Merger<Core, PatternGraph> merger;
    auto [allCores, coreMap, coreGroups, inputCliques] = merger.initialize(from);

    std::set<std::unique_ptr<PatternGraph>, std::function<bool(const std::unique_ptr<PatternGraph>&, const std::unique_ptr<PatternGraph>&)>> results(
        [](const std::unique_ptr<PatternGraph> &a, const std::unique_ptr<PatternGraph> &b) {
          return *a < *b;
        });
    for (auto &[coreKey, coreGroup] : coreGroups) {
      for (auto iiter = coreGroup.begin(); iiter < coreGroup.end(); ++iiter) {
        Core<PatternGraph> &core_i = **iiter;
        auto [core_i_check_cliques, cliqueGroup] = merger.checkCliques(core_i, coreGroups, coreMap);
        for (auto jiter = iiter; jiter < coreGroup.end(); ++jiter) {
          Core<PatternGraph> &core_j = **jiter;
          std::vector<std::shared_ptr<Core<PatternGraph>>> fixers;
          fixers.push_back(*iiter);
          fixers.push_back(*jiter);
          std::vector<Permutation> automorphisms = coreKey->findAutomorphisms(fixers);
          for (auto &automorphism : automorphisms) {
            std::unique_ptr<PatternGraph> newPattern = core_i.merge(*jiter->get(), automorphism);
            results.insert(std::move(newPattern));
            prunedPatternsCount++;
            if (core_i_check_cliques && core_j.parent()->isClique()) {
              uint matchedVertex = core_i.parent()->cliqueCanonPermutation()[core_i.numVertices() - 2];
              Core<PatternGraph> &cliqueBase = *(coreMap[core_i.parent()][matchedVertex]);
              for (auto& clique3 : cliqueGroup->second) {
                if (!clique3->parent()->isClique()) {
                  continue;
                }
                if (clique3->label() != core_j.label()) {
                  continue;
                }
                std::vector<Permutation> automorphisms2 = clique3->findAutomorphisms(
                    std::vector<std::shared_ptr<Core<PatternGraph>>>(0));
                for (auto& automorphism2 : automorphisms2) {
                  std::unique_ptr<PatternGraph> p = core_i.mergeClique(core_j, automorphism, matchedVertex, cliqueBase, *clique3, automorphism2);
                  if (p) {
                    bool allChildrenGood = true;
                    for (uint exclude = 0; exclude < p->nodes(); ++exclude) {
                      std::vector<unsigned int> allBut(p->nodes() - 1);
                      std::generate(allBut.begin(), allBut.end(),
                          [exclude, n=-1]() mutable { return (++n == exclude) ? ++n : n; } );
                      const CanonicalShorthand cs = GraphHelpers::getCanonicalShorthand(p->boostGraph(), allBut);
                      PatternGraph p2(cs);
                      if (inputCliques.find(p2.canonicalShorthand()) == inputCliques.end()) {
                        allChildrenGood = false;
                        break;
                      }
                    }
                    if (allChildrenGood) {
                      prunedPatternsCount++;
                      results.insert(std::move(p));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

#if 0
    std::vector<PatternGraph> newPatterns;
    newPatterns.reserve(results.size());
    for (auto &result : results) {
      newPatterns.push_back(std::move(*result));
    }
#else
    std::set<PatternGraph> newPatterns;
    for (auto &result : results) {
      newPatterns.insert(std::move(*result));
    }
#endif
    return std::make_pair(newPatterns, prunedPatternsCount);
  }
 
  template <class PatternGraph>
  std::pair<std::set<PatternGraph>, uint> extendTiming(std::vector<PatternGraph> &from) {
    uint prunedPatternsCount = 0;
    Merger<Core, PatternGraph> merger;
    auto [allCores, coreMap, coreGroups, inputCliques] = merger.initialize(from);

    std::set<std::unique_ptr<PatternGraph>, std::function<bool(const std::unique_ptr<PatternGraph>&, const std::unique_ptr<PatternGraph>&)>> results(
        [](const std::unique_ptr<PatternGraph> &a, const std::unique_ptr<PatternGraph> &b) {
          return *a < *b;
        });
    for (auto &[coreKey, coreGroup] : coreGroups) {
      for (auto iiter = coreGroup.begin(); iiter < coreGroup.end(); ++iiter) {
        Core<PatternGraph> &core_i = **iiter;
        auto [core_i_check_cliques, cliqueGroup] = merger.checkCliques(core_i, coreGroups, coreMap);
        for (auto jiter = iiter; jiter < coreGroup.end(); ++jiter) {
          Core<PatternGraph> &core_j = **jiter;
          std::vector<std::shared_ptr<Core<PatternGraph>>> fixers;
          fixers.push_back(*iiter);
          fixers.push_back(*jiter);
        BlissTimer bt;
        bt.start();
          std::vector<Permutation> automorphisms = coreKey->findAutomorphisms(fixers);
        bt.end();
          
          for (auto &automorphism : automorphisms) {
            std::unique_ptr<PatternGraph> newPattern = core_i.merge(*jiter->get(), automorphism);
            results.insert(std::move(newPattern));
            prunedPatternsCount++;
            if (core_i_check_cliques && core_j.parent()->isClique()) {
              uint matchedVertex = core_i.parent()->cliqueCanonPermutation()[core_i.numVertices() - 2];
              Core<PatternGraph> &cliqueBase = *(coreMap[core_i.parent()][matchedVertex]);
              for (auto& clique3 : cliqueGroup->second) {
                if (!clique3->parent()->isClique()) {
                  continue;
                }
                if (clique3->label() != core_j.label()) {
                  continue;
                }
              BlissTimer bt2;
              bt2.start();
                std::vector<Permutation> automorphisms2 = clique3->findAutomorphisms(
                    std::vector<std::shared_ptr<Core<PatternGraph>>>(0));
              bt2.end();
                for (auto& automorphism2 : automorphisms2) {
                  std::unique_ptr<PatternGraph> p = core_i.mergeClique(core_j, automorphism, matchedVertex, cliqueBase, *clique3, automorphism2);
                  if (p) {
                    bool allChildrenGood = true;
                    for (uint exclude = 0; exclude < p->nodes(); ++exclude) {
                      std::vector<unsigned int> allBut(p->nodes() - 1);
                      std::generate(allBut.begin(), allBut.end(),
                          [exclude, n=-1]() mutable { return (++n == exclude) ? ++n : n; } );
                      const CanonicalShorthand cs = GraphHelpers::getCanonicalShorthand(p->boostGraph(), allBut);
                      PatternGraph p2(cs);
                      if (inputCliques.find(p2.canonicalShorthand()) == inputCliques.end()) {
                        allChildrenGood = false;
                        break;
                      }
                    }
                    if (allChildrenGood) {
                      prunedPatternsCount++;
                      results.insert(std::move(p));
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

#if 0
    std::vector<PatternGraph> newPatterns;
    newPatterns.reserve(results.size());
    for (auto &result : results) {
      newPatterns.push_back(std::move(*result));
    }
#else
    std::set<PatternGraph> newPatterns;
    for (auto &result : results) {
      newPatterns.insert(std::move(*result));
    }
#endif
    return std::make_pair(newPatterns, prunedPatternsCount);
  }
};
