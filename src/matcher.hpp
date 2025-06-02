#pragma once


#include <algorithm>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <tuple>
#include <map>
#include <vector>
#include <optional>
#include <functional>

#include <iostream>

#include <vf3lib/VFLib.h>
#include <vf3lib/VF3LightSubState.hpp>
#include <vf3lib/VF3NodeSorter.hpp>

#include <fmt/format.h>

#include "NodeClassifier.hpp"
#include "instruction.hpp"
#include "matchState.hpp"

namespace Vf3Frequency {

  template <typename DataGraph>
  struct DataGraphStruct {
    const DataGraph mDataGraph;
    vflib::ARGraph<uint32_t, uint32_t> mDataGraphARG;
    vflib::NodeClassifier<uint32_t, uint32_t> mDataGraphClassifier;
    vflib::VF3NodeSorter<data_t, data_t, vflib::SubIsoNodeProbability<data_t, data_t>> mDataGraphSorter;
    DataGraphStruct() = default;
    DataGraphStruct(DataGraphStruct<DataGraph> &&) = default;
    DataGraphStruct(const DataGraph &dataGraph) : mDataGraph(dataGraph),
      mDataGraphARG(vflib::ARGraph<uint32_t, uint32_t>(&mDataGraph)),
      mDataGraphClassifier(vflib::NodeClassifier<uint32_t, uint32_t>(&mDataGraphARG)),
      mDataGraphSorter(vflib::VF3NodeSorter<data_t, data_t, vflib::SubIsoNodeProbability<data_t, data_t>>(&mDataGraphARG)) {}
  };

  template<typename Node1, typename Node2,
    typename Edge1, typename Edge2,
    typename NodeComparisonFunctor = vflib::EqualityComparator<Node1, Node2>,
    typename EdgeComparisonFunctor = vflib::EqualityComparator<Edge1, Edge2>>
      class VF3LightMining : public vflib::VF3LightSubState<Node1, Node2, Edge1, Edge2, NodeComparisonFunctor, EdgeComparisonFunctor> {

        public:

          VF3LightMining(const vflib::ARGraph<Node1, Edge1> *g1, const vflib::ARGraph<Node2, Edge2> *g2,
              uint32_t* class_1, uint32_t* class_2, uint32_t nclass,
              std::shared_ptr<std::vector<bool>> used,
              uint32_t* order,
              uint32_t pThreshold,
              bool induced = false // edge induced
              )
            : vflib::VF3LightSubState<Node1, Node2, Edge1, Edge2, NodeComparisonFunctor, EdgeComparisonFunctor>(g1, g2, class_1, class_2, nclass, order, induced) {
              mThreshold = std::make_shared<uint32_t>(pThreshold);
              mCount = std::make_shared<uint32_t>(0);
              mUsed = used; // std::make_shared<std::vector<bool>>(g2->NodeCount(), false);

          }
          VF3LightMining(const VF3LightMining& other) = default;
          virtual ~VF3LightMining() = default;

          bool accept() {
            mCount.get()[0]++;
            for (int i = 0; i < this->n1; i++) {
              assert(this->core_1[i] != vflib::NULL_NODE);
              assert(this->core_1[i] < this->g2->NodeCount());
              mUsed->at(this->core_1[i]) = true;
            }
            //fmt::print("{}, core_1[0] {} = {}, core_1[1] {} = {}\n", this->core_len, this->core_1[0], mUsed->at(this->core_1[0]), this->core_1[1], mUsed->at(this->core_1[1]));
            return mThreshold.get()[0] != 0 && mCount.get()[0] >= mThreshold.get()[0];
          }

          bool MiningPossible(uint32_t n1, uint32_t n2) {
            //fmt::print("{}, MiningPossible({}, {}) : {}, {}, count:{}\n", this->core_len, n1, n2, mUsed->at(n1), mUsed->at(n2), mCount.get()[0]);
            assert(n2 < this->g2->NodeCount());
            return !mUsed->at(n2);
          }

        private:
          std::shared_ptr<uint32_t> mThreshold;
          std::shared_ptr<uint32_t> mCount;

          std::shared_ptr<std::vector<bool>> mUsed;

      };

  template<typename VFState>
  class MiningVisitor : public vflib::MatchingVisitor<VFState>  {
    public:
      virtual bool operator()(VFState &s) {
        return s.accept();
      }
      virtual ~MiningVisitor() = default;
  };

  template <typename VFState>
  class MiningMatchingEngine : public vflib::MatchingEngine<VFState> {
    public:
      MiningMatchingEngine(bool storeSolutions = false) : vflib::MatchingEngine<VFState>(storeSolutions) {}
      MiningMatchingEngine(vflib::MatchingVisitor<VFState> *visit, bool storeSolutions = false) : vflib::MatchingEngine<VFState>(visit, storeSolutions) {}
      virtual ~MiningMatchingEngine() = default;
      bool mExit = false;
      virtual bool FindAllMatchingsRec(VFState &s) {
        if (s.IsGoal()) {
          solCount++;
          if (storeSolutions) {
            vflib::MatchingSolution sol;
            s.GetCoreSet(sol);
            solutions.push_back(sol);
          }
          mExit = true;
          if (visit) {
            return (*visit)(s);
          }
          return false;
        }
        if (s.IsDead()) {
          return false;
        }
        uint32_t n1 = vflib::NULL_NODE, n2 = vflib::NULL_NODE;
        while(s.NextPair(&n1, &n2, n1, n2)) {
          if (!s.MiningPossible(n1, n2)) continue;
          if (s.IsFeasiblePair(n1, n2)) {
            VFState s1(s);
            s1.AddPair(n1, n2);
            if (FindAllMatchingsRec(s1)) {
              return true;
            }
            if (mExit) return false;
          }
        }
        return false;
      }
      virtual bool FindAllMatchings(VFState &s) override {
        if (s.IsGoal()) {
          solCount++;
          if (storeSolutions) {
            vflib::MatchingSolution sol;
            s.GetCoreSet(sol);
            solutions.push_back(sol);
          }
          if (visit) {
            return (*visit)(s);
          }
          return false;
        }
        if (s.IsDead()) {
          return false;
        }
        uint32_t n1 = vflib::NULL_NODE, n2 = vflib::NULL_NODE;
        while(s.NextPair(&n1, &n2, n1, n2)) {
          if (!s.MiningPossible(n1, n2)) continue;
          if (s.IsFeasiblePair(n1, n2)) {
            VFState s1(s);
            s1.AddPair(n1, n2);
            if (FindAllMatchingsRec(s1)) {
              return true;
            }
            mExit = false;
          }
        }
        return false;
      }
      inline void GetSolutions(std::vector<MatchState::VertexEmbedding> &pSols) {
        pSols.reserve(solutions.size());
        for (auto &sol : solutions) {
          MatchState::VertexEmbedding embedding;
          for (auto &pair : sol) {
            embedding.push_back(pair.second);
          }
          pSols.push_back(embedding);
        }
      }
    protected:
      using vflib::MatchingEngine<VFState>::solCount;
      using vflib::MatchingEngine<VFState>::storeSolutions;
      using vflib::MatchingEngine<VFState>::solutions;
      using vflib::MatchingEngine<VFState>::visit;
  };

  template <class PatternGraph, class DataGraph>
  class DFSMatcher {

    public:
      DFSMatcher(DataGraphStruct<DataGraph> &dataGraphStruct) : mDataGraphStruct(dataGraphStruct)
    {
      mUsed = std::make_shared<std::vector<bool>>(dataGraphStruct.mDataGraph.nodes(), false);
      }

      template <bool stepEmbeddings>
      MatchState::Vf3MatchState matchState(const PatternGraph &pPattern, int pThreshold, bool pRandom = false) {
        MatchState::Vf3MatchState state(mDataGraphStruct.mDataGraph.nodes(), pPattern.nodes(), pThreshold, stepEmbeddings);
        vf3matcher(pPattern, state, pRandom, stepEmbeddings);
        return state;
      }

      MatchState::Vf3MatchState vf3matcher(const PatternGraph &pattern, MatchState::Vf3MatchState &matchState, bool pRandom, bool pRecordEmbedding) {
        const vflib::ARGLoader<uint32_t, uint32_t>* patt_loader = &pattern;
        vflib::ARGraph<uint32_t, uint32_t> patt_graph(patt_loader);
        auto n1 = patt_graph.NodeCount();
        auto n2 = mDataGraphStruct.mDataGraphARG.NodeCount();
        using state_t = VF3LightMining<uint32_t, uint32_t, uint32_t, uint32_t>;
        MiningVisitor<state_t> visitor;
        MiningMatchingEngine<state_t> me(&visitor, pRecordEmbedding);
        std::fill(mUsed->begin(), mUsed->end(), false);

        std::vector<uint32_t> class_patt;
        std::vector<uint32_t> class_data;
        uint32_t classes_count = 0;

        vflib::FastCheck<uint32_t, uint32_t, uint32_t, uint32_t> check(&patt_graph, &mDataGraphStruct.mDataGraphARG);
        if (check.CheckSubgraphIsomorphism()) {
          vflib::NodeClassifier<uint32_t, uint32_t> classifier2(&patt_graph, mDataGraphStruct.mDataGraphClassifier);
          class_patt = classifier2.GetClasses();
          class_data = mDataGraphStruct.mDataGraphClassifier.GetClasses();
          classes_count =  mDataGraphStruct.mDataGraphClassifier.CountClasses();
        }

        bool statified = false;
        if (check.CheckSubgraphIsomorphism()) {
          std::vector<uint32_t> sorted = mDataGraphStruct.mDataGraphSorter.SortNodes(&patt_graph);

          state_t s0(&patt_graph, &mDataGraphStruct.mDataGraphARG, class_patt.data(), class_data.data(), classes_count, mUsed, sorted.data(), matchState.mThreshold);
          statified = me.FindAllMatchings(s0);
        }

        matchState.mCount = me.GetSolutionsCount();
        if (matchState.mRecordEmbedding) {
          me.GetSolutions(matchState.mEmbeddings);
        }

        return matchState;
      }

    private:
      DataGraphStruct<DataGraph> &mDataGraphStruct;
      std::shared_ptr<std::vector<bool>> mUsed;
  };

  template <bool stepEmbeddings, bool trace, class DataGraph, class PatternGraph>
  std::vector<PatternGraph> matchThres(DataGraphStruct<DataGraph> &dataGraphStruct, const std::set<PatternGraph> &patterns, int maxThreshold, float accuracy) {
    DFSMatcher<PatternGraph, DataGraph> matcher(dataGraphStruct);
    std::vector<PatternGraph> frequentPatterns;
    auto firstPatternNodes = patterns.begin()->nodes();
    int minThreshold = maxThreshold / firstPatternNodes;
    int threshold = std::lerp(minThreshold, maxThreshold, accuracy);
    std::copy_if(patterns.begin(), patterns.end(), std::back_inserter(frequentPatterns), [&matcher, threshold](const PatternGraph &pattern) {
      auto matchState = matcher.template matchState<stepEmbeddings>(pattern, threshold);
      if (trace) {
        MatchState::printEmbedding(matchState, pattern);
      }
      return matchState.mCount >= threshold;
    });
    return frequentPatterns;
  }

}; // namespace Vf3Frequency

