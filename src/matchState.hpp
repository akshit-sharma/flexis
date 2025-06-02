#pragma once

#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>

namespace MatchState {

  using LocalNodeID = std::uint32_t;
  using VertexEmbedding = std::vector<LocalNodeID>;

  struct Vf3MatchState {
    VertexEmbedding mEmbedding;

    uint32_t mThreshold{0};

    bool mRecordEmbedding{false};
    std::vector<VertexEmbedding> mEmbeddings;

    bool mEarlyStop{false};
    uint32_t mCount{0};
    //std::vector<bool> mMarked;
    //std::vector<bool> mUsed;

    Vf3MatchState(int dataSize, int patternSize, int pThreashold, bool pRecordEmbedding = false)
      : mEmbedding(VertexEmbedding(patternSize)),
      mThreshold(pThreashold),
      //mMarked(std::vector<bool>(dataSize, false)),
      //mUsed(std::vector<bool>(dataSize, false)),
      mRecordEmbedding(pRecordEmbedding)
    {}

    void used(LocalNodeID node) {
      //mUsed[node] = true;
    }

    bool isUsed(LocalNodeID node) const {
      //return mUsed[node];
      return true;
    }

    void mark(LocalNodeID node) {
     // mMarked[node] = true;
    }

    void unmark(LocalNodeID node) {
      //mMarked[node] = false;
    }

    bool isMarked(LocalNodeID node) const {
      //return mUsed[node] || mMarked[node];
      return true;
    }

  };

  struct MatchState {
    VertexEmbedding mEmbedding;

    bool mEarlyStop{false};
    std::vector<bool> mMarked;
    std::vector<bool> mUsed;

    MatchState(int dataSize, int patternSize)
      : mEmbedding(VertexEmbedding(patternSize)),
      mMarked(std::vector<bool>(dataSize, false)),
      mUsed(std::vector<bool>(dataSize, false))
    {}

    void used(LocalNodeID node) {
      mUsed[node] = true;
    }

    bool isUsed(LocalNodeID node) const {
      return mUsed[node];
    }

    void mark(LocalNodeID node) {
      mMarked[node] = true;
    }

    void unmark(LocalNodeID node) {
      mMarked[node] = false;
    }

    bool isMarked(LocalNodeID node) const {
      return mUsed[node] || mMarked[node];
    }

  };


  struct MatchStruct {
    uint32_t mCount{0};
    uint32_t mThreshold{0};

    std::vector<VertexEmbedding> mEmbeddings;
  };

  template<class PatternGraph>
  void printEmbedding(Vf3MatchState &matchState, const PatternGraph &pattern) {
    std::cout << "Count " << matchState.mCount << " : ( ";
    pattern.visualize(std::cout);
    std::cout << " ) | ";
    for (auto &embedding : matchState.mEmbeddings) {
      std::copy(embedding.begin(), embedding.end(), std::ostream_iterator<int>(std::cout, " "));
    }
    std::cout << std::endl;
  }

};
