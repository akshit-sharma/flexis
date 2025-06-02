#pragma once

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <set>
#include <unordered_map>
#include <vector>

#include "common.h"

struct NeighborCSR {
  std::vector<uint> compressedRow;
  std::vector<uint> vec;
  std::vector<uint> edgeLabel;
};
struct OutEdgeIterator{
  std::vector<uint>::const_iterator vec;
  std::vector<uint>::const_iterator vecEnd;
  std::vector<uint>::const_iterator edgeLabel;
};

template<typename T>
struct CSRIterator {
  using value_type = T;
  using difference_type = ptrdiff_t;
  using iterator = typename std::vector<T>::iterator;
  using size_type = size_t;
  using const_reference = const typename std::vector<T>::const_reference;

  CSRIterator(const_reference start, const_reference finish) : mStart(start), mFinish(finish), mSize(finish - start) {}

  iterator begin() const {
    return mStart;
  }

  iterator end() const {
    return mFinish;
  }

  value_type operator[](size_type i) const {
    assert(i < mSize);
    return mStart[i];
  }

  private:
    iterator mStart;
    iterator mFinish;
    ptrdiff_t mSize;
};

template<typename T>
struct VertexCSR {
  using value_type = T;
  std::vector<T> compressedRow;
  std::vector<T> vec;
};

struct NeighborLabelCSR {
  std::vector<uint> compressedLabels;
  std::vector<uint> vec;
  std::vector<uint> label;
  std::vector<uint> edgeLabel;
};

struct NeighborLabelIterator {
  std::vector<uint>::const_iterator vec;
  std::vector<uint>::const_iterator vecEnd;
  std::vector<uint>::const_iterator label;
  std::vector<uint>::const_iterator edgeLabel;
};

namespace DataHelpers {
    using Labels = Type::Labels;
    using Edges = Type::LabelledEdges;
    using LabelsVector = Type::LabelsVector;
    using Remap = Type::RemapEdges;
    using RemapReverse = Type::RemapEdgesReverse;

    uint numLabels(LabelsVector &vertexLabels) { // NOLINT(misc-definitions-in-headers)
      std::set<uint> uniqueLabels(vertexLabels.begin(), vertexLabels.end());
      return uniqueLabels.size();
    }

    template <template<typename, typename> typename Map, typename Key, typename Value>
    Map<Key, Value> pack(const std::vector<Key> &keys, const std::vector<Value> &values) {
      assert(keys.size() == values.size());
      Map<Key, Value> map;
      for (size_t i = 0; i < keys.size(); ++i) {
        map[keys[i]] = values[i];
      }
      return map;
    }

    template <class Map>
    std::tuple<std::vector<typename Map::key_type>, std::vector<typename Map::mapped_type>> unpack(const Map &edgeMap) { // NOLINT(misc-definitions-in-headers)
      std::vector<typename Map::key_type> neighbors;
      std::vector<typename Map::mapped_type> edgeLabels;
      for (auto &edge : edgeMap) {
        neighbors.emplace_back(edge.first);
        edgeLabels.emplace_back(edge.second);
      }
      return std::make_tuple(neighbors, edgeLabels);
    }

    std::tuple<std::vector<uint>, std::vector<uint>> unpack(std::vector<std::tuple<uint, uint>> &edgeMap) { // NOLINT(misc-definitions-in-headers)
      std::vector<uint> neighbors;
      std::vector<uint> edgeLabels;
      for (auto &edge : edgeMap) {
        neighbors.emplace_back(std::get<0>(edge));
        edgeLabels.emplace_back(std::get<1>(edge));
      }
      return std::make_tuple(neighbors, edgeLabels);
    }
    NeighborCSR processEdge(size_t n, Edges &edges, const LabelsVector &vertexLabels) { // NOLINT(misc-definitions-in-headers)
      assert(n == vertexLabels.size());
      NeighborCSR edgesCSR;
      edgesCSR.compressedRow.reserve(n+1);
      edgesCSR.compressedRow.emplace_back(edgesCSR.vec.size());
      for (int i = 0; i < n; i++) {
        auto [neighbors, edgeLabels] = unpack(edges[i]);
        assert(std::is_sorted(neighbors.begin(), neighbors.end()));
        std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(edgesCSR.vec));
        std::copy(edgeLabels.begin(), edgeLabels.end(), std::back_inserter(edgesCSR.edgeLabel));
        edgesCSR.compressedRow.emplace_back(edgesCSR.vec.size());
      }
      assert(edgesCSR.compressedRow.size() == n+1);
      return edgesCSR;
    }
    VertexCSR<uint> csrVertexLabels(const LabelsVector &vertexLabels) { // NOLINT(misc-definitions-in-headers)
      std::unordered_map<uint, std::set<uint>> labelToVertices;
      for (int i = 0; i < vertexLabels.size(); i++) {
        labelToVertices[vertexLabels[i]].insert(i);
      }
      VertexCSR<uint> vertexCSR;
      vertexCSR.compressedRow.reserve(labelToVertices.size()+1);
      vertexCSR.compressedRow.emplace_back(vertexCSR.vec.size());
      for (uint i = 0; i < labelToVertices.size(); i++) {
        std::copy(labelToVertices[i].begin(), labelToVertices[i].end(), std::back_inserter(vertexCSR.vec));
        vertexCSR.compressedRow.emplace_back(vertexCSR.vec.size());
      }
      assert(vertexCSR.compressedRow.size() == labelToVertices.size()+1);
      return vertexCSR;
    }
    NeighborLabelCSR csrNeighborLabels(uint n, Edges &edges, const LabelsVector &vertexLabels) { // NOLINT(misc-definitions-in-headers)
      assert(n == vertexLabels.size());
      std::unordered_map<uint, std::set<uint>> labelToVertices;
      uint maxLabels = 0;
      for (int i = 0; i < vertexLabels.size(); i++) {
        labelToVertices[vertexLabels[i]].insert(i);
        maxLabels = std::max(maxLabels, vertexLabels[i]+1);
      }
      NeighborLabelCSR neighborLabelCSR;
  //    neighborLabelCSR.compressedRow.reserve(n+1);
      neighborLabelCSR.compressedLabels.reserve(maxLabels * n + 1);
      auto map2vecTup = [](std::map<uint, uint> &map) {
        std::vector<std::tuple<uint, uint>> vec;
        for (auto &entry : map) {
          vec.emplace_back(entry.first, entry.second);
        }
        return vec;
      };

  //    neighborLabelCSR.compressedRow.emplace_back(neighborLabelCSR.compressedLabels.size());
      for (int i = 0; i < n; i++) {
        std::vector<std::tuple<uint, uint>> neighborsTup = map2vecTup(edges[i]);
        std::stable_sort(neighborsTup.begin(), neighborsTup.end(),
            [&vertexLabels](auto a, auto b) { return vertexLabels[std::get<0>(a)] < vertexLabels[std::get<0>(b)]; });
        auto currentSize = neighborLabelCSR.vec.size();
        auto [neighbors, edgeLabels] = unpack(neighborsTup);
        std::copy(neighbors.begin(), neighbors.end(), std::back_inserter(neighborLabelCSR.vec));
        std::copy(edgeLabels.begin(), edgeLabels.end(), std::back_inserter(neighborLabelCSR.edgeLabel));
        std::vector<uint> neighborLabelVector;
        for (auto label : neighbors) {
          neighborLabelVector.emplace_back(vertexLabels[label]);
        }
        for (auto label = 0; label < maxLabels; label++) {
          auto it = std::lower_bound(neighborLabelVector.begin(), neighborLabelVector.end(), label);
          auto distance = it - neighborLabelVector.begin();
          neighborLabelCSR.compressedLabels.emplace_back(currentSize + distance);
        }
 //       neighborLabelCSR.compressedRow.emplace_back(neighborLabelCSR.vec.size());
      }
      neighborLabelCSR.compressedLabels.emplace_back(neighborLabelCSR.vec.size());
//      assert(neighborLabelCSR.compressedRow.size() == n+1);
      assert(neighborLabelCSR.compressedLabels.size() == maxLabels * n + 1);
      return neighborLabelCSR;
    }

    RemapReverse reverseMap(const Remap &remap) { // NOLINT(misc-definitions-in-headers)
      RemapReverse reverse;
      for (auto &entry : remap) {
        reverse[entry.second] = entry.first;
      }
      return reverse;
    }

    struct GraphConstructor{
      uint n;
      uint nz;
      Edges outgoingEdges;
      Edges incomingEdges;
      LabelsVector vertexLabels;
      Remap remapEdgeLabels;
    };
};
