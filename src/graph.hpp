#pragma once

#include <algorithm>
#include <cassert>
#include <set>
#include <fstream>
#include <sstream>
#include <iterator>
#include <map>
#include <unordered_map>
#include <vector>

#include <vf3lib/ARGraph.hpp>
#include <vf3lib/VFLib.h>

#include "common.h"
#include "dataHelpers.hpp"

class Vf3Graph : public vflib::ARGLoader<uint32_t, uint32_t> {
  public:
    uint nodes() const { return n_; };
    uint nonZeros() const { return nz_; };
    uint labels() const { return nl_; };

    virtual ~Vf3Graph() = default;
    virtual uint32_t NodeCount() const override { return n_; }
    virtual uint32_t GetNodeAttr(vflib::nodeID_t node) const override {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)

      return mVertexLabels.at(node);
    }
    virtual uint32_t OutEdgeCount(vflib::nodeID_t node) const override {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(outgoingcsr_.compressedRow[node+1] >= outgoingcsr_.compressedRow[node]);
      return outgoingcsr_.compressedRow[node+1] - outgoingcsr_.compressedRow[node];
    }
    virtual vflib::nodeID_t GetOutEdge(vflib::nodeID_t node, uint32_t i, uint32_t *pattr) const override {
      auto [neigh, edgeLabel] = GetOutEdgeWeight(node, i);
      pattr = nullptr;
      return neigh;
    }
    virtual std::pair<vflib::nodeID_t, uint32_t> GetOutEdgeWeight(vflib::nodeID_t node, uint32_t i) const override {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(i >= 0);
      assert(i < OutEdgeCount(node));
      uint index = outgoingcsr_.compressedRow[node] + i;
      assert(index < outgoingcsr_.vec.size());
      assert(index < outgoingcsr_.edgeLabel.size());
      uint neigh = outgoingcsr_.vec[index];
      uint edgeLabel = outgoingcsr_.edgeLabel[index];
      return std::make_pair(neigh, edgeLabel);
    }
    OutEdgeIterator outEdges(uint node) const {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      OutEdgeIterator it;
      it.vec = outgoingcsr_.vec.cbegin() + outgoingcsr_.compressedRow[node];
      it.vecEnd = outgoingcsr_.vec.cbegin() + outgoingcsr_.compressedRow[node+1];
      it.edgeLabel = outgoingcsr_.edgeLabel.cbegin() + outgoingcsr_.compressedRow[node];
      return it;
    }
    uint label(uint i) const {
      return GetNodeAttr(i);
    }
    uint outDegree(uint i) const {
      return OutEdgeCount(i);
    }
    uint outEdge(uint i, uint j) const {
      return GetOutEdgeWeight(i, j).first;
    }
    std::optional<uint> containsEdge(uint x, uint y) const {
      // binary search for y in outgoing edges of x
      // return edge label if found, else std::nullopt
      auto [begin, end] = neighbors(x);
      auto it = std::lower_bound(begin, end, y);
      if (it != end && *it == y) {
        auto distance = std::distance(outgoingcsr_.vec.cbegin(), it);
        auto edgeLabel = outgoingcsr_.edgeLabel[distance];
        return edgeLabel;
      }
      return std::nullopt;
    }

    uint edgeCount() const {
      uint count = 0;
      for (uint i = 0; i < n_; i++) {
        count += outDegree(i);
      }
      return count;
    }

    template<class PatternGraph>
    std::set<PatternGraph> distinctEdges() const {
      std::set<PatternGraph> patterns;
      std::set<std::tuple<uint32_t, uint32_t, uint32_t, std::optional<uint32_t>>> edgeSet;
      for (uint v1 = 0; v1 < n_; v1++) {
        auto vp1 = label(v1);
        auto [vec, vecEnd, edgeLabel] = outEdges(v1);
        for (auto dist = 0; dist < std::distance(vec, vecEnd); dist++) {
          auto v2 = vec[dist];
          auto vp2 = label(v2);
          auto elabel = edgeLabel[dist];
          auto relabel = containsEdge(v2, v1);
          if (relabel && v1 > v2) { // already done via v1 -> v2
            continue;
          }
          auto ele = std::make_tuple(vp1, vp2, elabel, relabel);
          if (edgeSet.insert(ele).second) {
            Type::LabelsVector labels;
            labels.push_back(vp1);
            labels.push_back(vp2);
            if (relabel) {
              std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edges;
              edges.push_back(std::make_tuple(0, 1, elabel));
              if (edgeSet.insert(std::make_tuple(vp1, vp2, elabel, std::nullopt)).second)
                patterns.insert(PatternGraph(labels, edges));
              edges.push_back(std::make_tuple(1, 0, *relabel));
              patterns.insert(PatternGraph(labels, edges));
              if (edgeSet.insert(std::make_tuple(vp2, vp1, *relabel, std::nullopt)).second) {
                edges.erase(edges.cbegin());
                patterns.insert(PatternGraph(labels, edges));
              }
            } else {
              std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edges;
              edges.push_back(std::make_tuple(0, 1, elabel));
              patterns.insert(PatternGraph(labels, edges));
            }
          }
        }
      }
      return patterns;
    }

    uint edgeLabels() const {
      return mEdgeLabelMap.size();
    }

  private:
    uint n_{};
    uint nz_{};
    uint nl_{};

    NeighborCSR outgoingcsr_;

    std::vector<uint> mVertexLabels;

    Type::RemapEdges mEdgeLabelMap;
    Type::RemapEdgesReverse mEdgeLabelMapReverse;

  public:
    using NeighborIterator = std::vector<uint>::const_iterator;

    std::tuple<NeighborIterator, NeighborIterator> neighbors(int node) const {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      return std::make_tuple(neighborsBegin(node), neighborsEnd(node));
    }

    NeighborIterator neighborsBegin(int node) const {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      return outgoingcsr_.vec.begin() + outgoingcsr_.compressedRow[node];
    }

    NeighborIterator neighborsEnd(int node) const {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      return outgoingcsr_.vec.begin() + outgoingcsr_.compressedRow[node+1];
    }

    friend std::ostream& operator<<(std::ostream &oss, const Vf3Graph &graph);
    friend std::ofstream& operator<<(std::ofstream &oss, const Vf3Graph &graph);
    friend std::ifstream& operator>>(std::ifstream& iss, Vf3Graph &graph);

  public:
    Vf3Graph() = default;
    Vf3Graph(DataHelpers::GraphConstructor &gC) :
      n_(gC.n), nz_(gC.nz), nl_(DataHelpers::numLabels(gC.vertexLabels)),
      outgoingcsr_(DataHelpers::processEdge(n_, gC.outgoingEdges, gC.vertexLabels)),
      mVertexLabels(gC.vertexLabels), mEdgeLabelMap(gC.remapEdgeLabels),
      mEdgeLabelMapReverse(DataHelpers::reverseMap(gC.remapEdgeLabels))
    {
      assert(outgoingcsr_.vec.size() == outgoingcsr_.edgeLabel.size());
    }
};

class SortedGraph {
  public:

    using NeighborIterator = std::vector<uint>::const_iterator;
    using IteratorPair = std::pair<NeighborIterator, NeighborIterator>;

    uint nodes() const { return n_; };
    uint nonZeros() const { return nz_; };
    uint labels() const { return nl_; };

    IteratorPair label(int l) const {
      assert(l >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(l < maxLabel_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      auto startIdx = mVertexLabels.compressedRow.cbegin() + l;
      auto endIdx = mVertexLabels.compressedRow.cbegin() + l + 1;
      assert(*startIdx >= 0);
      assert(*startIdx <= *endIdx);
      assert(*endIdx <= mVertexLabels.vec.size());
      auto start = mVertexLabels.vec.cbegin() + *startIdx;
      auto end = mVertexLabels.vec.cbegin() + *endIdx;
      return std::make_pair(start, end);
    }

    IteratorPair neighbors(uint node, uint label) const {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(label >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(label < nl_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      auto labelStart = outgoingcsr_.compressedLabels.cbegin() + node * (maxLabel_ + 1) + label;
      auto labelEnd = outgoingcsr_.compressedLabels.cbegin() + node * (maxLabel_ + 1) + label + 1;
      assert(*labelStart >= 0);
      assert(*labelStart <= *labelEnd);
      assert(*labelEnd <= outgoingcsr_.vec.size());
      auto start = outgoingcsr_.vec.cbegin() + *labelStart;
      auto end = outgoingcsr_.vec.cbegin() + *labelEnd;
      return std::make_pair(start, end);
    }
    uint edgeLabels() const {
      return mEdgeLabelMap.size();
    }
    std::optional<uint> containsEdge(uint x, uint y) const {
      assert(x >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(x < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(y >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(y < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      uint distance = std::numeric_limits<uint>::max();
      for (uint i = 0; i < maxLabel_+1; ++i) {
        auto [start, end] = neighbors(x, i);
        auto it = std::lower_bound(start, end, y);
        if (it != end && *it == y) {
          distance = std::distance(outgoingcsr_.vec.cbegin(), it);
          return outgoingcsr_.edgeLabel[distance];
        }
      }
      return std::nullopt;
    }
    NeighborLabelIterator outEdges(uint node) const {
      assert(node >= 0); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      assert(node < n_); // NOLINT(cppcoreguidelines-pro-bounds-array-to-pointer-decay)
      uint offset = outgoingcsr_.compressedLabels[node];
      NeighborLabelIterator it;
      return it;
    }

    template<class PatternGraph>
    std::set<PatternGraph> distinctEdges() const {
      std::set<PatternGraph> patterns;
      std::set<std::tuple<uint32_t, uint32_t, uint32_t, std::optional<uint32_t>>> edgeSet;
      for (uint vp1 = 0; vp1 < maxLabel_+1; ++vp1) {
        auto [start1, end1] = label(vp1);
        for (auto v1 = start1; v1 != end1; ++v1) {
          auto [start2, end2, label, edgeLabel] = outEdges(*v1);
          for (auto v2 = start2; v2 != end2; ++v2) {
            auto distance = std::distance(outgoingcsr_.vec.cbegin(), v2);
            auto vp2 = outgoingcsr_.label[distance];
            auto elabel = outgoingcsr_.edgeLabel[distance];
            auto relabel = containsEdge(*v2, *v1);
            if (relabel && v1 > v2) { // already done via v1 -> v2
              continue;
            }
            auto ele = std::make_tuple(vp1, vp2, elabel, relabel);
            if (edgeSet.insert(ele).second) {
              Type::LabelsVector labels;
              labels.push_back(vp1);
              labels.push_back(vp2);
              if (relabel) {
                std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edges;
                edges.push_back(std::make_tuple(0, 1, elabel));
                edges.push_back(std::make_tuple(1, 0, *relabel));
                patterns.insert(PatternGraph(labels, edges));
              } else {
                std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edges;
                edges.push_back(std::make_tuple(0, 1, elabel));
                patterns.insert(PatternGraph(labels, edges));
              }
            }
          }
        }
      }
      return patterns;
    }

    friend std::ostream& operator<<(std::ostream &oss, const SortedGraph &graph);
    friend std::ofstream& operator<<(std::ofstream &oss, const SortedGraph &graph);
    friend std::ifstream& operator>>(std::ifstream& iss, SortedGraph &graph);

    SortedGraph() = default;
    SortedGraph(DataHelpers::GraphConstructor &gC) :
      n_(gC.n), nz_(gC.nz), nl_(DataHelpers::numLabels(gC.vertexLabels)),
      maxLabel_(*std::max_element(gC.vertexLabels.begin(), gC.vertexLabels.end())),
      outgoingcsr_(DataHelpers::csrNeighborLabels(n_, gC.outgoingEdges, gC.vertexLabels)),
      mVertexLabels(DataHelpers::csrVertexLabels(gC.vertexLabels)),
      mEdgeLabelMap(gC.remapEdgeLabels),
      mEdgeLabelMapReverse(DataHelpers::reverseMap(gC.remapEdgeLabels))
    {
      assert(mVertexLabels.compressedRow.size() == (maxLabel_ + 1) + 1);
      assert(outgoingcsr_.compressedLabels.size() == n_ * (maxLabel_ + 1) + 1);
    }

  private:
    uint n_{};
    uint nz_{};
    uint nl_{};
    uint maxLabel_{};

    struct NeighborLabelCSR outgoingcsr_;

    VertexCSR<uint> mVertexLabels;
    Type::RemapEdges mEdgeLabelMap;
    Type::RemapEdgesReverse mEdgeLabelMapReverse;

};

template <typename T>
std::ostream& vectorPrint(std::ostream &oss, const std::string &prefix, const std::vector<T> &vec) {
    size_t size = vec.size();
    oss << prefix << " " << size << " : ";
    auto end = vec.size() > 20 ? vec.begin() + 10 : vec.end();
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(oss, " "));
    if (vec.size() > 20) { oss << "... "; }
    oss << "\n";
    return oss;
}

std::ostream& operator<<(std::ostream &oss, const Vf3Graph &graph) { // NOLINT(misc-definitions-in-headers)
  oss << "[ " << graph.n_ << " x " << graph.n_ << " ] ( " << graph.nz_ << " ) \n";
  vectorPrint(oss, "outgoing compressed ", graph.outgoingcsr_.compressedRow);
  vectorPrint(oss, "outgoing vector     ", graph.outgoingcsr_.vec);
  vectorPrint(oss, "vertex labels       ", graph.mVertexLabels);
  return oss;
}

std::ostream& operator<<(std::ostream &oss, const SortedGraph &graph) { // NOLINT(misc-definitions-in-headers)
  oss << "[ " << graph.n_ << " x " << graph.n_ << " ] ( " << graph.nz_ << " ) \n";
  vectorPrint(oss, "outgoing neigh with labels", graph.outgoingcsr_.compressedLabels);
  vectorPrint(oss, "outgoing vector     ", graph.outgoingcsr_.vec);
  vectorPrint(oss, "vertex labels compressed       ", graph.mVertexLabels.compressedRow);
  vectorPrint(oss, "vertex labels       ", graph.mVertexLabels.vec);

  return oss;
}

namespace graphStream {
  template <typename T>
    std::ofstream& serialize(std::ofstream &oss, const std::vector<T> &vec) { // NOLINT(misc-definitions-in-headers)
      size_t size = vec.size();
      oss.write(reinterpret_cast<const char*>(&size), sizeof(size));
      oss.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(T));
      return oss;
    }

  template <typename T>
    std::ifstream& deserialize(std::ifstream &iss, std::vector<T> &vec) { // NOLINT(misc-definitions-in-headers)
      size_t size;
      iss.read(reinterpret_cast<char*>(&size), sizeof(size));
      vec.resize(size);
      iss.read(reinterpret_cast<char*>(vec.data()), size * sizeof(T));
      return iss;
    }

  template <template<typename, typename> typename Map, typename T, typename U> // Container can be map or unordered_map
    std::ifstream& deserialize(std::ifstream &iss, Map<T, U> &map) { // NOLINT(misc-definitions-in-headers)
      std::vector<T> first;
      std::vector<U> second;
      deserialize(iss, first);
      deserialize(iss, second);
      map = DataHelpers::pack<Map>(first, second);
      return iss;
    }

  template <template<typename, typename> typename Map, typename T, typename U> // Container can be map or unordered_map
    std::ofstream& serialize(std::ofstream &oss, const Map<T, U> &map) { // NOLINT(misc-definitions-in-headers)
      std::vector<T> first;
      std::vector<U> second;
      std::tie(first, second) = DataHelpers::unpack(map);
      serialize(oss, first);
      serialize(oss, second);
      return oss;
    }

  template <typename T>
    std::ofstream& serialize(std::ofstream &oss, const T &data) { // NOLINT(misc-definitions-in-headers)
      oss.write(reinterpret_cast<const char*>(&data), sizeof(data));
      return oss;
    }

  template <typename T>
    std::ifstream& deserialize(std::ifstream &iss, T &data) { // NOLINT(misc-definitions-in-headers)
      iss.read(reinterpret_cast<char*>(&data), sizeof(data));
      return iss;
    }
}; // namespace graphStream

std::ofstream& operator<<(std::ofstream &oss, const Vf3Graph &graph) { // NOLINT(misc-definitions-in-headers)
  assert(graph.n_ == graph.outgoingcsr_.compressedRow.size() - 1);
  assert(graph.outgoingcsr_.vec.size() == graph.outgoingcsr_.edgeLabel.size());

  std::set<uint> edgeLabels(graph.outgoingcsr_.edgeLabel.begin(), graph.outgoingcsr_.edgeLabel.end());
  assert(edgeLabels.size() == graph.mEdgeLabelMap.size());

  graphStream::serialize(oss, graph.n_);
  graphStream::serialize(oss, graph.nz_);
  graphStream::serialize(oss, graph.nl_);

  graphStream::serialize(oss, graph.outgoingcsr_.compressedRow);
  graphStream::serialize(oss, graph.outgoingcsr_.vec);
  graphStream::serialize(oss, graph.outgoingcsr_.edgeLabel);

  graphStream::serialize(oss, graph.mVertexLabels);

  graphStream::serialize(oss, graph.mEdgeLabelMap);
  graphStream::serialize(oss, graph.mEdgeLabelMapReverse);

  return oss;
}

std::ifstream& operator>>(std::ifstream& iss, Vf3Graph &graph) { // NOLINT(misc-definitions-in-headers)
  graphStream::deserialize(iss, graph.n_);
  graphStream::deserialize(iss, graph.nz_);
  graphStream::deserialize(iss, graph.nl_);

  graphStream::deserialize(iss, graph.outgoingcsr_.compressedRow);
  graphStream::deserialize(iss, graph.outgoingcsr_.vec);
  graphStream::deserialize(iss, graph.outgoingcsr_.edgeLabel);

  assert(graph.outgoingcsr_.compressedRow.size() == graph.n_ + 1);
  assert(graph.outgoingcsr_.vec.size() == graph.outgoingcsr_.edgeLabel.size());

  graphStream::deserialize(iss, graph.mVertexLabels);

  graphStream::deserialize(iss, graph.mEdgeLabelMap);
  graphStream::deserialize(iss, graph.mEdgeLabelMapReverse);

  std::set<uint> edgeLabels(graph.outgoingcsr_.edgeLabel.begin(), graph.outgoingcsr_.edgeLabel.end());
  assert(edgeLabels.size() == graph.mEdgeLabelMap.size());
  assert(edgeLabels.size() == graph.mEdgeLabelMapReverse.size());

  return iss;
}

std::ofstream& operator<<(std::ofstream& oss, const SortedGraph &graph) { // NOLINT(misc-definitions-in-headers)

  graphStream::serialize(oss, graph.n_);
  graphStream::serialize(oss, graph.nz_);
  graphStream::serialize(oss, graph.nl_);
  graphStream::serialize(oss, graph.maxLabel_);

  assert(graph.outgoingcsr_.compressedLabels.size() == (graph.n_ * (graph.maxLabel_ + 1))+1);

  graphStream::serialize(oss, graph.outgoingcsr_.compressedLabels);
  graphStream::serialize(oss, graph.outgoingcsr_.vec);
  graphStream::serialize(oss, graph.outgoingcsr_.label);
  graphStream::serialize(oss, graph.outgoingcsr_.edgeLabel);

  assert(graph.mVertexLabels.compressedRow.size() == graph.nl_ + 1);

  graphStream::serialize(oss, graph.mVertexLabels.compressedRow);
  graphStream::serialize(oss, graph.mVertexLabels.vec);

  graphStream::serialize(oss, graph.mEdgeLabelMap);
  graphStream::serialize(oss, graph.mEdgeLabelMapReverse);

  return oss;
}

std::ifstream& operator>>(std::ifstream& iss, SortedGraph &graph) { // NOLINT(misc-definitions-in-headers)

  graphStream::deserialize(iss, graph.n_);
  graphStream::deserialize(iss, graph.nz_);
  graphStream::deserialize(iss, graph.nl_);
  graphStream::deserialize(iss, graph.maxLabel_);

  graphStream::deserialize(iss, graph.outgoingcsr_.compressedLabels);
  graphStream::deserialize(iss, graph.outgoingcsr_.vec);
  graphStream::deserialize(iss, graph.outgoingcsr_.label);
  graphStream::deserialize(iss, graph.outgoingcsr_.edgeLabel);

  graphStream::deserialize(iss, graph.mVertexLabels.compressedRow);
  graphStream::deserialize(iss, graph.mVertexLabels.vec);

  graphStream::deserialize(iss, graph.mEdgeLabelMap);
  graphStream::deserialize(iss, graph.mEdgeLabelMapReverse);

  return iss;
}

template <typename DataGraph>
DataGraph loadBinGraph(const std::string &filename) {
  // if DataGraph is Vf3Graph, then filename should end with .vf3.bin
  // if DataGraph is SortedGraph, then filename should end with .srt.bin
  if (std::is_same<DataGraph, Vf3Graph>::value) {
    assert(filename.substr(filename.size() - 8) == ".vf3.bin");
  } else if (std::is_same<DataGraph, SortedGraph>::value) {
    assert(filename.substr(filename.size() - 8) == ".srt.bin");
  } else {
    std::cerr << "Unknown graph type" << std::endl;
    assert(false);
    std::exit(1);
  }

  std::ifstream fstream(filename, std::ios::binary);
  if (!fstream.is_open()) {
    std::cerr << "Could not open file " << filename << std::endl;
    assert(false);
    std::exit(1);
  }
  DataGraph graph;
  fstream >> graph;
  return graph;
}

