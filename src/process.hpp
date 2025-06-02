#pragma once

#include <algorithm>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <unordered_map>
#include <iterator>
#include <random>
#include <map>
#include <set>

#include <cassert>
#include <cstdlib>

#include "graph.hpp"

// code below is used to load .mtx file and convert it to Graph object.

void expectValue(const std::string &value, const std::vector<std::string> &expectVec) { // NOLINT(misc-definitions-in-headers)
  for (const auto& expect : expectVec) {
    if (value == expect) {
      return;
    }
  }
  std::cerr << "Expected " << value << " to be ";
  bool first = true;
  for (const auto &expect : expectVec) {
    if (first) {
      std::cerr << expect;
    } else {
      std::cerr << " or " << expect;
    }
    first = false;
  }
  std::cerr << "\n";
  std::abort();
}

bool readFirstLine(const std::string &line) { // NOLINT(misc-definitions-in-headers)
  std::string first, object, format, field, symmetry;  // NOLINT(readability-isolate-declaration)
  std::istringstream iss(line);
  if (!(iss >> first >> object >> format >> field >> symmetry)) {
    std::cerr << "Error getting first line, got " << first << ", " << object
              << ", " << format << ", " << field << ", " << symmetry << std::endl;
    std::abort();
  }
  expectValue(first, { "%%MatrixMarket", "MatrixMarket" });
  expectValue(object,  { "matrix" });
  expectValue(format, { "array", "coordinate" });
  expectValue(field, { "real", "double", "integer", "pattern" });
  expectValue(symmetry, {"general", "symmetric"});
  return symmetry == "symmetry";
}

std::vector<uint> randomLabels(int pSize, int pLabels, int pSeed) { // NOLINT(misc-definitions-in-headers)
  std::vector<uint> labels(pSize);
  std::mt19937 gen(pSeed);
  std::uniform_int_distribution<> dis(0, pLabels - 1);
  std::generate(labels.begin(), labels.end(), [&] { return dis(gen); });

  return labels;
}



DataHelpers::GraphConstructor loadMtxGraph(const std::string &pFilename, int seed) { // NOLINT(misc-definitions-in-headers)
  std::ifstream fstream(pFilename);
  if (!fstream.good()) {
    std::cerr << "Error opening " << pFilename << std::endl;
    std::abort();
  }
  std::string line;
  if (!std::getline(fstream, line)) {
    std::cerr << "error reading first line" << std::endl;
    std::abort();
  }

  bool isSymmetric = readFirstLine(line);
  while(std::getline(fstream, line)) {
    if (line[0] != '%') break;
  }

  uint r = 0, c = 0, nz = 0;  // NOLINT(readability-isolate-declaration)
  std::istringstream iss(line);
  if (!(iss >> r >> c >> nz)) {
    std::cerr << "Error getting dimensions of matrix, line("<<line<<")" << std::endl;
    std::abort();
  }
  if (r != c) {
    std::cerr << "Expected matrix to be square but got " << r << "x" << c << std::endl;
    std::abort();
  }

  std::unordered_map<uint, std::map<uint, uint>> outgoingEdges;
  std::unordered_map<uint, std::map<uint, uint>> incomingEdges;

  std::unordered_map<uint, uint> edgeLabelRemap { { kDefaultEdgeLabel, 0 } };

  int actualLines = 0;
  while(std::getline(fstream, line)) {
    if (line[0] == '%') continue;
    if (line.empty()) continue;

    int m = 0, n = 0; // NOLINT(readability-isolate-declaration)
    std::istringstream iss(line);
    if (!(iss >> m >> n )) {
      std::cerr << "Error reading line " << line << std::endl;
      std::abort();
    }
    actualLines++;
    if (isSymmetric && m != n) {
      outgoingEdges[n-1].insert(std::make_pair(m-1, kDefaultEdgeLabel));
      incomingEdges[m-1].insert(std::make_pair(n-1, kDefaultEdgeLabel));
    }
    outgoingEdges[m-1].insert(std::make_pair(n-1, kDefaultEdgeLabel));
    incomingEdges[n-1].insert(std::make_pair(m-1, kDefaultEdgeLabel));
  }
  if (actualLines != nz) {
    std::cerr << "nonzeros " << nz << " does not match acutalLines " << actualLines << std::endl;
    std::abort();
  }
  auto labels = randomLabels(r, 10, seed);

  auto min = std::min_element(labels.begin(), labels.end())[0];
  if (min > 0) {
    std::transform(labels.begin(), labels.end(), labels.begin(), [min](uint label) { return label - min; });
  }
  return { r, nz, outgoingEdges, incomingEdges, labels, edgeLabelRemap };
}

DataHelpers::GraphConstructor loadDimacsGraph(const std::string &pFilename, int seed) { // NOLINT(misc-definitions-in-headers)
  std::ifstream fstream(pFilename);
  if (!fstream.good()) {
    std::cerr << "Error opening " << pFilename << std::endl;
    std::abort();
  }
  std::string line;
  do {
    if (!std::getline(fstream, line)) {
      std::cerr << "error reading first line" << std::endl;
      std::abort();
    }
  } while(line[0] == 'c');
  uint n = 0, edges = 0; // NOLINT
  if (line[0] == 'p') {
    std::string f, s, t, fo; // NOLINT
    std::istringstream iss(line);
    if (!(iss >> f >> s >> t >> fo )) {
      std::cerr << "Error reading line " << line << std::endl;
      std::abort();
    }
    n = std::atoi(t.c_str()); // NOLINT
    edges = std::atoi(fo.c_str()); // NOLINT
  }
  std::unordered_map<uint, std::map<uint,uint>> outgoingEdges;
  std::unordered_map<uint, std::map<uint,uint>> incomingEdges;
  std::unordered_map<uint, uint> edgeLabelRemap { { kDefaultEdgeLabel, 0 } };
  int actualLines = 0;
  while (std::getline(fstream, line)) {
    if (line[0] != 'a') continue;
    if (line.empty()) continue;

    actualLines++;
    int m = 0, n = 0; // NOLINT(readability-isolate-declaration)
    std::string f, s, t; // NOLINT
    std::istringstream iss(line);
    if (!(iss >> f >> m >> n)) {
      std::cerr << "Error reading line " << line << std::endl;
      std::abort();
    }
    outgoingEdges[m-1].insert(std::make_pair(n-1, kDefaultEdgeLabel));
    incomingEdges[n-1].insert(std::make_pair(m-1, kDefaultEdgeLabel));
  }
  auto labels = randomLabels(n, 10, seed);

  auto min = std::min_element(labels.begin(), labels.end())[0];
  if (min > 0) {
    std::transform(labels.begin(), labels.end(), labels.begin(), [min](uint label) { return label - min; });
  }
  return { n, edges, outgoingEdges, incomingEdges, labels, edgeLabelRemap };
}

DataHelpers::GraphConstructor loadKCitPatentsTxt(const std::string &pFilename, int seed) { // NOLINT(misc-definitions-in-headers)
  std::ifstream fstream(pFilename);
  if (!fstream.good()) {
    std::cerr << "Error opening " << pFilename << std::endl;
    std::abort();
  }
  std::string line;

  std::unordered_map<uint, std::map<uint,uint>> outgoingEdges;
  std::unordered_map<uint, std::map<uint,uint>> incomingEdges;

  std::map<int, int> nodeMap;
  uint nodeIdx = 0;
  uint numEdges = 0;
  std::unordered_map<uint, uint> edgeLabelRemap { { kDefaultEdgeLabel, 0 } };
  while(std::getline(fstream, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    int m = 0, n = 0; // NOLINT
    std::istringstream iss(line);
    if (!(iss >> m >> n )) {
      std::cerr << "Error reading line " << line << std::endl;
      std::abort();
    }
    numEdges++;
    m--;
    n--;
    // if not already in nodeMap when add nodeMap with nodeIdx value
    if (nodeMap.find(m) == nodeMap.end()) {
      nodeMap[m] = nodeIdx++;
    }
    if (nodeMap.find(n) == nodeMap.end()) {
      nodeMap[n] = nodeIdx++;
    }
    outgoingEdges[nodeMap[m]].insert(std::make_pair(nodeMap[n], kDefaultEdgeLabel));
    incomingEdges[nodeMap[n]].insert(std::make_pair(nodeMap[m], kDefaultEdgeLabel));
  }
  //assert(numEdges == edges);
  //assert(nodeIdx == nodes);
  assert(nodeMap.size() == nodeIdx);
  assert(outgoingEdges.size() < nodeIdx);
  assert(incomingEdges.size() < nodeIdx);

  auto labels = randomLabels(nodeIdx, 10, seed);

  auto min = std::min_element(labels.begin(), labels.end())[0];
  if (min > 0) {
    std::transform(labels.begin(), labels.end(), labels.begin(), [min](uint label) { return label - min; });
  }
  return { nodeIdx, numEdges, outgoingEdges, incomingEdges, labels, edgeLabelRemap };
}

DataHelpers::GraphConstructor loadKCustomTxt(const std::string &pFilename, int seed) { // NOLINT(misc-definitions-in-headers)
  std::ifstream fstream(pFilename);
  if (!fstream.good()) {
    std::cerr << "Error opening " << pFilename << std::endl;
    std::abort();
  }
  std::string line;
  if (!std::getline(fstream, line)) {
    std::cerr << "error reading first line" << std::endl;
    std::abort();
  }
  if (line.compare(0, 8, "# Nodes:") != 0) {
    std::cerr << "error reading third line (no '# Nodes' found) " << line << std::endl;
    std::abort();
  }
  uint nodes = 0, edges = 0; // NOLINT
  std::istringstream iss(line);
  std::string dump;
  if (!(iss >> dump >> dump >> nodes >> dump >> edges)) {
    std::cerr << "Error reading line for nodes and edges: " << line << std::endl;
    std::abort();
  }
  bool first = false;
  std::vector<uint> labels;
  labels.reserve(nodes);
  std::unordered_map<uint, std::map<uint,uint>> outgoingEdges;
  std::unordered_map<uint, std::map<uint,uint>> incomingEdges;

  int numEdges = 0;
  std::set<uint> nodeSet;
  std::unordered_map<uint, uint> edgeLabelRemap { { kDefaultEdgeLabel, 0 } };
  while(std::getline(fstream, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    if (!first) {
      std::istringstream iss(line);
      for (int i = 0; i < nodes; ++i) {
        int label = 0;
        if (!(iss >> label)) {
          std::cerr << "Error reading line " << line << std::endl;
          std::abort();
        }
        labels.push_back(label);
      }
      first = true;
      continue;
    }

    int m = 0, n = 0; // NOLINT
    std::istringstream iss(line);
    if (!(iss >> m >> n )) {
      std::cerr << "Error reading line " << line << std::endl;
      std::abort();
    }
    nodeSet.insert(m);
    nodeSet.insert(n);
    numEdges++;
    outgoingEdges[m].insert(std::make_pair(n, kDefaultEdgeLabel));
    incomingEdges[n].insert(std::make_pair(m, kDefaultEdgeLabel));
  }
  assert(nodeSet.size() == nodes);
  assert(labels.size() == nodes);
  assert(numEdges == edges);
  assert(outgoingEdges.size() < nodes);
  assert(incomingEdges.size() < nodes);

  auto min = std::min_element(labels.begin(), labels.end())[0];
  if (min > 0) {
    std::transform(labels.begin(), labels.end(), labels.begin(), [min](uint label) { return label - min; });
  }

  return { nodes, edges, outgoingEdges, incomingEdges, labels, edgeLabelRemap};
}

DataHelpers::GraphConstructor loadLGGraph(const std::string &pFilename, int seed) { // NOLINT(misc-definitions-in-headers)
  std::ifstream fstream(pFilename);
  if (!fstream.good()) {
    std::cerr << "Error opening " << pFilename << std::endl;
    std::abort();
  }
  std::string line;
  std::vector<uint> labels;
  std::unordered_map<uint, std::map<uint,uint>> outgoingEdges;
  std::unordered_map<uint, std::map<uint,uint>> incomingEdges;

  std::unordered_map<uint, uint> edgeLabelRemap;
  auto edgeRemap = [&edgeLabelRemap] (uint label) {
    if (edgeLabelRemap.find(label) == edgeLabelRemap.end()) {
      edgeLabelRemap[label] = edgeLabelRemap.size();
    }
    return edgeLabelRemap[label];
  };

  uint nodes = 0, edges = 0; // NOLINT
  while (std::getline(fstream, line)) {
    if (line[0] == '#') continue;
    if (line.empty()) continue;

    if (line[0] == 't') continue;

    if (line[0] == 'v') {
      std::istringstream iss(line);
      std::string dump;
      int node = 0, label = 0; // NOLINT
      if (!(iss >> dump >> node >> label)) {
        std::cerr << "Error reading line " << line << std::endl;
        std::abort();
      }
      if (node != labels.size()) {
        std::cerr << "Error reading line " << line << std::endl;
        std::abort();
      }
      nodes++;
      labels.push_back(label);
      continue;
    } else if (line[0] == 'e') {
      std::istringstream iss(line);
      std::string dump;
      uint m = 0, n = 0; // NOLINT
      uint weight = 0; // NOLINT
      if (!(iss >> dump >> m >> n >> weight)) {
        std::cerr << "Error reading line " << line << std::endl;
        std::abort();
      }
      assert(dump == "e");
      uint newWeight = edgeRemap(weight);
      edges++;
      outgoingEdges[m].insert(std::make_pair(n, newWeight));
      incomingEdges[n].insert(std::make_pair(m, newWeight));
      continue;
    }

    throw std::runtime_error("Unknown line type for line: " + line);
  }

  auto min = std::min_element(labels.begin(), labels.end())[0];
  if (min > 0) {
    std::transform(labels.begin(), labels.end(), labels.begin(), [min](uint label) { return label - min; });
  }

  return { nodes, edges, outgoingEdges, incomingEdges, labels, edgeLabelRemap };

}

enum class Format {
  KMtx, KDimacs, KCitPatentsTxt, KCustom, KLG
};

struct Arguments {
  std::string mInputFile;
  std::string mOutputFile;
  std::filesystem::path mDotFolder;
  int mSeed;
  Format mFormat{Format::KMtx};
  bool mVerbose{false};
  bool mProgress{false};
  bool mPrintVf3Graph{false};
  bool mPrintSortedGraph{false};
  bool mSkipExec{false};
  bool mMaxDegree{false};
};
