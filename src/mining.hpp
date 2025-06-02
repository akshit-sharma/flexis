#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <map>
#include <set>
#include <numeric>

#include <cassert>
#include <cstdlib>

#include "graph.hpp"

#include <filesystem>

struct Arguments {
  std::string mFilename;
  std::filesystem::path mJsonFilename;
  int mSupport;
  int mMaxSize;
  float mAccuracy;
  bool mVerbose{false};
  bool mPrintOutput{false};
  bool mPrintPatterns{false};
};

struct Metrics {
  std::chrono::milliseconds mStepMining;
  std::chrono::milliseconds mStepGeneration;
  size_t mStepPatternCount;
  size_t mStepFrequentPatternCount;
  size_t mStepDuplicateCount;
};

using MetricsVec = std::vector<Metrics>;

MetricsVec runProgram(const Arguments &args);
