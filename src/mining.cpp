#include "mining.hpp"

#include "pattern.hpp"
#include "matcher.hpp"
#include "generator.hpp"

#include <CLI/CLI.hpp>
#include <fmt/core.h>

#include <chrono>
#include <iostream>
#include <optional>
#include <filesystem>
#include <nlohmann/json.hpp>

template<auto Func, typename... Args>
auto timeWrapper(Args &&... args) {
  auto t1 = std::chrono::high_resolution_clock::now();
  auto res = Func(std::forward<Args>(args)...);
  auto t2 = std::chrono::high_resolution_clock::now();
  return std::make_pair(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1), res);
}

template<typename PatternGraph, typename DataGraph>
auto timeWrapperEdges(const DataGraph &dataGraph) {
  auto t1 = std::chrono::high_resolution_clock::now();
  auto res = dataGraph.template distinctEdges<PatternGraph>();
  auto t2 = std::chrono::high_resolution_clock::now();
  return std::make_pair(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1), res);
}

MetricsVec runProgram(const Arguments &args) {
  using namespace std::chrono_literals;
#if D_CUSTOM_FREQUENCY
  auto dataGraph = loadBinGraph<SortedGraph>(args.mFilename);
#else
  auto dataGraph = loadBinGraph<Vf3Graph>(args.mFilename);
  auto dataGraphStruct = Vf3Frequency::DataGraphStruct(dataGraph);
#endif
  auto [initialGeneration, patternVec] = timeWrapperEdges<Vf3Pattern>(dataGraph);
  //auto [initialGeneration, patternVec] = timeWrapper<PatternGenerator::edges<Vf3Pattern>>(dataGraph.labels(), dataGraph.edgeLabels());
  if (args.mVerbose || args.mPrintOutput) flushedPrint("Initial generation took {} ms, generated {} patterns\n", initialGeneration.count(), patternVec.size());
  if (args.mPrintPatterns) {
    std::cout << "Merged patterns: ";
    std::cout << patternVec;
  }
  int step = 2;
  MetricsVec metricsVec;
  metricsVec.emplace_back(Metrics{0ms, initialGeneration, patternVec.size(), 0, dataGraph.edgeCount()});
  while(step <= args.mMaxSize && !patternVec.empty()) {
    if (args.mVerbose || args.mPrintOutput) flushedPrint("Step {}: starting with {} patterns\n", step, patternVec.size());
#ifdef D_CUSTOM_FREQUENCY
    auto [stepMining, freqPatterns] = timeWrapper<customFrequency::matchThres<0, 0, SortedGraph, Vf3Pattern>>(dataGraph, patternVec, args.mSupport, args.mAccuracy);
#else
    auto [stepMining, freqPatterns] = timeWrapper<Vf3Frequency::matchThres<0, 0, Vf3Graph, Vf3Pattern>>(dataGraphStruct, patternVec, args.mSupport, args.mAccuracy);
#endif
    if (args.mVerbose || args.mPrintOutput) flushedPrint("Matching took {} ms, reduced {} patterns to {} frequent patterns\n", stepMining.count(), patternVec.size(), freqPatterns.size());
    if (args.mPrintPatterns) {
      std::cout << "Frequent patterns: ";
      std::cout << freqPatterns;
    }
    if (step == args.mMaxSize)  {
      metricsVec.emplace_back(Metrics{stepMining, 0ms, freqPatterns.size(), 0, 0});
      break;
    }
    auto [stepGeneration, combinedPair] = timeWrapper<PatternGenerator::extend<Vf3Pattern>>(freqPatterns);
    auto [newPatternVecs, prunedPatternsCount] = combinedPair;
    if (args.mVerbose || args.mPrintOutput) flushedPrint("Generation took {} ms, extended {} patterns to {} new patterns\n", stepGeneration.count(), freqPatterns.size(), newPatternVecs.size());
    if (args.mPrintPatterns) {
      std::cout << "Merged patterns: ";
      std::cout << newPatternVecs;
    }
    metricsVec.emplace_back(Metrics{stepMining, stepGeneration, freqPatterns.size(), newPatternVecs.size(), prunedPatternsCount});
    std::swap(patternVec, newPatternVecs);
    step++;
  }
  if (args.mVerbose || args.mPrintOutput) flushedPrint("Finished mining {} steps out of {} steps, {} patterns left\n", step - 1, args.mMaxSize, patternVec.size());
  return metricsVec;
}

nlohmann::json metricsToJson(const MetricsVec &metricsVec) {
  nlohmann::json json;
  for (int i = 0; i < metricsVec.size(); i++) {
    auto [miningTime, generationTime, patternCount, freqPatternCount, prunedPatternsCount] = metricsVec[i];
    json.push_back({
        {"level", i+2},
        {"miningTime (ms)", miningTime.count()},
        {"generationTime (ms)", generationTime.count()},
        {"generatedPatterns", patternCount},
        {"freqPatternCount", freqPatternCount},
        {"prunedPatterns", prunedPatternsCount}
        });
  }
  return json;
}

#ifndef NO_MAIN

std::optional<Arguments> commandParser(int argc, char *argv[]) {
  CLI::App app{"flexis"};

  std::string filename;
  std::string jsonFilename;
  int support(10);
  int maxSize(6);
  float accuracy(0.98);
  bool verbose{false};
  bool printOutput{false};
  bool outputJson{false};
  bool printPatterns{false};
  app.add_option("--input", filename, "Input bin file")->check(CLI::ExistingFile)->required();
  app.add_option("--output", jsonFilename, "Output json file")->default_val("");
  app.add_option("-s,--support", support, "support")->default_val(support);
  app.add_option("-k,--max-size", maxSize, "max size of candidate pattern")->check(CLI::Range(2, 10))->default_val(maxSize);
  app.add_option("-a,--slider", accuracy, "slider")->check(CLI::Range(0.0, 1.0))->default_val(accuracy);
  app.add_flag("--verbose", verbose, "verbose mode")->default_val(verbose);
  app.add_flag("--print-output", printOutput, "print output")->default_val(printOutput);
  app.add_flag("--output-json", outputJson, "output json")->default_val(outputJson);
  app.add_flag("--print-patterns", printPatterns, "print patterns")->default_val(printPatterns); 

  bool skipExec = false;
  try{
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    if (app.exit(e) != 0) {
      throw e;
    }
    return {};
  }

  using namespace std::filesystem;
  if (outputJson && jsonFilename.empty()) {
    auto datagraphName = path(filename).stem().stem().string();
    if (printPatterns)
      fmt::print("fileVersion: 3\n");
    else
      fmt::print("fileVersion: 1\n");
    fmt::print("datagraphName: {}\n", datagraphName);
#ifdef D_CUSTOM_FREQUENCY
    std::string type = "custom";
#else
    std::string type = "vf3";
#endif
    // change accuracy to 4 digits
    jsonFilename = path("timings") / path(fmt::format("{}-{}-{:4f}-{}-{}.time.json", datagraphName, support, accuracy, type, maxSize));
  }
  if (verbose) {
    std::cout << "Filename " << filename << std::endl;
  }
  return { { filename, jsonFilename, support, maxSize, accuracy, verbose, printOutput, printPatterns } };
}

int main(int argc, char *argv[]) {
  auto args = commandParser(argc, argv);
  if (!args) return 0;
  if (args->mVerbose) {
    std::cout << "threshold is set to " << args->mSupport << std::endl;
  }
  auto res = runProgram(args.value());
  if (!args->mJsonFilename.empty()) {
    auto json = metricsToJson(res);
    nlohmann::json jsonData;
    jsonData["support"] = args->mSupport;
    jsonData["accuracy"] = args->mAccuracy;
    jsonData["maxSize"] = args->mMaxSize;
#ifdef D_CUSTOM_FREQUENCY
    jsonData["type"] = "custom";
#else
    jsonData["type"] = "vf3";
#endif
    jsonData["data"] = json;
    std::ofstream jsonFile(args->mJsonFilename);
    jsonFile << jsonData.dump(2) << std::endl;
    jsonFile.close();
  }
  return 0;
}

#endif
