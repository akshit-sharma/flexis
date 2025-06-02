#include "process.hpp"
#include "graph.hpp"

#include <CLI/CLI.hpp>

#include <chrono>
#include <fmt/core.h>
#include <string>
#include <filesystem>
#include <tuple>

#include <indicators/progress_bar.hpp>
#include <indicators//cursor_control.hpp>

Arguments commandParser(int argc, char *argv[]) {
  CLI::App app{"Preprocess Graph"};

  std::filesystem::path mtxFilename;
  std::filesystem::path binFilename;
  std::filesystem::path dotDirectory;
  int seed;
  Format kFormat;
  bool verbose;
  bool progress;
  bool printVf3Graph;
  bool printSortedGraph;
  std::map<std::string, Format> formatMat{{"mtx", Format::KMtx}, {"dimacs", Format::KDimacs}, {"citPatentsTxt", Format::KCitPatentsTxt}, {"custom", Format::KCustom}, {"lg", Format::KLG}};
  app.add_option("-i,--input", mtxFilename, "Input mtx file")->required()->check(CLI::ExistingFile);
  app.add_option("-o,--output", binFilename, "Output bin prefix (without extension)")->required();
  app.add_flag("--verbose", verbose, "verbose mode")->default_val(false);
  app.add_flag("--progress", progress, "show progress")->default_val(true);
  app.add_option("--format", kFormat, "file format")->default_val(Format::KCustom)
    ->transform(CLI::CheckedTransformer(formatMat, CLI::ignore_case));
  app.add_option("--seed", seed, "seed for random number generator")->default_val(0);
  app.add_flag("--print-vf3-graph", printVf3Graph, "print vf3 graph")->default_val(false);
  app.add_flag("--print-sorted-graph", printSortedGraph, "print sorted graph")->default_val(false);
  app.add_option("--dot-files-output", dotDirectory, "dot files output")->default_val("./data/actual");

  binFilename = binFilename.stem();

  bool skipExec = false;
  try{
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    skipExec = true;
    if (app.exit(e) != 0) {
      abort();
    }
  }

  return {mtxFilename, binFilename, dotDirectory, seed, kFormat, verbose, progress, printVf3Graph, printSortedGraph, skipExec};
}

DataHelpers::GraphConstructor loadGraph(const std::string &filename,
                Format format, int seed) {
  if (format == Format::KCitPatentsTxt)
    return loadKCitPatentsTxt(filename, seed);
  if (format == Format::KMtx)
    return loadMtxGraph(filename, seed);
  if (format == Format::KCustom)
    return loadKCustomTxt(filename, seed);
  if (format == Format::KLG)
    return loadLGGraph(filename, seed);
  if (format != Format::KDimacs)
    throw std::runtime_error("Unknown format");
  return loadDimacsGraph(filename, seed);
}

template<typename DataGraph>
DataGraph saveGraph(DataHelpers::GraphConstructor &graphContructor,
               const std::string &outputFilename) {
  DataGraph graph(graphContructor);
  std::ofstream out(outputFilename.c_str(), std::ios::binary);
  out << graph;
  out.close();
  return graph;
}
/*
template<typename DataGraph>
void saveDotGraph(DataGraph graph, std::ofstream &stream) {

}

template<typename DataGraph>
std::string saveDotGraph(DataGraph graph, const std::filesystem::path &outDotFolder,
    const std::string &inputFilename) {
  std::filesystem::path inputFilepath(inputFilename);
  if (inputFilepath.extension() != ".lg") {
    return "";
  }
  std::filesystem::path outDotPath = outDotFolder / inputFilepath.filename();
  std::ofstream stream(outDotPath);
  saveDotGraph(graph, stream);
  return outDotPath;
}
*/
int main(int argc, char *argv[]) {
  auto t1 = std::chrono::high_resolution_clock::now();

  auto args = commandParser(argc, argv);

  if (args.mSkipExec) return 236;

  using namespace indicators;
  show_console_cursor(false);

  ProgressBar bar {
    option::BarWidth{50},
    option::Start{"["},
    option::Fill{"■"},
    option::Lead{"■"},
    option::Remainder{"-"},
    option::End{" ]"},
    option::PostfixText{"Loading "+args.mInputFile+" 1/5"},
    option::ForegroundColor{Color::cyan},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}}
  };

  if (args.mProgress) {
    bar.set_progress(0);
  }

  auto graphContructor = loadGraph(args.mInputFile, args.mFormat, args.mSeed);

  if (args.mProgress) {
    bar.set_option(option::PostfixText{"saving "+args.mInputFile+".vf3.bin 2/5"});
  }

  auto vf3Graph = saveGraph<Vf3Graph>(graphContructor, args.mOutputFile + ".vf3.bin");

  std::cout << " #V_labels: " << vf3Graph.labels()
    << " | #E_labels: " << vf3Graph.edgeLabels() << " " << std::endl;

  if (args.mPrintVf3Graph) {
    std::cout << vf3Graph << std::endl;
  }

  if (args.mProgress) {
    bar.set_option(option::PostfixText{"saving "+args.mInputFile+".srt.bin 3/5"});
  }

  auto sortedGraph = saveGraph<SortedGraph>(graphContructor, args.mOutputFile + ".srt.bin");

  if (args.mPrintSortedGraph) {
    std::cout << sortedGraph << std::endl;
  }

 // auto dotFilepath = saveDotGraph(vf3Graph, args.mDotFolder, args.mInputFile);

  //if (args.mProgress) {
 //   bar.set_option(option::PostfixText{"saving "+dotFilepath+" 4/5"});
  //}

  auto t2 = std::chrono::high_resolution_clock::now();
  auto diff = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();

  auto timeSecs = static_cast<float>(diff)/1000000;

  if (args.mProgress) {
    bar.set_option(option::PostfixText{"time taken "+std::to_string(timeSecs)+" seconds 5/5"});
  }

  if (args.mVerbose) {
    std::cout <<"Time taken " << timeSecs << " seconds" << std::endl;
  }

  show_console_cursor(true);

  return 0;
}

