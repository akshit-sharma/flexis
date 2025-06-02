#include "graph.hpp"

#include "colorScheme.hpp"

#include <CLI/CLI.hpp>

#include <chrono>
#include <cstdint>
#include <optional>
#include <string>
#include <tuple>

#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

#include <nlohmann/json.hpp>

using DataGraphClass = Vf3Graph;

struct Arguments {
  std::filesystem::path mBinFile;
  DataGraphClass mGraph;
  std::filesystem::path mOutputFile;
  std::filesystem::path mGraphFile;
  std::string mColorScheme;
  bool mVerbose;
  bool mProgress;
  Arguments(std::string graphFile, bool outputFile, std::filesystem::path dotDirectory, bool verbose, bool progress)
    : mBinFile(graphFile),
      mGraph(loadBinGraph<DataGraphClass>(graphFile)),
      mOutputFile(""),
      mGraphFile(""),
      mColorScheme("pastel19"), // set312, set28
      mVerbose(verbose),
      mProgress(progress)
  {
    using std::filesystem::path;
    path filename = path(graphFile).filename();
    if (outputFile){
      mOutputFile = path("./data/info");
      mOutputFile = mOutputFile / filename.filename();
      mOutputFile.replace_extension(""); // remove first extension
      mOutputFile.replace_extension(".info.json");
    }
    if (filename.filename().string().find("test") != std::string::npos) {
      mGraphFile = dotDirectory / filename.filename();
      mGraphFile.replace_extension(""); // remove first extension
      mGraphFile.replace_extension(".graph.dot");
    }
  }
};

std::optional<Arguments> commandParser(int argc, char *argv[]) {
  CLI::App app{"Graph info"};

  std::string binFilename = "./data/extracted/wiki-Vote.txt.bin";
  bool infoFilename;
  std::string dotDirectory;
  bool verbose = false;
  bool progress = false;
  app.add_option("--input", binFilename, "Input bin file")->check(CLI::ExistingFile);
  app.add_flag("--output-info", infoFilename, "Output info file");
  app.add_option("--dot-directory", dotDirectory, "Output dot graph directory")->default_val("./data/actual");
  app.add_flag("-v,--verbose", verbose, "Verbose output");
  app.add_flag("-p,--progress", progress, "Show progress")->default_val(true);

  try{
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    if (app.exit(e) != 0) {
      throw e;
    }
    return std::nullopt;
  }

  if (verbose) {
    std::cout << "Filename: " << binFilename << std::endl;
  }

  return { { binFilename, infoFilename, dotDirectory, verbose, progress} };
}

std::string colorHash(uint32_t key) {
  auto color = colorMap.find(key);
  if (color == colorMap.end()) {
    std::cerr << "Error: color not found for key " << key << std::endl;
    std::abort();
  } 
  return color->second;
}

void saveGraph(Vf3Graph &graph, const std::string &dotFileName, std::string colorscheme) {
  std::ofstream dotFile(dotFileName.c_str());
  if (!dotFile.good()) {
    std::cerr << "Error opening file " << dotFileName << std::endl;
    std::abort();
  }
  dotFile << "digraph G {" << std::endl;
  //dotFile << "node [colorscheme=" << colorscheme << ",style=filled];" << std::endl;
  dotFile << "node [style=filled];" << std::endl;
  for (int i = 0; i < graph.nodes(); ++i) {
    // label with vertex id as subscript
    dotFile << i << " [fillcolor="<<colorHash(graph.label(i)+1)<<"];" << std::endl;
  }
  auto bothDirEdge = [&graph] (int src, int dst) {
    bool present = false;
    for (int k = 0; k < graph.outDegree(src); ++k) {
      if (graph.outEdge(src, k) == dst) {
        present = true;
      }
    }
    if (!present) {
      return false;
    }
    for (int k = 0; k < graph.outDegree(dst); ++k) {
      if (graph.outEdge(dst, k) == src) {
        return true;
      }
    }
    return false;
  };
  //dotFile << "edge [colorscheme=" << colorscheme << ",style=filled];" << std::endl;
  dotFile << "edge [style=filled];" << std::endl;
  for (int i = 0; i < graph.nodes(); ++i) {
    for (int j = 0; j < graph.outDegree(i); ++j) {
      auto [dst, edgelabel] = graph.GetOutEdgeWeight(i, j);
      if (i > dst && bothDirEdge(i, dst))
        continue;
      std::string arrow = " -> ";
      dotFile << i << arrow << dst << "[fillcolor="<< colorHash(edgelabel+1)
        << ",dir=" << (bothDirEdge(i, dst) ? "both" : "forward")
        << "];" << std::endl;
    }
  }
  dotFile << "}" << std::endl;
  dotFile.close();
}

int main(int argc, char * argv[]) {
  auto argsOpt = commandParser(argc, argv);

  if (!argsOpt) {
    return 255;
  }

  auto args = argsOpt.value();

  auto graph = args.mGraph;

  auto filename = args.mBinFile.filename();

  auto maxProgress = 40;
  auto printIterval = 1000;

  using namespace indicators;
  //show_console_cursor(false);

  ProgressBar bar {
    option::BarWidth{maxProgress},
    option::Start{"["},
    option::Fill{"■"},
    option::Lead{"■"},
    option::Remainder{"-"},
    option::End{" ]"},
    option::PostfixText{"checking maxDegree of " + filename.string()},
    option::ForegroundColor{Color::cyan},
    option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
    option::MaxProgress{graph.nodes()},
    option::ShowElapsedTime{true},
    option::ShowRemainingTime{true},
  };

  uint maxDegree = 0;
  for (int i = 0; i < graph.nodes(); ++i) {
    maxDegree = std::max(maxDegree, graph.outDegree(i));
    if (args.mProgress) {
      if (graph.nodes() > printIterval) {
        if (i % (graph.nodes() / printIterval) == 0 || i+1 == graph.nodes()) {
        bar.set_option(option::PostfixText{
          std::to_string(i+1) + "/" + std::to_string(graph.nodes()) + " [" + filename.string() + "]"
          });
        }
      } else {
        bar.set_option(option::PostfixText{
            std::to_string(i+1) + "/" + std::to_string(graph.nodes()) + " [" + filename.string() + "]"
            });
      }
      bar.tick();
    }
  }
  bar.mark_as_completed();

  uint edges = 0;
  std::set<uint> edgeLabels;
  for (int i = 0; i < graph.nodes(); ++i) {
    edges += graph.outDegree(i);
    for (int j = 0; j < graph.outDegree(i); ++j) {
      edgeLabels.insert(graph.GetOutEdgeWeight(i, j).second);
    }
  }
  assert(edgeLabels.size() == graph.edgeLabels());
  if (edgeLabels.size() != graph.edgeLabels()) {
    std::cerr << "Error: edgeLabels.size() ["<<edgeLabels.size()<<"] != graph.edgeLabels() ["<<graph.edgeLabels()<<"]" << std::endl;
    std::abort();
  }

  if (args.mVerbose) {
    std::cout << "Num of nodes : " << graph.nodes() << " labels : " << graph.labels() << std::endl;
    std::cout << "Num of edges : " << edges << " labels : " << graph.edgeLabels() << std::endl;
    std::cout << "Max degree : " << maxDegree << std::endl;
  }

  using json = nlohmann::json;

  if (args.mOutputFile != "") {
    std::ofstream out(args.mOutputFile.c_str());
    if (!out.good()) {
      std::cerr << "Error opening file " << args.mOutputFile << std::endl;
      std::abort();
    }
    json j;
    j["nodes"] = graph.nodes();
    j["labels"] = graph.labels();
    j["edges"] = edges;
    j["edgeLabels"] = graph.edgeLabels();
    j["maxDegree"] = maxDegree;
    out << j.dump(2) << std::endl;
    out.close();
  }

  if (args.mGraphFile != "") {
    saveGraph(graph, args.mGraphFile, args.mColorScheme);
  }

  show_console_cursor(true);

  return 0;

}
