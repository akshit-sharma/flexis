#pragma once

#include "def.h"
#include "instruction.hpp"
#include <cassert>
#include <cstdio>
#include <iostream>
#include <regex>
#include <string>
#include <fstream>

namespace generatorInstructions {

  using instruction::Instruction;
  using Instructions = std::vector<instruction::Instruction>;
  using instruction::EInstructionType;
  using instruction::VidType;
  using instruction::LabelType;
  using instruction::VertexSetIdxType;
  using instruction::SizeType;
  using instruction::InstructionIdxType;

  using instruction::callback;
  using instruction::init;
  using instruction::label;
  using instruction::neigh;
  using instruction::intersection;
  using instruction::difference;
  using instruction::loop;

  using PositionMap = std::map<std::string, VertexSetIdxType>;

  constexpr bool mTrace = false;

  template<class RestSet>
  void generateRecursiveDifferenceInstruction(
      PositionMap &positionMap,
      Instructions &instructions,
      const RestSet &restSet,
      const size_t recursiveStartingIndex,
      const RestSet &par,
      size_t depth_remaining)
  {

    assert(depth_remaining <= restSet.out.size());
    assert(restSet.label != -1);

    std::stringstream ssout;
    ssout << "y" << restSet.out[depth_remaining] << "l" << restSet.label;
    std::string secondvar = ssout.str();

    std::string firstvar = restSet.varname;

    if (depth_remaining == 0) {
      std::stringstream ss;
      ss << "y" << restSet.ins[0] << "l" << restSet.label;
      firstvar = ss.str();

      if (mTrace)
        std::cout << "firstvar: " << firstvar << std::endl;

      if (restSet.res_chain[restSet.depth] != -1) {
        std::stringstream ss;
        ss << "y" << restSet.res_chain[restSet.depth] << "l" << restSet.label;
        std::string secondvar = ss.str();

        instructions.emplace_back(difference(
            positionMap.at(restSet.varname),
            positionMap.at(firstvar),
            positionMap.at(secondvar)));
      }
      return;
      //assert(restSet.res_chain[restSet.depth] == -1);// not implemented
    }

    assert(positionMap.find(firstvar) != positionMap.end());
    assert(positionMap.find(secondvar) != positionMap.end());

    if (mTrace)
      std::cout << positionMap.at(restSet.varname) << " = "
        << positionMap.at(firstvar) << " - "
        << positionMap.at(secondvar) << std::endl;

    instructions.emplace(
        instructions.begin() + recursiveStartingIndex,
        difference(
          positionMap.at(restSet.varname),
          positionMap.at(firstvar),
          positionMap.at(secondvar)));

    if (depth_remaining > 0) {
      generateRecursiveDifferenceInstruction(
          positionMap,
          instructions,
          restSet,
          recursiveStartingIndex,
          par,
          depth_remaining - 1);
    }
  }

  template <class RestSet>
  void generateDifferenceInstruction(
      PositionMap &positionMap,
      instruction::Instructions &instructions,
      const RestSet &restSet,
      const size_t vsEffectiveIndex,
      const RestSet &par)
  {
    assert(restSet.out.size() != par.out.size());
    std::stringstream ss;
    ss << "y" << restSet.out[restSet.out.size() - 1] << "l" << restSet.label;
    std::string prevVarName = ss.str();
    assert(restSet.depth < restSet.res_chain.size());
    assert(positionMap.find(par.varname) != positionMap.end());
    assert(positionMap.find(prevVarName) != positionMap.end());
    positionMap[restSet.varname] = vsEffectiveIndex;
    if (mTrace)
      std::cout << "difference: BS[" << vsEffectiveIndex << "] = "
        << positionMap.at(par.varname) << " - " << positionMap.at(prevVarName) << std::endl;
    instructions.emplace_back(difference(
          VertexSetIdxType(vsEffectiveIndex),
          positionMap.at(par.varname),
          positionMap.at(prevVarName)));
  }

  template <class RestSet>
  void generateIntersectonInstruction(
      PositionMap &positionMap,
      Instructions &instructions,
      const RestSet &restSet,
      const size_t vsEffectiveIndex,
      const RestSet &par)
  {
    assert(restSet.ins.size() != par.ins.size());
    std::stringstream ss;
    ss << "y" << restSet.ins[restSet.ins.size() - 1] << "l" << restSet.label;
    std::string prevVarName = ss.str();
    assert(restSet.depth < restSet.res_chain.size());
    assert(positionMap.find(par.varname) != positionMap.end());
    assert(positionMap.find(prevVarName) != positionMap.end());
    positionMap[restSet.varname] = vsEffectiveIndex;
    if (mTrace)
      std::cout << "intersection: BS[" << vsEffectiveIndex << "] = "
        << positionMap.at(par.varname) << " - " << positionMap.at(prevVarName) << std::endl;
    instructions.emplace_back(intersection(
          VertexSetIdxType(vsEffectiveIndex),
          positionMap.at(par.varname),
          positionMap.at(prevVarName)));
  }

  template<class RestSet, class MultiRestPlan>
  size_t generatePlanInstructions(
      Instructions &instructions,
      size_t pForLoop,
      PositionMap &positionMap,
      size_t pNewVSGenerated,
      const MultiRestPlan *mrp,
      size_t nextVSId)
  {
    //	assert(mrp.size() == 1 || (mrp.size() > 1 && mrp.children.empty()));
    for (const auto &restSet : mrp->atlev) {
      assert(positionMap.find(restSet.varname) == positionMap.end());
      assert(instructions.size() - pForLoop >= 0);
      size_t vsEffectiveIndex = nextVSId++;
      if (restSet.ins.empty() && restSet.out.empty()) {
        positionMap[restSet.varname] = vsEffectiveIndex;
        if (mTrace)
          std::cout << "label: BS[" << vsEffectiveIndex << "] = g.L("
            << restSet.label << ")" << std::endl;
        instructions.emplace_back(label(VertexSetIdxType(vsEffectiveIndex), LabelType(restSet.label), SizeType(0)));
        continue;
      }
      if (restSet.ins.size() == 1 && restSet.out.empty()) {
        if (restSet.res_chain[restSet.depth] != -1) {
          std::stringstream ss;
          ss << "y" << restSet.ins[0] << "l" << restSet.label;
          std::string prevVarName = ss.str();
          assert(restSet.varname != prevVarName);
          assert(positionMap.find(prevVarName) != positionMap.end());
          if (mTrace)
            std::cout << "remapping " << restSet.varname << " to " << prevVarName << "(" << positionMap.at(prevVarName) << ")" << std::endl;
          positionMap[restSet.varname] = positionMap.at(prevVarName);
          nextVSId--;
          continue;
        }
        positionMap[restSet.varname] = vsEffectiveIndex;
        if (mTrace)
          std::cout << "neigh: BS[" << vsEffectiveIndex << "] = g.N(v"
            << restSet.ins[0] << ", l" << restSet.label << ")" << std::endl;
        instructions.emplace_back(neigh(
              VertexSetIdxType(vsEffectiveIndex),
              VidType(restSet.depth),
              LabelType(restSet.label)));
        continue;
      }
      RestSet par = restSet.parent();
      if (par.ins.empty()) {
        assert(par.ins.empty());
        if (par.out.empty()) { assert(false); }
        assert(!restSet.out.empty());
        assert(positionMap.find(restSet.varname) == positionMap.end());
        positionMap[restSet.varname] = vsEffectiveIndex;
        generateRecursiveDifferenceInstruction(
            positionMap,
            instructions,
            restSet,
            instructions.size(),
            par,
            restSet.out.size() - 1);
        assert(positionMap.find(restSet.varname) != positionMap.end());
        assert(positionMap[restSet.varname] == vsEffectiveIndex);
        pNewVSGenerated++;
      } else if (restSet.out.size() != par.out.size()) {
        generateDifferenceInstruction(
            positionMap, instructions, restSet, vsEffectiveIndex, par);
        pNewVSGenerated++;
      } else if (restSet.ins.size() != par.ins.size()) {
        generateIntersectonInstruction(
            positionMap, instructions, restSet, vsEffectiveIndex, par);
        pNewVSGenerated++;
      }
    }

    const auto &restset = mrp->children.empty() ? mrp->counters.begin()->first
      : mrp->children.begin()->first;
    assert(positionMap.find(restset.varname) != positionMap.end());
    size_t num_instructions = instructions.size();
    instructions.emplace_back(loop(
          positionMap[restset.varname],
          InstructionIdxType(instructions.size() + 1)));

    if (mrp->children.empty()) {
      assert(!mrp->counters.empty());
      assert(instructions[0].type() == EInstructionType::INIT);
      instructions[0] = init(
          VidType(instructions[0].first().vid),
          SizeType(instructions.size() - pForLoop - 1)); // 1 is for init instruction
      instructions[1] = label(
          instructions[1].dst(),
          LabelType(instructions[1].first().label),
          SizeType(pNewVSGenerated));
      instructions.emplace_back(callback());
      return num_instructions;
    }

    pForLoop++;

    MultiRestPlan *childrenMrp = mrp->children.begin()->second;

    instructions[num_instructions] = loop(
        positionMap[restset.varname],
        InstructionIdxType(generatePlanInstructions<RestSet>(
            instructions, pForLoop, positionMap,
            pNewVSGenerated, childrenMrp, nextVSId)));

    return num_instructions;
  }

  template <typename RestSet, typename MultiRestPlan>
  Instructions generateInstructions(const MultiRestPlan &mrp)
  {
    PositionMap positionMap;
    Instructions instructions;
    assert(mrp.size() == 1);
    instructions.emplace_back(init(VidType(mrp.getPlan(0).vertices()), SizeType(0)));
    generatePlanInstructions<RestSet>(instructions, 1, positionMap, 0, &mrp, 0);

    return instructions;
  }

  std::vector<std::string> SplitString(std::string str, std::string delimeter)
  {

    auto ltrim = [](std::string &s) {
      s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
            }));
    };

    // trim from end (in place)
    auto rtrim = [](std::string &s) {
      s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
            }).base(), s.end());
    };

    // trim from both ends (in place)
    auto trim = [&ltrim, &rtrim](std::string s) {
      ltrim(s);
      rtrim(s);
      return s;
    };

    std::vector<std::string> splittedStrings = {};
    size_t pos = 0;

    while ((pos = str.find(delimeter)) != std::string::npos){
      std::string token = str.substr(0, pos);
      if (token.length() > 0)  splittedStrings.push_back(trim(token));
      str.erase(0, pos + delimeter.length());
    }

    if (str.length() > 0)  splittedStrings.push_back(trim(str));
    return splittedStrings;
  }

  Instructions generateInstructions(const std::string &planfile)
  {
    Instructions instructions;
    std::ifstream file(planfile);
    assert(file.is_open());
    std::string line;
    int forNum = 0;
    int variableNum = 0;
    int variableSetNum = 0;
    std::map<std::string, int> variableMap;
    std::map<std::string, int> variableLabel;
    std::map<std::string, int> forLoopVar;
    int lastLoopIdx = -1;
    std::regex init_regex("std::vector<vidType>\\([0-9]+,[0-9]+\\)");
    std::regex var_regex("auto [[:alnum:]]+");
    int var_prefix = std::string("auto ").size();
    std::regex labelLine_regex("auto [[:alnum:]]+ = g.L\\([[:digit:]]+\\)");
    std::regex for_regex("for\\(vidType [[:alnum:]]+:[[:alnum:]|_]+");
    int for_prefix = std::string("for(vidType ").size();
    std::regex neighbor_regex("auto [[:alnum:]]+ = g.N\\([[:alnum:]]+,[[:blank:]]*[[:digit:]]+");
    std::regex neighbor_regex_2("g.N\\([[:alnum:]]+,[[:blank:]]*[[:digit:]]+");
    int neighbor_prefix_2 = std::string("g.N(").size();
    std::regex var_alias_regex("auto &[[:alnum:]]+ = [[:alnum:]]+");
    int var_alias_prefix = std::string("auto &").size();
    std::regex is_regex("auto [[:alnum:]|_]+ = SetOp::intersection_set MATCHING ");
    std::regex ds_regex("auto [[:alnum:]|_]+ = SetOp::difference_set MATCHING ");
    std::regex ds_rec_regex("auto [[:alnum:]|_]+ = MATCHING_DIFF MATCHING");
    std::regex set_regex("\\([[:alnum:]|_]+,[[:blank:]]+[[:alnum:]]+");
    while(std::getline(file, line)) {
      // std::cout  << std::endl << line << " ";
      std::smatch init_match;
      std::smatch var_match;
      std::smatch label_match;
      std::smatch for_match;
      std::smatch neighbor_match;
      std::smatch var_alias_match;
      std::smatch set_match;
      if (std::regex_search(line, std::regex("mMatchingVec.at\\(tid\\)"))) continue;
      if (std::regex_search(line, std::regex("callbackVec\\[tid\\]"))) { instructions.emplace_back(callback()); continue; }
      if (std::regex_search(line, set_match, ds_rec_regex)) {
        std::string match = set_match.str();
        std::string variable = match.substr(var_prefix, match.find("=")-1-var_prefix);
        assert(variableMap.find(variable) == variableMap.end());
        if (std::regex_search(line, set_match, set_regex)) {
          std::vector<std::string> param_vecs = SplitString(set_match.str().substr(1), ",");
          assert(param_vecs.size() == 2);
          assert(variableMap.find(param_vecs[0]) != variableMap.end());
          assert(variableMap.find(param_vecs[1]) != variableMap.end());
          variableMap[variable] = variableMap[param_vecs[0]];
          // std::cout << "DS REC " << variable << " " << variableMap[variable];
          instructions.emplace_back(difference(VertexSetIdxType(variableMap[variable]),
                VertexSetIdxType(variableMap[param_vecs[0]]),
                VertexSetIdxType(variableMap[param_vecs[1]])));
        }
      } else if (std::regex_search(line, set_match, is_regex)) {
        std::string match = set_match.str();
        std::string variable = match.substr(var_prefix, match.find("=")-1-var_prefix);
        assert(variableMap.find(variable) == variableMap.end());
        variableMap[variable] = variableNum++;
        // std::cout << "IS " << variable << " " << variableMap[variable];
        if (std::regex_search(line, set_match, set_regex)) {
          std::vector<std::string> param_vecs = SplitString(set_match.str().substr(1), ",");
          assert(param_vecs.size() == 2);
          assert(variableMap.find(param_vecs[0]) != variableMap.end());
          assert(variableMap.find(param_vecs[1]) != variableMap.end());
          instructions.emplace_back(intersection(
                VertexSetIdxType(variableMap[variable]),
                VertexSetIdxType(variableMap[param_vecs[0]]),
                VertexSetIdxType(variableMap[param_vecs[1]])));
          variableSetNum++;
        }
      } else if (std::regex_search(line, set_match, ds_regex)) {
        std::string match = set_match.str();
        std::string variable = match.substr(var_prefix, match.find("=")-1-var_prefix);
        assert(variableMap.find(variable) == variableMap.end());
        variableMap[variable] = variableNum++;
        // std::cout << "DS " << variable << " " << variableMap[variable];
        if (std::regex_search(line, set_match, set_regex)) {
          std::vector<std::string> param_vecs = SplitString(set_match.str().substr(1), ",");
          assert(param_vecs.size() == 2);
          assert(variableMap.find(param_vecs[0]) != variableMap.end());
          assert(variableMap.find(param_vecs[1]) != variableMap.end());
          instructions.emplace_back(intersection(
                VertexSetIdxType(variableMap[variable]),
                VertexSetIdxType(variableMap[param_vecs[0]]),
                VertexSetIdxType(variableMap[param_vecs[1]])));
          variableSetNum++;
        }
      } else if (std::regex_search(line, var_alias_match, var_alias_regex)) {
        // std::cout << "alias ";
        std::string var_alias_match_string = var_alias_match.str();
        std::string alias_variable = var_alias_match_string.substr(var_alias_prefix, var_alias_match_string.find("=")-1-var_alias_prefix);
        std::string actual_variable = var_alias_match_string.substr(var_alias_match_string.find("=")+2);
        assert(variableMap.find(actual_variable) != variableMap.end());
        variableMap[alias_variable] = variableMap[actual_variable];
        // std::cout << "variable " << alias_variable << " " << actual_variable << " " << variableMap[alias_variable] << " " << variableMap[actual_variable];
      } else if (std::regex_search(line, neighbor_match, neighbor_regex)) {
        // std::cout << "neigh ";
        std::string variable = neighbor_match.str();
        variable = variable.substr(var_prefix, variable.find(" ", var_prefix+1)-var_prefix);
        variableMap[variable] = variableNum++;
        // std::cout << variable << " " << variableMap[variable];
        std::string neighbor_string = neighbor_match.str();
        std::smatch neighbor_2_match;
        if (std::regex_search(line, neighbor_2_match, neighbor_regex_2)) {
          std::string neigh_parameters = neighbor_2_match.str().substr(neighbor_prefix_2);
          std::vector<std::string> neigh_parameters_vec = SplitString(neigh_parameters, ",");
          assert(neigh_parameters_vec.size() == 2);
          variableLabel[variable] = atoi(neigh_parameters_vec[1].c_str());
          instructions.emplace_back(neigh(
                VertexSetIdxType(variableMap[variable]),
                VidType(forLoopVar[neigh_parameters_vec[0]]),
                LabelType(variableLabel[variable])));
        } else {
          assert(false);
        }
      } else if (std::regex_search(line, for_match, for_regex)) {
        // std::cout << "for ";
        std::string for_match_string = for_match.str().substr(for_prefix);
        if (lastLoopIdx != -1) instructions[lastLoopIdx] = loop(
            VertexSetIdxType(instructions[lastLoopIdx].first().vertexSetIdx),
            InstructionIdxType(instructions.size()));
        std::vector<std::string> loop_vec = SplitString(for_match_string, ":");
        assert(loop_vec.size() == 2);
        forLoopVar[loop_vec[0]] = forNum++;
        // std::cout << loop_vec[0] << " : " << loop_vec[1] << "(" << variableMap[loop_vec[1]] << ")" <<  " with " << forNum;
        instructions.emplace_back(loop(
              VertexSetIdxType(variableMap[loop_vec[1]]),
              InstructionIdxType(0)));
        lastLoopIdx = instructions.size()-1;
      } else if (std::regex_search(line, var_match, var_regex) && std::regex_search(line, label_match, labelLine_regex)) {
        // std::cout << "label ";
        std::regex label_regex("\\([0-9]+");
        std::string parameters_string = label_match.str();
        std::smatch parameters_match;
        if (std::regex_search(parameters_string, parameters_match, label_regex)) {
          std::string variable = var_match.str().substr(var_prefix);
          variableMap[variable] = variableNum++;
          // std::cout << "variable " << variable << " " << variableMap[variable];
          int label_var = atoi(parameters_match.str().substr(1).c_str());
          variableLabel[variable] = label_var;
          instructions.emplace_back(label(VertexSetIdxType(0), LabelType(label_var), SizeType(0)));
          assert(variableMap[variable] == 0);
        }
      } else if (std::regex_search(line, init_match, init_regex)) {
        // std::cout << "init ";
        assert(init_match.size() == 1);
        std::regex parameters_regex("[0-9]+,[0-9]+");
        std::smatch parameters_match;
        std::string parameters_string = init_match.str();
        if (std::regex_search(parameters_string, parameters_match, parameters_regex)) {
          std::vector<std::string> parameters_vec = SplitString(parameters_match.str(),",");
          assert(parameters_vec.size() == 2);
          instructions.emplace_back(init(
                VidType(atoi(parameters_vec[0].c_str())),
                SizeType(2)));
        }
      }
    }
    file.close();
    for (auto [name, idx] : variableMap) {
      assert(name.find(" ") == std::string::npos);
    }
    instructions[0] = init(VidType(instructions[0].first().vid), SizeType(variableNum));
    instructions[1] = label(VertexSetIdxType(0), LabelType(instructions[1].first().label), SizeType(variableSetNum));
    instructions[lastLoopIdx] = loop(VertexSetIdxType(instructions[lastLoopIdx].first().vertexSetIdx), InstructionIdxType(instructions.size()-1));
    return instructions;
  }

} // namespace generator

