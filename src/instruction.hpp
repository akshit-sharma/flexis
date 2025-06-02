#pragma once

#include <cassert>
#include <iomanip>
#include <ostream>
#include <sys/types.h>
#include <tuple>
#include <vector>

namespace instruction {

  class Instruction;
  using Instructions = std::vector<Instruction>;

  template<typename T, typename Parameter>
  class StrongType {
    public:
      using value_type = T;
      explicit StrongType() : mValue(0) {}
      explicit StrongType(T value) : mValue(value) {}
      value_type &operator()() { return mValue; }
      template <typename U>
      friend bool operator==(const StrongType &lhs, const U &rhs) {
        return lhs.mValue == rhs;
      }
      // typename U, where U is equal to size_t type
      template <typename U, typename = typename std::enable_if<std::is_same<U, size_t>::value>::type>
      StrongType &operator=(U value);
      friend std::ostream& operator<<(std::ostream& os, const StrongType& strongType) {
        os << strongType.mValue;
        return os;
      }
    private:
      value_type mValue;
  };


  enum class EInstructionType {
    CALLBACK,
    INIT,
    LABEL,
    NEIGHBOR,
    INTERSECTION,
    DIFFERENCE,
    FORLOOP,
  };

  using VidType = StrongType<uint, struct VidTypeTag>;
  using LabelType = StrongType<uint, struct LabelTypeTag>;
  using VertexSetIdxType = StrongType<uint, struct VertexSetIdxTypeTag>;
  using SizeType = StrongType<uint, struct SizeTypeTag>;
  using InstructionIdxType = StrongType<uint, struct InstructionIdxTypeTag>;

  template<> template <>
  VertexSetIdxType &VertexSetIdxType::operator=<size_t>(size_t value) { // NOLINT(misc-definitions-in-headers)
    mValue = value;
    return *this;
  }

  class Instruction {
    public:
      union FirstArg {
        LabelType::value_type label;
        VidType::value_type vid;
        VertexSetIdxType::value_type vertexSetIdx;
        explicit FirstArg() : vertexSetIdx{0} {}
        explicit FirstArg(LabelType label) : label(label()) {}
        explicit FirstArg(VidType vid) : vid(vid()) {}
        explicit FirstArg(VertexSetIdxType vertexSetIdx) : vertexSetIdx(vertexSetIdx()) {}
      };
      union SecondArg {
        LabelType::value_type label;
        VertexSetIdxType::value_type vertexSetIdx;
        SizeType::value_type size;
        InstructionIdxType::value_type instructionIdx;
        explicit SecondArg() : vertexSetIdx{0} {}
        explicit SecondArg(LabelType label) : label(label()) {}
        explicit SecondArg(VertexSetIdxType vertexSetIdx) : vertexSetIdx(vertexSetIdx()) {}
        explicit SecondArg(SizeType size) : size(size()) {}
        explicit SecondArg(InstructionIdxType instructionIdx) : instructionIdx(instructionIdx()) {}
      };
      union ReturnType {
        VertexSetIdxType vertexSetIdx;
        SizeType::value_type size;
        explicit ReturnType() : vertexSetIdx{0} {}
        explicit ReturnType(VertexSetIdxType vertexSetIdx) : vertexSetIdx(vertexSetIdx()) {}
        explicit ReturnType(SizeType size) : size(size()) {}
      };

      Instruction() : mType(EInstructionType::CALLBACK) {}

      EInstructionType type() const { return mType; }
      VertexSetIdxType dst() const {
        assert(mType != EInstructionType::CALLBACK);
        assert(mType != EInstructionType::INIT);
        return mDstIndex.vertexSetIdx;
      }

      FirstArg first() const { return mFirstArg; }
      SecondArg second() const { return mSecondArg; }

      std::tuple<EInstructionType, VertexSetIdxType, FirstArg, SecondArg> operator()() const {
        return std::make_tuple(mType, mDstIndex.vertexSetIdx, mFirstArg, mSecondArg);
      }

    private:
      EInstructionType mType;
      ReturnType mDstIndex;
      FirstArg mFirstArg;
      SecondArg mSecondArg;

      Instruction(VidType firstArg, SizeType bufferSize)
        : mType(EInstructionType::INIT), mDstIndex(bufferSize), mFirstArg(firstArg), mSecondArg(InstructionIdxType(2)) {}
      Instruction(VertexSetIdxType dstIdx, LabelType label, SizeType extraNeeded)
        : mType(EInstructionType::LABEL), mDstIndex(dstIdx), mFirstArg(label), mSecondArg(extraNeeded) {}
      Instruction(VertexSetIdxType dstIdx, VidType vid, LabelType label)
        : mType(EInstructionType::NEIGHBOR), mDstIndex(dstIdx), mFirstArg(vid), mSecondArg(label) {}
      Instruction(VertexSetIdxType vsIdx, InstructionIdxType instructionIdx)
        : mType(EInstructionType::FORLOOP), mDstIndex(VertexSetIdxType(0)), mFirstArg(vsIdx), mSecondArg(instructionIdx) {}
      Instruction(EInstructionType type, VertexSetIdxType dstIdx, VertexSetIdxType firstArg, VertexSetIdxType secondArg)
        : mType(type), mDstIndex(dstIdx), mFirstArg(firstArg), mSecondArg(secondArg) {}

    public:
      friend Instruction init(VidType vertices, SizeType bufferSize);
      friend Instruction label(VertexSetIdxType dstIdx, LabelType label, SizeType extraNeeded);
      friend Instruction neigh(VertexSetIdxType dstIdx, VidType vid, LabelType label);
      friend Instruction loop(VertexSetIdxType vsIdx, InstructionIdxType instructionIdx);
      friend Instruction intersection(VertexSetIdxType dstIdx, VertexSetIdxType srcIdx1, VertexSetIdxType srcIdx2);
      friend Instruction difference(VertexSetIdxType dstIdx, VertexSetIdxType srcIdx1, VertexSetIdxType srcIdx2);
      friend Instruction callback();

      friend std::ostream& operator<<(std::ostream& os, const Instruction& instruction);

      friend uint bufferExpected(Instructions instructions);
  };

  Instruction init(VidType vertices, SizeType bufferSize) { // NOLINT(misc-definitions-in-headers)
    return Instruction(vertices, bufferSize);
  }
  Instruction label(VertexSetIdxType dstIdx, LabelType label, SizeType extraNeeded) { // NOLINT(misc-definitions-in-headers)
    return Instruction(dstIdx, label, extraNeeded);
  }
  Instruction neigh(VertexSetIdxType dstIdx, VidType vid, LabelType label) { // NOLINT(misc-definitions-in-headers)
    return Instruction(dstIdx, vid, label);
  }
  Instruction loop(VertexSetIdxType vsIdx, InstructionIdxType instructionIdx) { // NOLINT(misc-definitions-in-headers)
    return Instruction(vsIdx, instructionIdx);
  }
  Instruction intersection(VertexSetIdxType dstIdx, VertexSetIdxType srcIdx1, VertexSetIdxType srcIdx2) { // NOLINT(misc-definitions-in-headers)
    return Instruction(EInstructionType::INTERSECTION, dstIdx, srcIdx1, srcIdx2);
  }
  Instruction difference(VertexSetIdxType dstIdx, VertexSetIdxType srcIdx1, VertexSetIdxType srcIdx2) { // NOLINT(misc-definitions-in-headers)
    return Instruction(EInstructionType::DIFFERENCE, dstIdx, srcIdx1, srcIdx2);
  }
  Instruction callback() { // NOLINT(misc-definitions-in-headers)
    return Instruction();
  }

  std::ostream& operator<<(std::ostream& os, const Instruction& instruction) { // NOLINT(misc-definitions-in-headers)
    using EInstructionType::CALLBACK;
    using EInstructionType::LABEL;
    using EInstructionType::NEIGHBOR;
    using EInstructionType::INTERSECTION;
    using EInstructionType::DIFFERENCE;
    using EInstructionType::FORLOOP;
    using EInstructionType::INIT;

    switch (instruction.type()) {
      case CALLBACK:
        return os << "CALLBACK";
      case DIFFERENCE:
        return os << "B[" << instruction.dst() << "] = B[" << instruction.first().vertexSetIdx << "] - B[" << instruction.second().vertexSetIdx << "]";
      case INTERSECTION:
        return os << "B[" << instruction.dst() << "] = B[" << instruction.first().vertexSetIdx << "] & B[" << instruction.second().vertexSetIdx << "]";
      case FORLOOP:
        return os << "For v" << instruction.first().vertexSetIdx << " IN B[" << instruction.dst() << "] -> next loop at #" << instruction.second().instructionIdx;
      case NEIGHBOR:
        return os << "B[" << instruction.dst() << "] = Neigh(v" << instruction.first().vid << ", l" << instruction.second().label << ")";
      case LABEL:
        return os << "B[" << instruction.dst() << "] = Label(l" << instruction.first().label << ") with extraNeeded " << instruction.second().size;
      case INIT:
        return os << "numVertices[" << instruction.first().vid << "] BS[" << instruction.mDstIndex.size << "] next loop at #[" << instruction.second().instructionIdx << "]";
    }

    return os;
  }

  std::ostream& operator<<(std::ostream& os, const Instructions& instructions) { // NOLINT(misc-definitions-in-headers)
    assert(instructions.size() > 2);
    assert(instructions[0].type() == EInstructionType::INIT);
    assert(instructions[1].type() == EInstructionType::LABEL);
    assert(instructions[2].type() == EInstructionType::FORLOOP);
    os << "Meta information (" << instructions.size() << "):\n";
    uint spaceIntend = 0;
    auto printTab = [&os](uint spaceIntend) {
      for (uint i = 0; i < spaceIntend; ++i) {
        os << "  ";
      }
    };
    uint lineNum = 0;
    for (const auto& instruction : instructions) {
      // format lineNum with two digits
      os << std::setw(2) << std::setfill('0') << lineNum << ": ";
      printTab(spaceIntend);
      os << instruction << std::endl;
      if (instruction.type() == EInstructionType::FORLOOP) {
        spaceIntend++;
      }
      lineNum++;
    }
    return os;
  }

  std::ostream &operator<<(std::ostream &os, const EInstructionType &type) { // NOLINT(misc-definitions-in-headers)
    switch (type) {
      case EInstructionType::CALLBACK:
        return os << "CALLBACK";
      case EInstructionType::DIFFERENCE:
        return os << "DIFFERENCE";
      case EInstructionType::FORLOOP:
        return os << "FORLOOP";
      case EInstructionType::INIT:
        return os << "INIT";
      case EInstructionType::INTERSECTION:
        return os << "INTERSECTION";
      case EInstructionType::LABEL:
        return os << "LABEL";
      case EInstructionType::NEIGHBOR:
        return os << "NEIGHBOR";
      default:
        return os << "UNKNOWN";
        assert(false);
    }
    return os;
  }

  uint bufferExpected(Instructions instructions) { // NOLINT(misc-definitions-in-headers)
    assert(instructions.size() > 2);
    assert(instructions[0].type() == EInstructionType::INIT);
    return instructions[0].mDstIndex.size;
  }

  uint extraNeeded(Instructions instructions) { // NOLINT(misc-definitions-in-headers)
    assert(instructions.size() > 2);
    assert(instructions[0].type() == EInstructionType::INIT);
    assert(instructions[1].type() == EInstructionType::LABEL);
    return instructions[1].second().size;
  }

  uint numVertices(Instructions instructions) { // NOLINT(misc-definitions-in-headers)
    assert(instructions.size() > 2);
    assert(instructions[0].type() == EInstructionType::INIT);
    return instructions[0].first().vid;
  }

};
