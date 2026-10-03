#include "graph.h"
#include "byte_io.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include <cstdio>

namespace NGraphCompressor {
    namespace {

        // uint64_t sentinel cannot collide with a valid uint32_t vertex ID.
        class TIdMap {
        public:
            explicit TIdMap(const std::vector<uint32_t>& originalIds) {
                size_t capacity = 1;
                while (capacity < originalIds.size() * 2 + 1) {
                    capacity *= 2;
                }
                Mask_ = capacity - 1;
                Keys_.assign(capacity, UINT64_MAX);
                Values_.resize(capacity);
                for (size_t i = 0; i < originalIds.size(); ++i) {
                    size_t slot = Hash(originalIds[i]) & Mask_;
                    while (Keys_[slot] != UINT64_MAX) {
                        slot = (slot + 1) & Mask_;
                    }
                    Keys_[slot] = originalIds[i];
                    Values_[slot] = i;
                }
            }

            uint32_t FindDenseId(uint32_t originalId) const {
                size_t slot = Hash(originalId) & Mask_;
                while (Keys_[slot] != originalId) {
                    Require(Keys_[slot] != UINT64_MAX, "missing ID");
                    slot = (slot + 1) & Mask_;
                }
                return Values_[slot];
            }

        private:
            static uint64_t Hash(uint64_t value) {
                value ^= value >> 30;
                value *= 0xbf58476d1ce4e5b9ULL;
                value ^= value >> 27;
                value *= 0x94d049bb133111ebULL;
                return value ^ (value >> 31);
            }

        private:
            std::vector<uint64_t> Keys_;
            std::vector<uint32_t> Values_;
            size_t Mask_ = 0;
        };

    } // anonymous namespace

    TGraph ReadGraphTsv(const std::string& path) {
        double phase = GetTimeSeconds();
        TBytes inputBytes = ReadFile(path);
        size_t position = 0;
        TGraph graph;
        graph.Edges.reserve(inputBytes.size() / 24);
        std::vector<uint32_t> originalIds;
        auto readNumber = [&]() {
            Require(position < inputBytes.size() && inputBytes[position] >= '0' && inputBytes[position] <= '9', "invalid TSV integer");
            uint64_t value = 0;
            while (position < inputBytes.size() && inputBytes[position] >= '0' && inputBytes[position] <= '9') {
                value = value * 10 + (inputBytes[position++] - '0');
                Require(value <= UINT32_MAX, "TSV overflow");
            }
            return uint32_t(value);
        };
        while (position < inputBytes.size()) {
            uint32_t sourceId = readNumber();
            Require(position < inputBytes.size() && inputBytes[position++] == '\t', "expected tab");
            uint32_t targetId = readNumber();
            Require(position < inputBytes.size() && inputBytes[position++] == '\t', "expected tab");
            uint32_t weight = readNumber();
            Require(weight <= 255, "weight overflow");
            if (position < inputBytes.size() && inputBytes[position] == '\r') {
                ++position;
            }
            if (position < inputBytes.size()) {
                Require(inputBytes[position++] == '\n', "expected newline");
            }
            if (sourceId > targetId) {
                std::swap(sourceId, targetId);
            }
            graph.Edges.push_back({sourceId, targetId, uint8_t(weight)});
            originalIds.push_back(sourceId);
            originalIds.push_back(targetId);
        }
        std::cerr << "TSV_read_parse_seconds=" << GetTimeSeconds() - phase << "\n";
        phase = GetTimeSeconds();
        std::sort(originalIds.begin(), originalIds.end());
        originalIds.erase(std::unique(originalIds.begin(), originalIds.end()), originalIds.end());
        Require(originalIds.size() <= UINT32_MAX, "too many vertices");
        graph.OriginalIds = std::move(originalIds);
        std::cerr << "ID_sort_unique_seconds=" << GetTimeSeconds() - phase << "\n";
        phase = GetTimeSeconds();
        TIdMap idMap(graph.OriginalIds);
        for (auto& edge : graph.Edges) {
            edge.Source = idMap.FindDenseId(edge.Source);
            edge.Target = idMap.FindDenseId(edge.Target);
        }
        std::cerr << "dense_mapping_seconds=" << GetTimeSeconds() - phase << "\n";
        phase = GetTimeSeconds();
        std::sort(graph.Edges.begin(), graph.Edges.end(), CompareEdges);
        graph.ForwardDegrees.assign(graph.OriginalIds.size(), 0);
        for (size_t i = 0; i < graph.Edges.size(); i++) {
            auto& edge = graph.Edges[i];
            if (i) {
                Require(edge.Source != graph.Edges[i - 1].Source || edge.Target != graph.Edges[i - 1].Target, "parallel edges");
            }
            ++graph.ForwardDegrees[edge.Source];
        }
        std::cerr << "edge_sort_seconds=" << GetTimeSeconds() - phase << "\n";
        return graph;
    }

    void WriteGraphTsv(const std::string& path, const TGraph& graph) {
        FILE* output = fopen(path.c_str(), "wb");
        Require(output, "cannot open TSV output");
        std::vector<char> outputBuffer;
        outputBuffer.reserve(1 << 20);
        auto flushBuffer = [&]() {
            if (fwrite(outputBuffer.data(), 1, outputBuffer.size(), output) != outputBuffer.size()) {
                fclose(output);
                throw std::runtime_error("TSV write failed");
            }
            outputBuffer.clear();
        };
        auto appendNumber = [&](uint32_t value) {
            char digits[10];
            unsigned digitCount = 0;
            do {
                digits[digitCount++] = char('0' + value % 10);
                value /= 10;
            } while (value);
            while (digitCount) {
                outputBuffer.push_back(digits[--digitCount]);
            }
        };
        for (auto& edge : graph.Edges) {
            appendNumber(graph.OriginalIds[edge.Source]);
            outputBuffer.push_back('\t');
            appendNumber(graph.OriginalIds[edge.Target]);
            outputBuffer.push_back('\t');
            appendNumber(edge.Weight);
            outputBuffer.push_back('\n');
            if (outputBuffer.size() > (1 << 20) - 40) {
                flushBuffer();
            }
        }
        flushBuffer();
        Require(!fclose(output), "TSV close failed");
    }

} // namespace NGraphCompressor
