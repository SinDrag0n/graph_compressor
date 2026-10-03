#pragma once

#include "common.h"

#include <algorithm>
#include <array>
#include <functional>
#include <queue>
#include <utility>

namespace NGraphCompressor {

    class THuffmanModel {
    public:
        explicit THuffmanModel(const std::array<uint32_t, 257>& frequencies) {
            struct TNode {
                uint64_t Frequency;
                int Left;
                int Right;
                int Symbol;
            };

            std::vector<TNode> nodes;
            using TQueueItem = std::pair<uint64_t, int>;
            std::priority_queue<TQueueItem, std::vector<TQueueItem>, std::greater<TQueueItem>> queue;
            for (unsigned i = 0; i < frequencies.size(); ++i) {
                if (frequencies[i] != 0) {
                    queue.push({frequencies[i], int(nodes.size())});
                    nodes.push_back({frequencies[i], -1, -1, int(i)});
                }
            }
            if (queue.empty()) {
                return;
            }

            while (queue.size() > 1) {
                const auto left = queue.top();
                queue.pop();
                const auto right = queue.top();
                queue.pop();
                const int parentIndex = nodes.size();
                nodes.push_back({left.first + right.first, left.second, right.second, -1});
                queue.push({left.first + right.first, parentIndex});
            }

            std::function<void(int, unsigned)> assignLengths = [&](int index, unsigned depth) {
                const auto& node = nodes[index];
                if (node.Symbol >= 0) {
                    Lengths[node.Symbol] = std::max(1u, depth);
                    return;
                }
                assignLengths(node.Left, depth + 1);
                assignLengths(node.Right, depth + 1);
            };
            assignLengths(queue.top().second, 0);
            BuildCanonicalCodes();
        }

        static THuffmanModel FromLengths(const std::array<unsigned, 257>& lengths) {
            THuffmanModel model;
            model.Lengths = lengths;
            model.BuildCanonicalCodes();
            return model;
        }

    private:
        THuffmanModel() = default;

        void BuildCanonicalCodes() {
            std::vector<std::pair<unsigned, unsigned>> symbolsByLength;
            for (unsigned i = 0; i < Lengths.size(); ++i) {
                if (Lengths[i] != 0) {
                    Require(Lengths[i] <= 63, "Huffman depth");
                    symbolsByLength.push_back({Lengths[i], i});
                }
            }
            std::sort(symbolsByLength.begin(), symbolsByLength.end());

            uint64_t nextCode = 0;
            unsigned previousLength = 0;
            for (auto [length, symbol] : symbolsByLength) {
                nextCode <<= length - previousLength;
                Require(nextCode < (uint64_t(1) << length), "invalid Huffman lengths");
                Codes[symbol] = nextCode++;
                previousLength = length;
            }
        }

    public:
        std::array<uint64_t, 257> Codes = {};
        std::array<unsigned, 257> Lengths = {};
    };

} // namespace NGraphCompressor
