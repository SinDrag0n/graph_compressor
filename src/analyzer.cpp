#include "analyzer.h"
#include "bitstream.h"
#include "codecs.h"
#include "format.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <utility>
#include <cmath>

namespace NGraphCompressor {

    void PrintStatistics(const char* name, std::vector<uint32_t> values, const std::string& histogramDirectory) {
        std::map<uint32_t, uint64_t> frequencies;
        for (auto value : values) {
            ++frequencies[value];
        }
        double entropy = 0;
        for (auto [value, frequency] : frequencies) {
            (void) value;
            double probability = double(frequency) / std::max(size_t(1), values.size());
            entropy -= probability * std::log2(probability);
        }
        if (!histogramDirectory.empty()) {
            std::ofstream out(histogramDirectory + "/" + name + ".csv");
            Require(bool(out), "histogram output");
            out << "value,count\n";
            for (auto [value, frequency] : frequencies) {
                out << value << "," << frequency << "\n";
            }
        }
        std::sort(values.begin(), values.end());
        std::cerr << name << " entropy=" << entropy << " max=" << (values.empty() ? 0 : values.back());
        for (double quantile : {.5,.9,.95,.99,.999 }) {
            std::cerr << " q" << quantile << "=" << (values.empty() ? 0 : values[std::min(values.size() -1, size_t(quantile * values.size()))]);
        }
        std::cerr << '\n';
    }

    TGraph ReorderGraph(const TGraph& graph, const std::vector<uint32_t>& order) {
        TGraph reorderedGraph;
        size_t vertexCount = graph.OriginalIds.size();
        std::vector<uint32_t> inverseOrder(vertexCount);
        reorderedGraph.OriginalIds.resize(vertexCount);
        reorderedGraph.ForwardDegrees.assign(vertexCount, 0);
        for (size_t i = 0; i < vertexCount; i++) {
            inverseOrder[order[i]] = i;
            reorderedGraph.OriginalIds[i] = graph.OriginalIds[order[i]];
        }
        reorderedGraph.Edges = graph.Edges;
        for (auto& edge : reorderedGraph.Edges) {
            edge.Source = inverseOrder[edge.Source];
            edge.Target = inverseOrder[edge.Target];
            if (edge.Source > edge.Target) {
                std::swap(edge.Source, edge.Target);
            }
            ++reorderedGraph.ForwardDegrees[edge.Source];
        }
        std::sort(reorderedGraph.Edges.begin(), reorderedGraph.Edges.end(), CompareEdges);
        return reorderedGraph;
    }

    struct TGraphAdjacency {
        std::vector<size_t> Offsets;
        std::vector<uint32_t> Vertices;
        std::vector<uint32_t> Degrees;
        explicit TGraphAdjacency(const TGraph& graph)
            : Offsets(graph.OriginalIds.size() + 1)
            , Vertices()
            , Degrees(graph.OriginalIds.size())
        {
            for (const auto& edge : graph.Edges) {
                if (edge.Source != edge.Target) {
                    ++Offsets[edge.Source + 1];
                    ++Offsets[edge.Target + 1];
                    ++Degrees[edge.Source];
                    ++Degrees[edge.Target];
                }
            }
            std::partial_sum(Offsets.begin(), Offsets.end(), Offsets.begin());
            auto positions = Offsets;
            Vertices.resize(Offsets.back());
            for (const auto& edge : graph.Edges) {
                if (edge.Source != edge.Target) {
                    Vertices[positions[edge.Source]++] = edge.Target;
                    Vertices[positions[edge.Target]++] = edge.Source;
                }
            }
        }
    };

    std::vector<uint32_t> BuildVertexOrder(const TGraph& graph,
        const std::string& mode,
        const TGraphAdjacency& adjacency,
        const std::string& histogramDirectory) {
        size_t vertexCount = graph.OriginalIds.size();
        std::vector<uint32_t> order(vertexCount);
        std::iota(order.begin(), order.end(), 0);
        if (mode == "degree") {
            std::sort(order.begin(), order.end(), [&](auto left, auto right) {
                if (adjacency.Degrees[left] != adjacency.Degrees[right]) {
                    return adjacency.Degrees[left] > adjacency.Degrees[right];
                }
                return left < right;
            });
            return order;
        }
        if (mode == "bfs") {
            std::vector<uint32_t> roots = order;
            std::sort(roots.begin(), roots.end(), [&](auto left, auto right) {
                if (adjacency.Degrees[left] != adjacency.Degrees[right]) {
                    return adjacency.Degrees[left] > adjacency.Degrees[right];
                }
                return left < right;
            });
            order.clear();
            std::vector<uint8_t> seen(vertexCount);
            size_t head = 0;
            for (auto root : roots) {
                if (!seen[root]) {
                    seen[root] = 1;
                    order.push_back(root);
                    while (head < order.size()) {
                        auto source = order[head++];
                        for (size_t neighborPosition = adjacency.Offsets[source]; neighborPosition < adjacency.Offsets[source + 1]; neighborPosition++) {
                            auto target = adjacency.Vertices[neighborPosition];
                            if (!seen[target]) {
                                seen[target] = 1;
                                order.push_back(target);
                            }
                        }
                    }
                }
            }
            return order;
        }
        if (mode == "peel") {
            order.clear();
            auto degrees = adjacency.Degrees;
            using TQueueItem = std::pair<uint32_t, uint32_t>;
            std::priority_queue<TQueueItem, std::vector<TQueueItem>, std::greater<TQueueItem>> queue;
            for (size_t i = 0; i < vertexCount; i++) {
                queue.push({degrees[i], i});
            }
            std::vector<uint8_t> gone(vertexCount);
            std::vector<uint32_t> removed;
            unsigned degeneracy = 0;
            while (!queue.empty()) {
                auto [degree, source] = queue.top();
                queue.pop();
                if (gone[source] || degree != degrees[source]) {
                    continue;
                }
                gone[source] = 1;
                order.push_back(source);
                removed.push_back(degree);
                degeneracy = std::max(degeneracy, degree);
                for (size_t neighborPosition = adjacency.Offsets[source]; neighborPosition < adjacency.Offsets[source + 1]; neighborPosition++) {
                    auto target = adjacency.Vertices[neighborPosition];
                    if (!gone[target]) {
                        --degrees[target];
                        queue.push({degrees[target], target});
                    }
                }
            }
            std::cerr << "degeneracy=" << degeneracy << '\n';
            PrintStatistics("removal-degree", removed, histogramDirectory);
            return order;
        }
        return order;
    }

    // All experimental candidates are decoded and checked before their sizes are reported.
    void RunAnalyzer(const TGraph& graph, bool allOrders, bool statisticsOnly, const std::string& histogramDirectory) {
        std::cout << "ordering,component,codec,bytes,bits_per_vertex,bits_per_edge,seconds\n";
        auto printRow = [&](const std::string& ordering, const char* component, const std::string& codec, size_t size, double time) {
            std::cout << ordering << ',' << component << ',' << codec << ',' << size
                      << ',' << 8.0 * size / std::max(size_t(1), graph.OriginalIds.size())
                      << ',' << 8.0 * size / std::max(size_t(1), graph.Edges.size())
                      << ',' << time << '\n';
        };
        std::cerr << "N=" << graph.OriginalIds.size() << " M=" << graph.Edges.size() << "\n";
        PrintStatistics("forward-degree", graph.ForwardDegrees, histogramDirectory);
        std::vector<uint32_t> gaps, weights;
        uint64_t last = 0;
        for (auto vertexId : graph.OriginalIds) {
            gaps.push_back(vertexId - last);
            last = vertexId;
        }
        PrintStatistics("ID-gaps", gaps, histogramDirectory);
        for (const auto& edge : graph.Edges) {
            weights.push_back(edge.Weight);
        }
        PrintStatistics("weights", weights, histogramDirectory);
        weights.clear();
        gaps.clear();
        size_t edgePosition = 0;
        for (size_t source = 0; source < graph.OriginalIds.size(); source++) {
            uint64_t nextMinimum = source;
            for (unsigned j = 0; j < graph.ForwardDegrees[source]; j++) {
                uint64_t target = graph.Edges[edgePosition++].Target;
                gaps.push_back(target - nextMinimum);
                nextMinimum = target + 1;
            }
        }
        PrintStatistics("endpoint-gaps", gaps, histogramDirectory);
        TGraphAdjacency adjacency(graph);
        PrintStatistics("undirected-degree", adjacency.Degrees, histogramDirectory);
        if (statisticsOnly) {
            return;
        }
        auto idSection = SelectIdCodec(graph);
        auto degreeSection = SelectDegreeCodec(graph);
        const char* idNames[] = {"raw32", "delta-varint", "Elias-Fano", "block128-gaps"};
        for (unsigned i = 0; i < 4; i++) {
            double startTime = GetTimeSeconds();
            auto section = EncodeIds(graph.OriginalIds, i);
            Require(DecodeIds(section, graph.OriginalIds.size()) == graph.OriginalIds, "ID roundtrip");
            printRow("sorted", "IDs", idNames[i], section.Data.size(), GetTimeSeconds() - startTime);
        }
        const char* degreeNames[] = {"varint", "gamma-prefix", "canonical-Huffman", "arithmetic", "Rice-k2"};
        for (unsigned i = 0; i < 5; i++) {
            double startTime = GetTimeSeconds();
            auto section = EncodeDegrees(graph.ForwardDegrees, i);
            Require(DecodeDegrees(section, graph.OriginalIds.size()) == graph.ForwardDegrees, "degree roundtrip");
            printRow("sorted", "degrees", degreeNames[i], section.Data.size(), GetTimeSeconds() - startTime);
        }
        const char* adjacencyNames[] = {"fixed", "delta-varint", "Rice-global", "Rice-degree-class", "block128", "interpolative", "Elias-delta", "Rice-BIC-switch", "degree-weighted-BIC", "inferred-weighted-BIC"};
        auto benchmarkAdjacency = [&](const TGraph& orderedGraph, const std::string& ordering, size_t idSize) {
            auto degreeSection = SelectDegreeCodec(orderedGraph);
            double ideal = 0;
            for (size_t source = 0; source < orderedGraph.OriginalIds.size(); source++) {
                double universe = orderedGraph.OriginalIds.size() - source, forwardDegree = orderedGraph.ForwardDegrees[source];
                ideal += (std::lgamma(universe + 1) - std::lgamma(forwardDegree + 1) - std::lgamma(universe - forwardDegree + 1)) / std::log(2.0);
            }
            printRow(ordering, "adjacency", "combinatorial-bound", size_t(std::ceil(ideal / 8)), 0);
            size_t best = SIZE_MAX;
            for (unsigned i = 0; i < 10; i++) {
                double startTime = GetTimeSeconds();
                auto section = EncodeAdjacency(orderedGraph, i);
                auto decodedEdges = DecodeAdjacency(orderedGraph, section, orderedGraph.Edges.size());
                for (size_t j = 0; j < decodedEdges.size(); j++) {
                    Require(decodedEdges[j].Source == orderedGraph.Edges[j].Source && decodedEdges[j].Target == orderedGraph.Edges[j].Target, "adjacency roundtrip");
                }
                printRow(ordering, "adjacency", adjacencyNames[i], section.Data.size(), GetTimeSeconds() - startTime);
                printRow(ordering, "TOTAL", adjacencyNames[i], HEADER_SIZE + idSize + degreeSection.Data.size() + section.Data.size() + orderedGraph.Edges.size(), 0);
                best = std::min(best, section.Data.size());
            }
            if (ordering == "sorted") {
                for (unsigned bias : {0u, 1u, 2u}) {
                    for (unsigned shift : {10u, 12u, 14u, 16u}) {
                        double startTime = GetTimeSeconds();
                        auto section = EncodeInferredAdjacency(orderedGraph, bias, shift);
                        auto decodedEdges = DecodeAdjacency(orderedGraph, section, orderedGraph.Edges.size());
                        for (size_t j = 0; j < decodedEdges.size(); j++) {
                            Require(decodedEdges[j].Source == orderedGraph.Edges[j].Source && decodedEdges[j].Target == orderedGraph.Edges[j].Target, "inferred roundtrip");
                        }
                        std::string name = "inferred-bias" + std::to_string(bias) + "-block" + std::to_string(shift);
                        printRow(ordering, "adjacency", name, section.Data.size(), GetTimeSeconds() - startTime);
                        printRow(ordering, "TOTAL", name, HEADER_SIZE + idSize + degreeSection.Data.size() + section.Data.size() + orderedGraph.Edges.size(), 0);
                        best = std::min(best, section.Data.size());
                    }
                }
            }
            printRow(ordering, "mapping", "raw-or-sorted", idSize, 0);
            printRow(ordering, "degrees", degreeNames[degreeSection.Codec], degreeSection.Data.size(), 0);
            printRow(ordering, "TOTAL", "best", HEADER_SIZE + idSize + degreeSection.Data.size() + best + orderedGraph.Edges.size(), 0);
        };
        benchmarkAdjacency(graph, "sorted", idSection.Data.size());
        if (allOrders) {
            for (std::string orderingMode : {"degree", "bfs", "peel", "construction"}) {
                double startTime = GetTimeSeconds();
                std::vector<uint32_t> ordering;
                if (orderingMode == "construction") {
                    ordering = BuildVertexOrder(graph, "peel", adjacency, histogramDirectory);
                    std::reverse(ordering.begin(), ordering.end());
                } else {
                    ordering = BuildVertexOrder(graph, orderingMode, adjacency, histogramDirectory);
                }
                auto orderedGraph = ReorderGraph(graph, ordering);
                std::cerr << orderingMode << " transformation seconds=" << GetTimeSeconds() - startTime << '\n';
                TBitWriter permutationWriter;
                unsigned bitCount = GetBitWidth(graph.OriginalIds.empty() ? 0 : graph.OriginalIds.size() -1);
                for (auto target : ordering) {
                    permutationWriter.WriteBits(target, bitCount);
                }
                printRow(orderingMode, "mapping", "sorted-set-plus-fixed-permutation", idSection.Data.size() + permutationWriter.GetBytes().size(), 0);
                benchmarkAdjacency(orderedGraph, orderingMode, 4 * graph.OriginalIds.size());
            }
        }
    }

    void TestBitPrimitives() {
        for (uint64_t range : std::vector<uint64_t> {1ULL, 2ULL, 3ULL, 4ULL, 7ULL, 8ULL, 9ULL, 31ULL, 32ULL, 33ULL, UINT64_C(4294967295), UINT64_C(4294967296), UINT64_C(4294967297), UINT64_MAX}) {
            TBitWriter writer;
            std::vector<uint64_t> testValues = {0, range / 2, range - 1};
            for (auto value : testValues) {
                writer.WriteBounded(value, range);
            }
            TBitReader reader{writer.GetBytes()};
            for (auto value : testValues) {
                Require(reader.ReadBounded(range) == value, "bounded primitive");
            }
        }
        for (unsigned bitCount : {0u, 1u, 32u, 63u, 64u}) {
            TBitWriter writer;
            writer.WriteBits(UINT64_MAX, bitCount);
            TBitReader reader{writer.GetBytes()};
            uint64_t expected = bitCount == 64 ? UINT64_MAX : bitCount ? (uint64_t(1) << bitCount) -1 : 0;
            Require(reader.ReadBits(bitCount) == expected, "bits primitive");
        }
        std::cerr << "primitives OK\n";
    }

} // namespace NGraphCompressor
