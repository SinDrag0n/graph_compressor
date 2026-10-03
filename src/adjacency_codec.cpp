#include "codecs.h"
#include "arithmetic.h"
#include "byte_io.h"
#include <algorithm>
#include <array>
#include <limits>
#include <utility>

namespace NGraphCompressor {

    // Sorted distinct neighbors reserve room for the elements on both sides of the median.
    // Bounds [source, N) include source itself, which represents a self-loop.
    void WriteInterpolative(TBitWriter& writer,
        const std::vector<TEdge>& edges,
        size_t begin,
        size_t end,
        uint64_t lowerBound,
        uint64_t upperBound) {
        if (begin == end) {
            return;
        }
        size_t middle = begin + (end - begin) / 2;
        uint64_t minimum = lowerBound + middle - begin;
        uint64_t maximum = upperBound - (end - middle - 1);
        Require(minimum <= maximum, "BIC bounds");
        uint64_t neighbor = edges[middle].Target;
        writer.WriteBounded(neighbor - minimum, maximum - minimum + 1);
        WriteInterpolative(writer, edges, begin, middle, lowerBound, neighbor ? neighbor - 1 : 0);
        WriteInterpolative(writer, edges, middle + 1, end, neighbor + 1, upperBound);
    }

    void ReadInterpolative(TBitReader& reader,
        std::vector<TEdge>& edges,
        size_t begin,
        size_t end,
        uint64_t lowerBound,
        uint64_t upperBound,
        uint32_t source) {
        if (begin == end) {
            return;
        }
        size_t middle = begin + (end - begin) / 2;
        uint64_t minimum = lowerBound + middle - begin;
        uint64_t maximum = upperBound - (end - middle - 1);
        Require(minimum <= maximum, "BIC bounds");
        uint64_t neighbor = minimum + reader.ReadBounded(maximum - minimum + 1);
        edges[middle].Source = source;
        edges[middle].Target = neighbor;
        ReadInterpolative(reader, edges, begin, middle, lowerBound, neighbor ? neighbor - 1 : 0, source);
        ReadInterpolative(reader, edges, middle + 1, end, neighbor + 1, upperBound, source);
    }

    struct TWeightedModel {
        std::vector<uint32_t> Cumulative;
        explicit TWeightedModel(const std::vector<uint32_t>& estimates)
            : Cumulative(estimates.size() + 1)
        {
            uint64_t totalFrequency = 0;
            for (auto estimate : estimates) {
                totalFrequency += std::max(1u, estimate);
            }
            uint64_t scale = std::max(uint64_t(1), (totalFrequency + MODEL_FREQUENCY_LIMIT - 1) / MODEL_FREQUENCY_LIMIT);
            for (size_t i = 0; i < estimates.size(); ++i) {
                uint64_t value = uint64_t(Cumulative[i]) + std::max(uint64_t(1), uint64_t(estimates[i]) / scale);
                Require(value <= UINT32_MAX, "weighted model total");
                Cumulative[i + 1] = value;
            }
        }
    };

    std::vector<uint32_t> CalculateFullDegrees(const TGraph& graph) {
        std::vector<uint32_t> degrees(graph.OriginalIds.size());
        for (const auto& edge : graph.Edges) {
            ++degrees[edge.Source];
            if (edge.Target != edge.Source) {
                ++degrees[edge.Target];
            }
        }
        return degrees;
    }

    void WriteWeightedInterpolative(TArithmeticEncoder& writer,
        const TWeightedModel& model,
        const std::vector<TEdge>& edges,
        size_t begin,
        size_t end,
        uint32_t lowerBound,
        uint32_t upperBound) {
        if (begin == end) {
            return;
        }
        size_t middle = begin + (end - begin) / 2;
        uint32_t minimum = lowerBound + middle - begin;
        uint32_t maximum = upperBound - (end - middle - 1);
        uint32_t neighbor = edges[middle].Target;
        uint32_t baseFrequency = model.Cumulative[minimum];
        uint32_t totalFrequency = model.Cumulative[size_t(maximum) + 1] - baseFrequency;
        writer.ProcessSymbol(model.Cumulative[neighbor] - baseFrequency, model.Cumulative[size_t(neighbor) + 1] - baseFrequency, totalFrequency);
        if (begin < middle) {
            WriteWeightedInterpolative(writer, model, edges, begin, middle, lowerBound, neighbor - 1);
        }
        if (middle + 1 < end) {
            WriteWeightedInterpolative(writer, model, edges, middle + 1, end, neighbor + 1, upperBound);
        }
    }

    void ReadWeightedInterpolative(TArithmeticDecoder& reader,
        const TWeightedModel& model,
        std::vector<TEdge>& edges,
        size_t begin,
        size_t end,
        uint32_t lowerBound,
        uint32_t upperBound,
        uint32_t source) {
        if (begin == end) {
            return;
        }
        size_t middle = begin + (end - begin) / 2;
        uint32_t minimum = lowerBound + middle - begin;
        uint32_t maximum = upperBound - (end - middle - 1);
        Require(minimum <= maximum, "weighted BIC bounds");
        uint32_t baseFrequency = model.Cumulative[minimum];
        uint32_t totalFrequency = model.Cumulative[size_t(maximum) + 1] - baseFrequency;
        uint32_t targetFrequency = reader.GetTarget(totalFrequency) + baseFrequency;
        auto neighborPosition = std::upper_bound(model.Cumulative.begin() + minimum, model.Cumulative.begin() + size_t(maximum) + 2, targetFrequency);
        uint32_t neighbor = neighborPosition - model.Cumulative.begin() -1;
        Require(neighbor >= minimum && neighbor <= maximum, "weighted endpoint");
        reader.ProcessSymbol(model.Cumulative[neighbor] - baseFrequency, model.Cumulative[size_t(neighbor) + 1] - baseFrequency, totalFrequency);
        edges[middle] = {source, neighbor, 0};
        if (begin < middle) {
            ReadWeightedInterpolative(reader, model, edges, begin, middle, lowerBound, neighbor - 1, source);
        }
        if (middle + 1 < end) {
            ReadWeightedInterpolative(reader, model, edges, middle + 1, end, neighbor + 1, upperBound, source);
        }
    }

    // These are probability estimates, not graph constraints. Every endpoint keeps positive mass.
    // Both processes know forward degrees and the backward counts from completed rows.
    std::vector<uint32_t> EstimateDegrees(const TGraph& graph,
        const std::vector<uint32_t>& backwardCounts,
        size_t blockBegin,
        unsigned bias) {
        std::vector<uint32_t> estimatedDegrees(graph.OriginalIds.size());
        uint64_t vertexCount = graph.OriginalIds.size();
        for (size_t target = 0; target < estimatedDegrees.size(); target++) {
            uint64_t knownDegree = uint64_t(graph.ForwardDegrees[target]) + backwardCounts[target];
            uint64_t numerator = std::max(uint64_t(1), knownDegree > bias ? knownDegree - bias : 0) * vertexCount;
            uint64_t observedFractionNumerator = vertexCount - target + blockBegin;
            estimatedDegrees[target] = std::min(MODEL_MAXIMUM_DEGREE, std::max(MODEL_MINIMUM_DEGREE, numerator / observedFractionNumerator));
        }
        return estimatedDegrees;
    }

    // Rebuild the model at block boundaries and freeze it while coding the block.
    // Update backward counts only after an entire row has been encoded or decoded.
    TSection EncodeInferredAdjacency(const TGraph& graph, unsigned bias, unsigned blockShift) {
        TSection section{ 9, {}, "inferred-weighted-BIC" };
        section.Data.push_back(bias);
        section.Data.push_back(blockShift);
        TArithmeticEncoder arithmeticCoder;
        std::vector<uint32_t> backwardCounts(graph.OriginalIds.size());
        TWeightedModel model(EstimateDegrees(graph, backwardCounts, 0, bias));
        size_t edgePosition = 0;
        for (size_t source = 0; source < graph.OriginalIds.size(); source++) {
            if (source && !(source % (size_t(1) << blockShift))) {
                model = TWeightedModel(EstimateDegrees(graph, backwardCounts, source, bias));
            }
            size_t end = edgePosition + graph.ForwardDegrees[source];
            WriteWeightedInterpolative(arithmeticCoder, model, graph.Edges, edgePosition, end, source, graph.OriginalIds.size() -1);
            while (edgePosition < end) {
                ++backwardCounts[graph.Edges[edgePosition++].Target];
            }
        }
        auto encodedBits = arithmeticCoder.Finish();
        section.Data.insert(section.Data.end(), encodedBits.begin(), encodedBits.end());
        return section;
    }

    std::vector<TEdge> DecodeInferredAdjacency(const TGraph& graph, const TSection& section, uint64_t edgeCount) {
        Require(section.Data.size() >= 2 && section.Data[0] <= 2 && section.Data[1] >= 8 && section.Data[1] <= 20, "inferred model variant");
        unsigned bias = section.Data[0], blockShift = section.Data[1];
        TBytes encodedBits(section.Data.begin() + 2, section.Data.end());
        TArithmeticDecoder arithmeticCoder(encodedBits);
        std::vector<TEdge> edges(edgeCount);
        std::vector<uint32_t> backwardCounts(graph.OriginalIds.size());
        TWeightedModel model(EstimateDegrees(graph, backwardCounts, 0, bias));
        size_t edgePosition = 0;
        for (size_t source = 0; source < graph.OriginalIds.size(); source++) {
            if (source && !(source % (size_t(1) << blockShift))) {
                model = TWeightedModel(EstimateDegrees(graph, backwardCounts, source, bias));
            }
            uint64_t end = edgePosition + uint64_t(graph.ForwardDegrees[source]);
            Require(end <= edgeCount && graph.ForwardDegrees[source] <= graph.OriginalIds.size() - source, "inferred degree count");
            ReadWeightedInterpolative(arithmeticCoder, model, edges, edgePosition, end, source, graph.OriginalIds.size() -1, source);
            while (edgePosition < end) {
                ++backwardCounts[edges[edgePosition++].Target];
            }
        }
        Require(edgePosition == edgeCount, "inferred edge count");
        return edges;
    }

    unsigned GetDegreeClass(uint32_t degree) {
        return GetBitWidth(degree);
    }

    TSection EncodeAdjacency(const TGraph& graph, unsigned codec) {
        if (codec == ADJ_INFERRED_MODEL) {
            return EncodeInferredAdjacency(graph, INFERRED_DEGREE_BIAS);
        }
        TSection section{ uint8_t(codec), {}, "adjacency" };
        TBitWriter writer;
        if (codec == ADJ_FULL_DEGREE_MODEL) {
            auto full = CalculateFullDegrees(graph);
            TSection degreeSection = EncodeDegrees(full, 3);
            section.Data.push_back(degreeSection.Codec);
            WriteVarint(section.Data, degreeSection.Data.size());
            section.Data.insert(section.Data.end(), degreeSection.Data.begin(), degreeSection.Data.end());
            TWeightedModel model(full);
            TArithmeticEncoder arithmeticCoder;
            size_t edgePosition = 0;
            for (size_t source = 0; source < graph.OriginalIds.size(); source++) {
                WriteWeightedInterpolative(arithmeticCoder, model, graph.Edges, edgePosition, edgePosition + graph.ForwardDegrees[source], source, graph.OriginalIds.size() -1);
                edgePosition += graph.ForwardDegrees[source];
            }
            auto encodedBits = arithmeticCoder.Finish();
            section.Data.insert(section.Data.end(), encodedBits.begin(), encodedBits.end());
            return section;
        }
        unsigned endpointWidth = GetBitWidth(graph.OriginalIds.empty() ? 0 : graph.OriginalIds.size() -1);
        std::array<unsigned, 33> riceParameters {};
        if (codec == ADJ_RICE_GLOBAL || codec == ADJ_RICE_DEGREE || codec == ADJ_RICE_INTERPOLATIVE) {
            std::array<std::array<uint64_t, 32>, 33> costs {};
            size_t edgePosition = 0;
            for (size_t source = 0; source < graph.OriginalIds.size(); source++) {
                unsigned degreeClass = codec == ADJ_RICE_GLOBAL ? 0 : GetDegreeClass(graph.ForwardDegrees[source]);
                uint64_t nextMinimum = source;
                for (unsigned j = 0; j < graph.ForwardDegrees[source]; j++) {
                    uint64_t target = graph.Edges[edgePosition++].Target, gap = target - nextMinimum;
                    nextMinimum = target + 1;
                    for (unsigned riceParameter = 0; riceParameter < 32; riceParameter++) {
                        costs[degreeClass][riceParameter] += (gap >> riceParameter) + 1 + riceParameter;
                    }
                }
            }
            for (unsigned degreeClass = 0; degreeClass < 33; degreeClass++) {
                riceParameters[degreeClass] = std::min_element(costs[degreeClass].begin(), costs[degreeClass].end()) - costs[degreeClass].begin();
            }
            if (codec == ADJ_RICE_GLOBAL) {
                writer.WriteBits(riceParameters[0], 8);
            } else {
                for (auto riceParameter : riceParameters) {
                    writer.WriteBits(riceParameter, 5);
                }
            }
        }
        size_t edgePosition = 0;
        for (size_t source = 0; source < graph.OriginalIds.size(); source++) {
            unsigned forwardDegree = graph.ForwardDegrees[source];
            if (codec == ADJ_INTERPOLATIVE) {
                WriteInterpolative(writer, graph.Edges, edgePosition, edgePosition + forwardDegree, source, graph.OriginalIds.size() -1);
                edgePosition += forwardDegree;
                continue;
            }
            uint64_t nextMinimum = source;
            unsigned riceParameter = riceParameters[codec == ADJ_RICE_GLOBAL ? 0 : GetDegreeClass(forwardDegree)];
            if (codec == ADJ_RICE_INTERPOLATIVE) {
                TBitWriter trialWriter;
                WriteInterpolative(trialWriter, graph.Edges, edgePosition, edgePosition + forwardDegree, source, graph.OriginalIds.size() -1);
                uint64_t riceBitCount = 0, previousInBlock = source;
                for (size_t i = edgePosition; i < edgePosition + forwardDegree; i++) {
                    uint64_t gap = graph.Edges[i].Target - previousInBlock;
                    previousInBlock = uint64_t(graph.Edges[i].Target) + 1;
                    riceBitCount += (gap >> riceParameter) + 1 + riceParameter;
                }
                bool useInterpolative = trialWriter.GetBitCount() < riceBitCount;
                writer.WriteBit(useInterpolative);
                if (useInterpolative) {
                    WriteInterpolative(writer, graph.Edges, edgePosition, edgePosition + forwardDegree, source, graph.OriginalIds.size() -1);
                    edgePosition += forwardDegree;
                    continue;
                }
            }
            if (codec == ADJ_BLOCK_GAPS) {
                for (size_t blockStart = edgePosition; blockStart < edgePosition + forwardDegree; blockStart += 128) {
                    size_t end = std::min(edgePosition + forwardDegree, blockStart + 128);
                    unsigned blockWidth = 0;
                    uint64_t previousInBlock = nextMinimum;
                    for (size_t i = blockStart; i < end; i++) {
                        blockWidth = std::max(blockWidth, GetBitWidth(uint64_t(graph.Edges[i].Target) - previousInBlock));
                        previousInBlock = uint64_t(graph.Edges[i].Target) + 1;
                    }
                    writer.WriteBits(blockWidth, 6);
                    for (size_t i = blockStart; i < end; i++) {
                        writer.WriteBits(uint64_t(graph.Edges[i].Target) - nextMinimum, blockWidth);
                        nextMinimum = uint64_t(graph.Edges[i].Target) + 1;
                    }
                }
                edgePosition += forwardDegree;
                continue;
            }
            for (unsigned j = 0; j < forwardDegree; j++) {
                uint32_t target = graph.Edges[edgePosition++].Target;
                uint64_t gap = target - nextMinimum;
                nextMinimum = uint64_t(target) + 1;
                if (codec == ADJ_FIXED) {
                    writer.WriteBits(target, endpointWidth);
                } else if (codec == ADJ_VARINT) {
                    WriteVarint(section.Data, gap);
                } else if (codec == ADJ_RICE_GLOBAL || codec == ADJ_RICE_DEGREE || codec == ADJ_RICE_INTERPOLATIVE) {
                    writer.WriteRice(gap, riceParameter);
                } else if (codec == ADJ_ELIAS_DELTA) {
                    writer.WriteDelta(gap + 1);
                } else {
                    throw std::runtime_error("adjacency mode");
                }
            }
        }
        if (codec != ADJ_VARINT) {
            section.Data = writer.TakeBytes();
        }
        return section;
    }

    std::vector<TEdge> DecodeAdjacency(const TGraph& graph, const TSection& section, uint64_t edgeCount) {
        Require(edgeCount <= std::numeric_limits<size_t>::max() / sizeof(TEdge), "edge size overflow");
        if (section.Codec == ADJ_INFERRED_MODEL) {
            return DecodeInferredAdjacency(graph, section, edgeCount);
        }
        std::vector<TEdge> edges(static_cast<size_t> (edgeCount));
        if (section.Codec == ADJ_FULL_DEGREE_MODEL) {
            size_t edgePosition = 0;
            TSection degreeSection;
            degreeSection.Codec = ReadInteger(section.Data, edgePosition, 1);
            uint64_t size = ReadVarint(section.Data, edgePosition);
            Require(size <= section.Data.size() - edgePosition, "weighted degree size");
            degreeSection.Data.assign(section.Data.begin() + edgePosition, section.Data.begin() + edgePosition + size);
            edgePosition += size;
            auto full = DecodeDegrees(degreeSection, graph.OriginalIds.size());
            TWeightedModel model(full);
            TBytes encodedBits(section.Data.begin() + edgePosition, section.Data.end());
            TArithmeticDecoder arithmeticCoder(encodedBits);
            edgePosition = 0;
            for (size_t source = 0; source < graph.OriginalIds.size(); source++) {
                uint64_t end = edgePosition + uint64_t(graph.ForwardDegrees[source]);
                Require(end <= edgeCount && graph.ForwardDegrees[source] <= graph.OriginalIds.size() - source, "weighted degree count");
                ReadWeightedInterpolative(arithmeticCoder, model, edges, edgePosition, end, source, graph.OriginalIds.size() -1, source);
                edgePosition = end;
            }
            Require(edgePosition == edgeCount, "weighted edge count");
            return edges;
        }
        TBitReader reader{section.Data};
        size_t edgePosition = 0, bytePosition = 0;
        unsigned endpointWidth = GetBitWidth(graph.OriginalIds.empty() ? 0 : graph.OriginalIds.size() -1);
        std::array<unsigned, 33> riceParameters {};
        if (section.Codec == ADJ_RICE_GLOBAL) {
            riceParameters[0] = reader.ReadBits(8);
        }
        if (section.Codec == ADJ_RICE_DEGREE || section.Codec == ADJ_RICE_INTERPOLATIVE) {
            for (auto& riceParameter : riceParameters) {
                riceParameter = reader.ReadBits(5);
            }
        }
        for (auto riceParameter : riceParameters) {
            Require(riceParameter < 32, "Rice parameter");
        }
        for (size_t source = 0; source < graph.OriginalIds.size(); source++) {
            uint64_t forwardDegree = graph.ForwardDegrees[source];
            Require(forwardDegree <= graph.OriginalIds.size() - source && forwardDegree <= edgeCount - edgePosition, "invalid degree");
            if (section.Codec == ADJ_INTERPOLATIVE || (section.Codec == ADJ_RICE_INTERPOLATIVE && reader.ReadBit())) {
                ReadInterpolative(reader, edges, edgePosition, edgePosition + forwardDegree, source, graph.OriginalIds.size() -1, source);
                edgePosition += forwardDegree;
                continue;
            }
            uint64_t nextMinimum = source;
            unsigned riceParameter = riceParameters[section.Codec == ADJ_RICE_GLOBAL ? 0 : GetDegreeClass(forwardDegree)];
            if (section.Codec == ADJ_BLOCK_GAPS) {
                for (size_t blockStart = 0; blockStart < forwardDegree; blockStart += 128) {
                    unsigned blockWidth = reader.ReadBits(6);
                    Require(blockWidth <= 32, "adjacency width");
                    for (size_t j = blockStart; j < std::min(forwardDegree, blockStart + 128); j++) {
                        uint64_t target = nextMinimum + reader.ReadBits(blockWidth);
                        Require(target < graph.OriginalIds.size(), "endpoint overflow");
                        edges[edgePosition++] = {uint32_t(source), uint32_t(target), 0};
                        nextMinimum = target + 1;
                    }
                }
                continue;
            }
            for (size_t j = 0; j < forwardDegree; j++) {
                uint64_t target;
                if (section.Codec == ADJ_FIXED) {
                    target = reader.ReadBits(endpointWidth);
                } else if (section.Codec == ADJ_VARINT) {
                    target = nextMinimum + ReadVarint(section.Data, bytePosition);
                } else if (section.Codec == ADJ_RICE_GLOBAL || section.Codec == ADJ_RICE_DEGREE || section.Codec == ADJ_RICE_INTERPOLATIVE) {
                    target = nextMinimum + reader.ReadRice(riceParameter);
                } else if (section.Codec == ADJ_ELIAS_DELTA) {
                    uint64_t bitCount = reader.ReadGamma();
                    Require(bitCount >= 1 && bitCount <= 33, "delta length");
                    target = nextMinimum + ((uint64_t(1) << (bitCount - 1)) | reader.ReadBits(bitCount - 1)) -1;
                } else {
                    throw std::runtime_error("adjacency codec");
                }
                Require(target >= nextMinimum && target < graph.OriginalIds.size(), "endpoint overflow/order");
                edges[edgePosition++] = {uint32_t(source), uint32_t(target), 0};
                nextMinimum = target + 1;
            }
        }
        Require(edgePosition == edgeCount, "edge count mismatch");
        return edges;
    }

} // namespace NGraphCompressor
