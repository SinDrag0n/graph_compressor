#pragma once

#include "graph.h"

#include <string>

namespace NGraphCompressor {

    constexpr unsigned ID_RAW32 = 0;
    constexpr unsigned ID_DELTA_VARINT = 1;
    constexpr unsigned ID_ELIAS_FANO = 2;
    constexpr unsigned ID_BLOCK_GAPS = 3;
    constexpr unsigned ID_CODEC_COUNT = 4;

    constexpr unsigned DEGREE_VARINT = 0;
    constexpr unsigned DEGREE_GAMMA = 1;
    constexpr unsigned DEGREE_HUFFMAN = 2;
    constexpr unsigned DEGREE_ARITHMETIC = 3;
    constexpr unsigned DEGREE_RICE = 4;
    constexpr unsigned DEGREE_CODEC_COUNT = 5;

    constexpr unsigned ADJ_FIXED = 0;
    constexpr unsigned ADJ_VARINT = 1;
    constexpr unsigned ADJ_RICE_GLOBAL = 2;
    constexpr unsigned ADJ_RICE_DEGREE = 3;
    constexpr unsigned ADJ_BLOCK_GAPS = 4;
    constexpr unsigned ADJ_INTERPOLATIVE = 5;
    constexpr unsigned ADJ_ELIAS_DELTA = 6;
    constexpr unsigned ADJ_RICE_INTERPOLATIVE = 7;
    constexpr unsigned ADJ_FULL_DEGREE_MODEL = 8;
    constexpr unsigned ADJ_INFERRED_MODEL = 9;
    constexpr unsigned ADJ_CODEC_COUNT = 10;

    constexpr unsigned INFERRED_DEGREE_BIAS = 2;
    constexpr unsigned INFERRED_BLOCK_SHIFT = 10;
    constexpr uint64_t MODEL_MINIMUM_DEGREE = 5;
    constexpr uint64_t MODEL_MAXIMUM_DEGREE = uint64_t(1) << 20;
    constexpr uint64_t MODEL_FREQUENCY_LIMIT = uint64_t(1) << 23;

    // Codec is a section-specific wire-format ID, not a shared enum of all codecs.
    struct TSection {
        uint8_t Codec = 0;
        TBytes Data;
        std::string Name;
    };

    TSection EncodeIds(const std::vector<uint32_t>& originalIds, unsigned codec);
    std::vector<uint32_t> DecodeIds(const TSection& section, size_t vertexCount);
    TSection SelectIdCodec(const TGraph& graph);

    TSection EncodeDegrees(const std::vector<uint32_t>& degrees, unsigned codec);
    std::vector<uint32_t> DecodeDegrees(const TSection& section, size_t vertexCount);
    TSection SelectDegreeCodec(const TGraph& graph);

    TSection EncodeAdjacency(const TGraph& graph, unsigned codec);
    std::vector<TEdge> DecodeAdjacency(const TGraph& graph, const TSection& section, uint64_t edgeCount);
    TSection EncodeInferredAdjacency(const TGraph& graph, unsigned bias,
                                    unsigned blockShift = INFERRED_BLOCK_SHIFT);

} // namespace NGraphCompressor
