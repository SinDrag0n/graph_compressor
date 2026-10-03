#pragma once

#include "common.h"

#include <string>
#include <tuple>

namespace NGraphCompressor {

    // Endpoints are dense IDs after parsing; weight always travels with its edge.
    struct TEdge {
        uint32_t Source = 0;
        uint32_t Target = 0;
        uint8_t Weight = 0;
    };

    struct TGraph {
        std::vector<uint32_t> OriginalIds;
        std::vector<TEdge> Edges;
        std::vector<uint32_t> ForwardDegrees;
    };

    inline bool CompareEdges(const TEdge& left, const TEdge& right) {
        return std::tie(left.Source, left.Target, left.Weight)
            < std::tie(right.Source, right.Target, right.Weight);
    }

    TGraph ReadGraphTsv(const std::string& path);
    void WriteGraphTsv(const std::string& path, const TGraph& graph);

} // namespace NGraphCompressor
