#pragma once

#include "graph.h"

#include <string>

namespace NGraphCompressor {

    constexpr size_t HEADER_SIZE = 60;
    constexpr unsigned FORMAT_VERSION = 1;

    void SerializeGraph(const TGraph& graph, const std::string& outputPath, int forcedCodec = -1);
    TGraph DeserializeGraph(const std::string& inputPath);

} // namespace NGraphCompressor
