#pragma once

#include "graph.h"

#include <string>

namespace NGraphCompressor {

    // Development-only interface; production compression does not run the analyzer.
    void RunAnalyzer(const TGraph& graph, bool allOrders, bool statisticsOnly,
                     const std::string& histogramDirectory);
    void TestBitPrimitives();

} // namespace NGraphCompressor
