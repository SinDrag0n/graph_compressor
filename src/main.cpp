#include "analyzer.h"
#include "format.h"

#include <iostream>
#include <string>

namespace NGraphCompressor {
    namespace {

        struct TOptions {
            std::string Operation;
            std::string InputPath;
            std::string OutputPath;
            std::string HistogramDirectory;
            bool AllOrders = false;
            bool StatisticsOnly = false;
            int ForcedAdjacencyCodec = -1;
        };

        TOptions ParseOptions(int argumentCount, char** arguments) {
            TOptions options;
            for (int i = 1; i < argumentCount; ++i) {
                const std::string argument = arguments[i];
                if (argument == "-s" || argument == "-d" || argument == "--analyze" || argument == "--selftest") {
                    options.Operation = argument;
                } else if (argument == "-i" && i + 1 < argumentCount) {
                    options.InputPath = arguments[++i];
                } else if (argument == "-o" && i + 1 < argumentCount) {
                    options.OutputPath = arguments[++i];
                } else if (argument == "--all-orders") {
                    options.AllOrders = true;
                } else if (argument == "--stats-only") {
                    options.StatisticsOnly = true;
                } else if (argument == "--histograms" && i + 1 < argumentCount) {
                    options.HistogramDirectory = arguments[++i];
                } else if (argument == "--adj" && i + 1 < argumentCount) {
                    options.ForcedAdjacencyCodec = std::stoi(arguments[++i]);
                    Require(options.ForcedAdjacencyCodec >= 0 && options.ForcedAdjacencyCodec <= 9, "adj mode");
                } else {
                    throw std::runtime_error("unknown argument");
                }
            }
            return options;
        }

        void Run(const TOptions& options) {
            if (options.Operation == "--selftest") {
                TestBitPrimitives();
                return;
            }
            Require(!options.InputPath.empty(), "input required");
            const double startTime = GetTimeSeconds();
            if (options.Operation == "-d") {
                Require(!options.OutputPath.empty() && options.InputPath != options.OutputPath,
                        "output required/different path");
                const auto graph = DeserializeGraph(options.InputPath);
                std::cerr << "binary_read_decode_seconds=" << GetTimeSeconds() - startTime << '\n';
                const double outputStartTime = GetTimeSeconds();
                WriteGraphTsv(options.OutputPath, graph);
                std::cerr << "TSV_format_write_seconds=" << GetTimeSeconds() - outputStartTime << '\n';
            } else {
                const auto graph = ReadGraphTsv(options.InputPath);
                std::cerr << "parse/map/sort seconds=" << GetTimeSeconds() - startTime << '\n';
                if (options.Operation == "--analyze") {
                    RunAnalyzer(graph, options.AllOrders, options.StatisticsOnly, options.HistogramDirectory);
                } else {
                    Require(options.Operation == "-s" && !options.OutputPath.empty() &&
                            options.InputPath != options.OutputPath, "usage: run -s|-d -i file -o file");
                    SerializeGraph(graph, options.OutputPath, options.ForcedAdjacencyCodec);
                }
            }
            std::cerr << "total_seconds=" << GetTimeSeconds() - startTime << '\n';
        }

    } // anonymous namespace
} // namespace NGraphCompressor

// The process entry point uses the name mandated by the C++ runtime.
int main(int argumentCount, char** arguments) {
    try {
        NGraphCompressor::Run(NGraphCompressor::ParseOptions(argumentCount, arguments));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
