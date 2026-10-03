#include "format.h"
#include "byte_io.h"
#include "codecs.h"

#include <array>
#include <iostream>
#include <utility>

#include <cstring>

namespace NGraphCompressor {
    namespace {

        constexpr uint64_t FNV_INITIAL_HASH = 14695981039346656037ULL;
        constexpr uint64_t FNV_PRIME = 1099511628211ULL;

        struct TEncodedSections {
            TSection Ids;
            TSection Degrees;
            TSection Adjacency;
        };

        struct TFileHeader {
            uint32_t VertexCount = 0;
            uint64_t EdgeCount = 0;
            std::array<uint8_t, 3> Codecs = {};
            std::array<uint64_t, 4> SectionSizes = {};
            uint64_t Checksum = 0;
        };

        void UpdateChecksum(uint64_t& checksum, uint8_t value) {
            checksum ^= value;
            checksum *= FNV_PRIME;
        }

        TSection SelectAdjacencyCodec(const TGraph& graph, int forcedCodec) {
            auto bestSection = EncodeAdjacency(graph, forcedCodec >= 0 ? forcedCodec : ADJ_RICE_DEGREE);
            if (forcedCodec >= 0) {
                return bestSection;
            }
            for (unsigned codec : {ADJ_RICE_GLOBAL, ADJ_INTERPOLATIVE, ADJ_RICE_INTERPOLATIVE,
                                   ADJ_FULL_DEGREE_MODEL, ADJ_INFERRED_MODEL}) {
                auto candidate = EncodeAdjacency(graph, codec);
                if (candidate.Data.size() < bestSection.Data.size()) {
                    bestSection = std::move(candidate);
                }
            }
            return bestSection;
        }

        TFileHeader ReadHeader(const TBytes& bytes) {
            Require(bytes.size() >= HEADER_SIZE && std::memcmp(bytes.data(), "GC01", 4) == 0, "bad header");
            size_t position = 4;
            Require(ReadInteger(bytes, position, 1) == FORMAT_VERSION, "unknown version");
            TFileHeader header;
            for (auto& codec : header.Codecs) {
                codec = ReadInteger(bytes, position, 1);
            }
            header.VertexCount = ReadInteger(bytes, position, 4);
            header.EdgeCount = ReadInteger(bytes, position, 8);

            uint64_t remainingBytes = bytes.size() - HEADER_SIZE;
            for (auto& sectionSize : header.SectionSizes) {
                sectionSize = ReadInteger(bytes, position, 8);
                Require(sectionSize <= remainingBytes, "section size");
                remainingBytes -= sectionSize;
            }
            Require(remainingBytes == 0 && header.SectionSizes[3] == header.EdgeCount, "payload size");
            header.Checksum = ReadInteger(bytes, position, 8);
            return header;
        }

    } // anonymous namespace

    void SerializeGraph(const TGraph& graph, const std::string& outputPath, int forcedCodec) {
        const double startTime = GetTimeSeconds();
        const TEncodedSections sections{
            SelectIdCodec(graph),
            SelectDegreeCodec(graph),
            SelectAdjacencyCodec(graph, forcedCodec)
        };

        TBytes bytes = {'G', 'C', '0', '1'};
        AppendInteger(bytes, FORMAT_VERSION, 1);
        AppendInteger(bytes, sections.Ids.Codec, 1);
        AppendInteger(bytes, sections.Degrees.Codec, 1);
        AppendInteger(bytes, sections.Adjacency.Codec, 1);
        AppendInteger(bytes, graph.OriginalIds.size(), 4);
        AppendInteger(bytes, graph.Edges.size(), 8);
        for (uint64_t size : {uint64_t(sections.Ids.Data.size()), uint64_t(sections.Degrees.Data.size()),
                              uint64_t(sections.Adjacency.Data.size()), uint64_t(graph.Edges.size())}) {
            AppendInteger(bytes, size, 8);
        }

        uint64_t checksum = FNV_INITIAL_HASH;
        for (const auto* payload : {&sections.Ids.Data, &sections.Degrees.Data, &sections.Adjacency.Data}) {
            for (uint8_t value : *payload) {
                UpdateChecksum(checksum, value);
            }
        }
        for (const auto& edge : graph.Edges) {
            UpdateChecksum(checksum, edge.Weight);
        }
        AppendInteger(bytes, checksum, 8);
        Require(bytes.size() == HEADER_SIZE, "header size");

        for (const auto* payload : {&sections.Ids.Data, &sections.Degrees.Data, &sections.Adjacency.Data}) {
            bytes.insert(bytes.end(), payload->begin(), payload->end());
        }
        for (const auto& edge : graph.Edges) {
            bytes.push_back(edge.Weight);
        }
        WriteFile(outputPath, bytes);

        std::cerr << "IDs=" << sections.Ids.Data.size() << " mode=" << unsigned(sections.Ids.Codec)
                  << " degrees=" << sections.Degrees.Data.size() << " mode=" << unsigned(sections.Degrees.Codec)
                  << " adjacency=" << sections.Adjacency.Data.size() << " mode=" << unsigned(sections.Adjacency.Codec)
                  << " weights=" << graph.Edges.size() << " header=" << HEADER_SIZE
                  << " TOTAL=" << bytes.size() << " compression_seconds=" << GetTimeSeconds() - startTime << '\n';
    }

    TGraph DeserializeGraph(const std::string& inputPath) {
        const auto bytes = ReadFile(inputPath);
        const auto header = ReadHeader(bytes);
        uint64_t checksum = FNV_INITIAL_HASH;
        for (size_t i = HEADER_SIZE; i < bytes.size(); ++i) {
            UpdateChecksum(checksum, bytes[i]);
        }
        Require(checksum == header.Checksum, "checksum mismatch");
        Require(header.VertexCount <= 2 * header.EdgeCount &&
                header.VertexCount <= SIZE_MAX / sizeof(uint32_t), "invalid vertex count");

        TEncodedSections sections;
        size_t position = HEADER_SIZE;
        unsigned sectionIndex = 0;
        for (auto* section : {&sections.Ids, &sections.Degrees, &sections.Adjacency}) {
            section->Codec = header.Codecs[sectionIndex];
            const uint64_t sectionSize = header.SectionSizes[sectionIndex++];
            section->Data.assign(bytes.begin() + position, bytes.begin() + position + sectionSize);
            position += sectionSize;
        }

        TGraph graph;
        graph.OriginalIds = DecodeIds(sections.Ids, header.VertexCount);
        graph.ForwardDegrees = DecodeDegrees(sections.Degrees, header.VertexCount);
        graph.Edges = DecodeAdjacency(graph, sections.Adjacency, header.EdgeCount);
        for (auto& edge : graph.Edges) {
            edge.Weight = bytes[position++];
        }
        return graph;
    }

} // namespace NGraphCompressor
