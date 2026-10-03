#include "codecs.h"
#include "bitstream.h"
#include "byte_io.h"
#include <algorithm>
#include <utility>

namespace NGraphCompressor {

    // Elias-Fano splits each sorted ID into a fixed-width low part and a unary high gap.
    TSection EncodeIds(const std::vector<uint32_t>& originalIds, unsigned codec) {
        TSection section{ uint8_t(codec), {}, "IDs" };
        TBitWriter writer;
        if (codec == ID_RAW32) {
            for (auto vertexId : originalIds) {
                AppendInteger(section.Data, vertexId, 4);
            }
            return section;
        }
        if (codec == ID_DELTA_VARINT) {
            uint64_t previousId = 0;
            for (auto vertexId : originalIds) {
                WriteVarint(section.Data, uint64_t(vertexId) - previousId);
                previousId = vertexId;
            }
            return section;
        }
        if (codec == ID_ELIAS_FANO) {
            unsigned lowBitCount = originalIds.empty() ? 0 : GetBitWidth((uint64_t(1) << 32) / originalIds.size()) -1;
            writer.WriteBits(lowBitCount, 8);
            for (auto vertexId : originalIds) {
                writer.WriteBits(vertexId, lowBitCount);
            }
            uint64_t previousId = 0;
            for (auto vertexId : originalIds) {
                uint64_t highPart = uint64_t(vertexId) >> lowBitCount;
                writer.WriteUnary(highPart - previousId);
                previousId = highPart;
            }
        }
        if (codec == ID_BLOCK_GAPS) {
            uint64_t previousId = 0;
            for (size_t blockBegin = 0; blockBegin < originalIds.size(); blockBegin += 128) {
                size_t blockEnd = std::min(blockBegin + 128, originalIds.size());
                unsigned bitCount = 0;
                uint64_t previousInBlock = previousId;
                for (size_t i = blockBegin; i < blockEnd; i++) {
                    bitCount = std::max(bitCount, GetBitWidth(uint64_t(originalIds[i]) - previousInBlock));
                    previousInBlock = originalIds[i];
                }
                writer.WriteBits(bitCount, 6);
                for (size_t i = blockBegin; i < blockEnd; i++) {
                    writer.WriteBits(uint64_t(originalIds[i]) - previousId, bitCount);
                    previousId = originalIds[i];
                }
            }
        }
        section.Data = writer.TakeBytes();
        return section;
    }

    std::vector<uint32_t> DecodeIds(const TSection& section, size_t vertexCount) {
        std::vector<uint32_t> originalIds(vertexCount);
        size_t position = 0;
        TBitReader reader{section.Data};
        uint64_t previousId = 0;
        if (section.Codec == ID_RAW32) {
            for (auto& vertexId : originalIds) {
                vertexId = uint32_t(ReadInteger(section.Data, position, 4));
            }
        } else if (section.Codec == ID_DELTA_VARINT) {
            for (auto& vertexId : originalIds) {
                previousId += ReadVarint(section.Data, position);
                Require(previousId <= UINT32_MAX, "ID overflow");
                vertexId = previousId;
            }
        } else if (section.Codec == ID_ELIAS_FANO) {
            unsigned lowBitCount = reader.ReadBits(8);
            Require(lowBitCount <= 32, "ID low width");
            for (auto& vertexId : originalIds) {
                vertexId = reader.ReadBits(lowBitCount);
            }
            for (auto& vertexId : originalIds) {
                previousId += reader.ReadUnary();
                uint64_t value = (previousId << lowBitCount) | vertexId;
                Require(value <= UINT32_MAX, "ID overflow");
                vertexId = value;
            }
        } else if (section.Codec == ID_BLOCK_GAPS) {
            for (size_t blockBegin = 0; blockBegin < vertexCount; blockBegin += 128) {
                unsigned bitCount = reader.ReadBits(6);
                Require(bitCount <= 32, "ID width");
                for (size_t i = blockBegin; i < std::min(vertexCount, blockBegin + 128); i++) {
                    previousId += reader.ReadBits(bitCount);
                    Require(previousId <= UINT32_MAX, "ID overflow");
                    originalIds[i] = previousId;
                }
            }
        } else {
            throw std::runtime_error("unknown ID codec");
        }
        if (section.Codec) {
            for (size_t i = 1; i < vertexCount; i++) {
                Require(originalIds[i] > originalIds[i - 1], "unordered IDs");
            }
        }
        return originalIds;
    }

    // The complete section size decides the winner, including its metadata.
    TSection SelectIdCodec(const TGraph& graph) {
        TSection bestSection = EncodeIds(graph.OriginalIds, ID_RAW32);
        for (unsigned codec = 1; codec < ID_CODEC_COUNT; codec++) {
            auto section = EncodeIds(graph.OriginalIds, codec);
            if (section.Data.size() < bestSection.Data.size()) {
                bestSection = std::move(section);
            }
        }
        return bestSection;
    }

} // namespace NGraphCompressor
