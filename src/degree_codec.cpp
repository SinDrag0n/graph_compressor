#include "codecs.h"
#include "arithmetic.h"
#include "byte_io.h"
#include "huffman.h"
#include <algorithm>
#include <array>
#include <utility>

namespace NGraphCompressor {

    // Symbols 0..255 are literal degrees; 256 escapes to a separate large-degree stream.
    TSection EncodeDegrees(const std::vector<uint32_t>& degrees, unsigned codec) {
        TSection section{ uint8_t(codec), {}, "degrees" };
        TBitWriter writer;
        TBytes escapeBytes;
        std::array<uint32_t, 257> frequencies {};
        for (auto degree : degrees) {
            ++frequencies[std::min(degree, 256u)];
        }
        if (codec == DEGREE_VARINT) {
            for (auto degree : degrees) {
                WriteVarint(section.Data, degree);
            }
            return section;
        }
        if (codec == DEGREE_GAMMA) {
            for (auto degree : degrees) {
                writer.WriteGamma(uint64_t(degree) + 1);
            }
            section.Data = writer.TakeBytes();
            return section;
        }
        if (codec == DEGREE_HUFFMAN) {
            THuffmanModel huffman(frequencies);
            for (auto length : huffman.Lengths) {
                section.Data.push_back(length);
            }
            for (auto degree : degrees) {
                unsigned symbol = std::min(degree, 256u);
                writer.WriteBits(huffman.Codes[symbol], huffman.Lengths[symbol]);
                if (symbol == 256) {
                    writer.WriteGamma(uint64_t(degree) - 255);
                }
            }
            section.Data.insert(section.Data.end(), writer.GetBytes().begin(), writer.GetBytes().end());
            return section;
        }
        if (codec == DEGREE_ARITHMETIC) {
            std::array<uint32_t, 258> cumulative {};
            for (unsigned i = 0; i < 257; i++) {
                const uint32_t scale = std::max(1u, uint32_t(degrees.size() / (1 << 22)));
                frequencies[i] = frequencies[i] ? std::max(1u, frequencies[i] / scale) : 0;
                WriteVarint(section.Data, frequencies[i]);
                cumulative[i + 1] = cumulative[i] + frequencies[i];
            }
            TArithmeticEncoder coder;
            for (auto degree : degrees) {
                unsigned symbol = std::min(degree, 256u);
                coder.ProcessSymbol(cumulative[symbol], cumulative[symbol + 1], cumulative[257]);
                if (symbol == 256) {
                    WriteVarint(escapeBytes, degree - 256);
                }
            }
            auto encodedBits = coder.Finish();
            WriteVarint(section.Data, escapeBytes.size());
            section.Data.insert(section.Data.end(), escapeBytes.begin(), escapeBytes.end());
            section.Data.insert(section.Data.end(), encodedBits.begin(), encodedBits.end());
            return section;
        }
        if (codec == DEGREE_RICE) {
            unsigned riceParameter = 2;
            writer.WriteBits(riceParameter, 8);
            for (auto degree : degrees) {
                writer.WriteRice(degree, riceParameter);
            }
            section.Data = writer.TakeBytes();
            return section;
        }
        throw std::runtime_error("unknown degree mode");
    }

    std::vector<uint32_t> DecodeDegrees(const TSection& section, size_t vertexCount) {
        std::vector<uint32_t> degrees(vertexCount);
        size_t position = 0;
        if (section.Codec == DEGREE_VARINT) {
            for (auto& degree : degrees) {
                auto decodedValue = ReadVarint(section.Data, position);
                Require(decodedValue <= UINT32_MAX, "degree overflow");
                degree = decodedValue;
            }
            return degrees;
        }
        if (section.Codec == DEGREE_ARITHMETIC) {
            std::array<uint32_t, 258> cumulative {};
            for (unsigned i = 0; i < 257; i++) {
                uint64_t frequency = ReadVarint(section.Data, position);
                Require(frequency <= (1 << 24) && uint64_t(cumulative[i]) + frequency <= (1 << 24), "degree model overflow");
                cumulative[i + 1] = cumulative[i] + frequency;
            }
            uint64_t escapeSize = ReadVarint(section.Data, position);
            Require(escapeSize <= section.Data.size() - position, "degree extra size");
            TBytes escapeBytes(section.Data.begin() + position, section.Data.begin() + position + escapeSize);
            position += escapeSize;
            TBytes encodedBytes(section.Data.begin() + position, section.Data.end());
            TArithmeticDecoder coder(encodedBytes);
            size_t escapePosition = 0;
            Require(cumulative[257] || !vertexCount, "empty degree model");
            for (auto& degree : degrees) {
                uint32_t symbol = coder.GetTarget(cumulative[257]);
                unsigned decodedValue = std::upper_bound(cumulative.begin(), cumulative.end(), symbol) - cumulative.begin() -1;
                Require(decodedValue < 257, "degree symbol");
                coder.ProcessSymbol(cumulative[decodedValue], cumulative[decodedValue + 1], cumulative[257]);
                uint64_t decodedDegree = decodedValue == 256 ? 256 + ReadVarint(escapeBytes, escapePosition) : decodedValue;
                Require(decodedDegree <= UINT32_MAX, "degree overflow");
                degree = decodedDegree;
            }
            return degrees;
        }
        if (section.Codec == DEGREE_HUFFMAN) {
            std::array<unsigned, 257> codeLengths {};
            for (auto& length : codeLengths) {
                length = ReadInteger(section.Data, position, 1);
            }
            auto huffman = THuffmanModel::FromLengths(codeLengths);
            TBytes encodedBytes(section.Data.begin() + position, section.Data.end());
            TBitReader reader{encodedBytes};
            struct TNode {
                int Children[2] = {-1, -1};
                int Symbol = -1;
            };
            std::vector<TNode> tree(1);
            for (unsigned i = 0; i < 257; i++) {
                if (codeLengths[i]) {
                    int node = 0;
                    for (unsigned bitIndex = codeLengths[i]; bitIndex; bitIndex--) {
                        Require(tree[node].Symbol < 0, "Huffman prefix");
                        unsigned bitValue = (huffman.Codes[i] >> (bitIndex - 1)) & 1;
                        int nextNode = tree[node].Children[bitValue];
                        if (nextNode < 0) {
                            nextNode = tree.size();
                            tree[node].Children[bitValue] = nextNode;
                            tree.emplace_back();
                        }
                        node = nextNode;
                    }
                    Require(tree[node].Symbol < 0 && tree[node].Children[0] < 0 && tree[node].Children[1] < 0, "Huffman collision");
                    tree[node].Symbol = i;
                }
            }
            for (auto& degree : degrees) {
                int node = 0;
                while (tree[node].Symbol < 0) {
                    node = tree[node].Children[reader.ReadBit()];
                    Require(node >= 0, "invalid Huffman code");
                }
                unsigned symbol = tree[node].Symbol;
                uint64_t decodedDegree = symbol == 256 ? 255 + reader.ReadGamma() : symbol;
                Require(decodedDegree <= UINT32_MAX, "degree overflow");
                degree = decodedDegree;
            }
            return degrees;
        }
        TBitReader reader{section.Data};
        if (section.Codec == DEGREE_GAMMA) {
            for (auto& degree : degrees) {
                uint64_t decodedValue = reader.ReadGamma() -1;
                Require(decodedValue <= UINT32_MAX, "degree overflow");
                degree = decodedValue;
            }
        } else if (section.Codec == DEGREE_RICE) {
            unsigned riceParameter = reader.ReadBits(8);
            Require(riceParameter < 32, "degree Rice k");
            for (auto& degree : degrees) {
                degree = reader.ReadRice(riceParameter);
            }
        } else {
            throw std::runtime_error("unknown degree codec");
        }
        return degrees;
    }

    TSection SelectDegreeCodec(const TGraph& graph) {
        TSection bestSection = EncodeDegrees(graph.ForwardDegrees, DEGREE_VARINT);
        for (unsigned codec = 1; codec < DEGREE_CODEC_COUNT; codec++) {
            auto section = EncodeDegrees(graph.ForwardDegrees, codec);
            if (section.Data.size() < bestSection.Data.size()) {
                bestSection = std::move(section);
            }
        }
        return bestSection;
    }

} // namespace NGraphCompressor
