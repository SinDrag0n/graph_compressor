#pragma once

#include "common.h"

#include <fstream>
#include <string>

namespace NGraphCompressor {

    inline void AppendInteger(TBytes& bytes, uint64_t value, unsigned byteCount) {
        for (unsigned i = 0; i < byteCount; ++i) {
            bytes.push_back(value & 255);
            value >>= 8;
        }
    }

    inline uint64_t ReadInteger(const TBytes& bytes, size_t& position, unsigned byteCount) {
        Require(byteCount <= 8 && position <= bytes.size() && byteCount <= bytes.size() - position,
                "truncated integer");
        uint64_t value = 0;
        for (unsigned i = 0; i < byteCount; ++i) {
            value |= uint64_t(bytes[position++]) << (8 * i);
        }
        return value;
    }

    inline void WriteVarint(TBytes& bytes, uint64_t value) {
        while (value >= 128) {
            bytes.push_back((value & 127) | 128);
            value >>= 7;
        }
        bytes.push_back(value);
    }

    inline uint64_t ReadVarint(const TBytes& bytes, size_t& position) {
        uint64_t value = 0;
        for (unsigned i = 0; i < 10; ++i) {
            const uint64_t nextByte = ReadInteger(bytes, position, 1);
            Require(i < 9 || nextByte < 2, "varint overflow");
            value |= (nextByte & 127) << (7 * i);
            if ((nextByte & 128) == 0) {
                return value;
            }
        }
        throw std::runtime_error("varint overflow");
    }

    inline TBytes ReadFile(const std::string& path) {
        std::ifstream input(path, std::ios::binary);
        Require(bool(input), "cannot open input");
        input.seekg(0, std::ios::end);
        const auto byteCount = input.tellg();
        Require(byteCount >= 0, "input size");
        input.seekg(0);
        TBytes bytes(static_cast<size_t>(byteCount));
        input.read(reinterpret_cast<char*>(bytes.data()), byteCount);
        Require(bool(input), "input read failed");
        return bytes;
    }

    inline void WriteFile(const std::string& path, const TBytes& bytes) {
        std::ofstream output(path, std::ios::binary);
        Require(bool(output), "cannot open output");
        output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
        Require(bool(output), "output write failed");
    }

} // namespace NGraphCompressor
