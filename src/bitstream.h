#pragma once

#include "common.h"

#include <utility>

namespace NGraphCompressor {

    // Packed bits are written MSB-first within each byte. Only the last byte is padded.
    class TBitWriter {
    public:
        void WriteBit(unsigned value) {
            if ((BitCount_ & 7) == 0) {
                Bytes_.push_back(0);
            }
            Bytes_.back() |= (value & 1) << (7 - (BitCount_ & 7));
            ++BitCount_;
        }

        void WriteBits(uint64_t value, unsigned bitCount) {
            Require(bitCount <= 64, "bit width");
            for (unsigned remaining = bitCount; remaining > 0; --remaining) {
                WriteBit(static_cast<unsigned>(value >> (remaining - 1)));
            }
        }

        // Encode value in [0, range) using floor(log2(range)) or one extra bit.
        void WriteBounded(uint64_t value, uint64_t range) {
            Require(range != 0 && value < range, "bounded value");
            const unsigned shortWidth = GetBitWidth(range) -1;
            // For width=63, unsigned wraparound represents 2^64 - range safely.
            const uint64_t shortCodeCount = shortWidth == 63
                ? uint64_t(0) - range
                : (uint64_t(1) << (shortWidth + 1)) - range;

            if (value < shortCodeCount) {
                WriteBits(value, shortWidth);
            } else {
                WriteBits(value + shortCodeCount, shortWidth + 1);
            }
        }

        void WriteUnary(uint64_t value) {
            for (uint64_t i = 0; i < value; ++i) {
                WriteBit(0);
            }
            WriteBit(1);
        }

        void WriteRice(uint64_t value, unsigned parameter) {
            WriteUnary(value >> parameter);
            WriteBits(value, parameter);
        }

        void WriteGamma(uint64_t value) {
            Require(value != 0, "gamma zero");
            const unsigned suffixWidth = GetBitWidth(value) -1;
            WriteUnary(suffixWidth);
            WriteBits(value, suffixWidth);
        }

        void WriteDelta(uint64_t value) {
            const unsigned bitCount = GetBitWidth(value);
            WriteGamma(bitCount);
            WriteBits(value, bitCount - 1);
        }

        const TBytes& GetBytes() const {
            return Bytes_;
        }

        uint64_t GetBitCount() const {
            return BitCount_;
        }

        TBytes TakeBytes() {
            return std::move(Bytes_);
        }

    private:
        TBytes Bytes_ = {};
        uint64_t BitCount_ = 0;
    };

    class TBitReader {
    public:
        explicit TBitReader(const TBytes& bytes)
            : Bytes_(bytes)
            , BitPosition_(0)
        {
        }

        unsigned ReadBit() {
            Require(BitPosition_ / 8 < Bytes_.size(), "truncated bitstream");
            const unsigned value = (Bytes_[BitPosition_ / 8] >> (7 - (BitPosition_ & 7))) & 1;
            ++BitPosition_;
            return value;
        }

        uint64_t ReadBits(unsigned bitCount) {
            Require(bitCount <= 64, "bit width");
            uint64_t value = 0;
            for (unsigned i = 0; i < bitCount; ++i) {
                value = (value << 1) | ReadBit();
            }
            return value;
        }

        uint64_t ReadBounded(uint64_t range) {
            Require(range != 0, "bounded zero");
            const unsigned shortWidth = GetBitWidth(range) -1;
            const uint64_t shortCodeCount = shortWidth == 63
                ? uint64_t(0) - range
                : (uint64_t(1) << (shortWidth + 1)) - range;
            uint64_t value = ReadBits(shortWidth);
            if (value >= shortCodeCount) {
                value = ((value << 1) | ReadBit()) - shortCodeCount;
            }
            Require(value < range, "bounded corrupt");
            return value;
        }

        uint64_t ReadUnary(uint64_t limit = uint64_t(1) << 32) {
            uint64_t value = 0;
            while (ReadBit() == 0) {
                Require(value < limit, "unary overflow");
                ++value;
            }
            return value;
        }

        uint64_t ReadRice(unsigned parameter) {
            const uint64_t quotient = ReadUnary();
            Require(parameter < 32 && quotient <= (UINT32_MAX >> parameter), "rice overflow");
            return (quotient << parameter) | ReadBits(parameter);
        }

        uint64_t ReadGamma() {
            const unsigned suffixWidth = static_cast<unsigned>(ReadUnary(63));
            return (uint64_t(1) << suffixWidth) | ReadBits(suffixWidth);
        }

    private:
        const TBytes& Bytes_;
        uint64_t BitPosition_;
    };

} // namespace NGraphCompressor
