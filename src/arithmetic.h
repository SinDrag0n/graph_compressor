#pragma once

#include "bitstream.h"

namespace NGraphCompressor {

    constexpr uint64_t ARITHMETIC_HALF = 0x80000000ULL;
    constexpr uint64_t ARITHMETIC_QUARTER = 0x40000000ULL;
    constexpr uint64_t ARITHMETIC_THREE_QUARTERS = 0xc0000000ULL;

    // Inclusive 32-bit interval; multiplication uses uint64_t to avoid overflow.
    class TArithmeticEncoder {
    public:
        void ProcessSymbol(uint32_t cumulativeLow, uint32_t cumulativeHigh, uint32_t total) {
            const uint64_t range = High_ - Low_ + 1;
            High_ = Low_ + range * cumulativeHigh / total - 1;
            Low_ += range * cumulativeLow / total;

            for (;;) {
                if (High_ < ARITHMETIC_HALF) {
                    EmitBit(0);
                } else if (Low_ >= ARITHMETIC_HALF) {
                    EmitBit(1);
                    Low_ -= ARITHMETIC_HALF;
                    High_ -= ARITHMETIC_HALF;
                } else if (Low_ >= ARITHMETIC_QUARTER && High_ < ARITHMETIC_THREE_QUARTERS) {
                    ++PendingBits_;
                    Low_ -= ARITHMETIC_QUARTER;
                    High_ -= ARITHMETIC_QUARTER;
                } else {
                    break;
                }
                Low_ <<= 1;
                High_ = (High_ << 1) | 1;
            }
        }

        TBytes Finish() {
            ++PendingBits_;
            EmitBit(Low_ < ARITHMETIC_QUARTER ? 0 : 1);
            Writer_.WriteBits(0, 32);
            return Writer_.TakeBytes();
        }

    private:
        void EmitBit(unsigned value) {
            Writer_.WriteBit(value);
            while (PendingBits_ != 0) {
                Writer_.WriteBit(1 - value);
                --PendingBits_;
            }
        }

    private:
        TBitWriter Writer_;
        uint64_t Low_ = 0;
        uint64_t High_ = UINT32_MAX;
        uint64_t PendingBits_ = 0;
    };

    class TArithmeticDecoder {
    public:
        explicit TArithmeticDecoder(const TBytes& bytes)
            : Reader_(bytes)
            , Low_(0)
            , High_(UINT32_MAX)
            , Value_(Reader_.ReadBits(32))
        {
        }

        uint32_t GetTarget(uint32_t total) const {
            return static_cast<uint32_t>(((Value_ - Low_ + 1) * total - 1) / (High_ - Low_ + 1));
        }

        void ProcessSymbol(uint32_t cumulativeLow, uint32_t cumulativeHigh, uint32_t total) {
            const uint64_t range = High_ - Low_ + 1;
            High_ = Low_ + range * cumulativeHigh / total - 1;
            Low_ += range * cumulativeLow / total;

            for (;;) {
                if (High_ < ARITHMETIC_HALF) {
                    // The interval is in the lower half; no offset is subtracted.
                } else if (Low_ >= ARITHMETIC_HALF) {
                    Low_ -= ARITHMETIC_HALF;
                    High_ -= ARITHMETIC_HALF;
                    Value_ -= ARITHMETIC_HALF;
                } else if (Low_ >= ARITHMETIC_QUARTER && High_ < ARITHMETIC_THREE_QUARTERS) {
                    Low_ -= ARITHMETIC_QUARTER;
                    High_ -= ARITHMETIC_QUARTER;
                    Value_ -= ARITHMETIC_QUARTER;
                } else {
                    break;
                }
                Low_ <<= 1;
                High_ = (High_ << 1) | 1;
                Value_ = (Value_ << 1) | Reader_.ReadBit();
            }
        }

    private:
        TBitReader Reader_;
        uint64_t Low_;
        uint64_t High_;
        uint64_t Value_;
    };

} // namespace NGraphCompressor
