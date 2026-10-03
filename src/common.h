#pragma once

#include <bit>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace NGraphCompressor {

    using TBytes = std::vector<uint8_t>;

    inline void Require(bool condition, const char* message) {
        if (!condition) {
            throw std::runtime_error(message);
        }
    }

    inline unsigned GetBitWidth(uint64_t value) {
        return static_cast<unsigned>(std::bit_width(value));
    }

    inline double GetTimeSeconds() {
        return std::chrono::duration<double>(
            std::chrono::steady_clock::now().time_since_epoch()).count();
    }

} // namespace NGraphCompressor
