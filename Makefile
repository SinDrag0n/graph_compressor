CXX ?= g++
CXXFLAGS ?= -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic
SOURCES := $(wildcard src/*.cpp)
HEADERS := $(wildcard src/*.h)

ifeq ($(OS),Windows_NT)
TARGET := run.exe
else
TARGET := run
endif

$(TARGET): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET)

clean:
	rm -f run run.exe

.PHONY: clean
