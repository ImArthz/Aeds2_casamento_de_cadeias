CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -Iinclude
TARGET = bin/casamento_cadeias
SRCS = src/main.cpp src/casamento_cadeias.cpp src/fasta_io.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	@mkdir -p bin benchmarks
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -rf bin/*.exe bin/casamento_cadeias
