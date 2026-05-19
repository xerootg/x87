# Build the x87 test harness on Linux.
#
#   make           - builds ./x87test
#   make run       - builds and runs the test
#   make clean
#
# On macOS, change FORMAT to macho64.

CXX      ?= g++
NASM     ?= nasm
CXXFLAGS ?= -std=c++20 -O2 -g -Wall -Wno-unused-but-set-variable -Wno-sign-compare
ASFLAGS  ?= -f elf64

TEST_BIN := x87test
TEST_OBJ := x87test.o
ASM_OBJ  := x87testasm.o

all: $(TEST_BIN)

$(TEST_BIN): $(TEST_OBJ) $(ASM_OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm

$(TEST_OBJ): x87test.cpp ../x87fp64.h ../x87fp80.h ../x87common.h ../x87fpext.h \
             ../x87fp64trans.cpp ../x87fp80.cpp ../x87fp80trans.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ASM_OBJ): x87testasm.asm
	$(NASM) $(ASFLAGS) -o $@ $<

run: $(TEST_BIN)
	./$(TEST_BIN)

clean:
	rm -f $(TEST_BIN) $(TEST_OBJ) $(ASM_OBJ)

.PHONY: all run clean
