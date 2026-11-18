# Self-contained build for the boolean-result consistency checker.
# No external test framework is required; only the Clipper2 library sources.
CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
ROOT      := ..
INC       := -I$(ROOT)/Clipper2Lib/include
LIB_SRC   := $(ROOT)/Clipper2Lib/src/clipper.engine.cpp \
             $(ROOT)/Clipper2Lib/src/clipper.offset.cpp \
             $(ROOT)/Clipper2Lib/src/clipper.rectclip.cpp \
             $(ROOT)/Clipper2Lib/src/clipper.triangulation.cpp

BIN  := consistency_checker
DBG  := diagnose
SRC  := main.cpp checker.cpp $(LIB_SRC)

.PHONY: all clean run diagnose
all: $(BIN)

$(BIN): $(SRC) checker.h pip_reference.h
	$(CXX) $(CXXFLAGS) $(INC) $(SRC) -o $(BIN)

# Standalone diagnostic: characterises where mismatches sit relative to
# input/solution edges (used to confirm boundary-band vs deep classification).
$(DBG): diagnose.cpp checker.cpp checker.h pip_reference.h $(LIB_SRC)
	$(CXX) $(CXXFLAGS) $(INC) diagnose.cpp checker.cpp $(LIB_SRC) -o $(DBG)

diagnose: $(DBG)

run: $(BIN)
	./$(BIN)

clean:
	rm -f $(BIN) $(DBG)

