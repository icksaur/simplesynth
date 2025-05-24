# Makefile for SDL Audio project
# Based on the Visual Studio project configuration

# Compiler and flags
CXX = g++
CC = gcc
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -g
CFLAGS = -Wall -Wextra -O2 -g

# Libraries and include paths
SDL2_CFLAGS = $(shell pkg-config --cflags sdl2)
SDL2_LIBS = $(shell pkg-config --libs sdl2)
OPENGL_LIBS = -lGL -lm

# Combine all flags
INCLUDES = $(SDL2_CFLAGS)
LIBS = $(SDL2_LIBS) $(OPENGL_LIBS)

# Source files
CPP_SOURCES = BitFont.cpp GLObject.cpp gui.cpp main.cpp SynthNode.cpp synth.cpp
C_SOURCES = glxw.c

# Object files
CPP_OBJECTS = $(CPP_SOURCES:.cpp=.o)
C_OBJECTS = $(C_SOURCES:.c=.o)
OBJECTS = $(CPP_OBJECTS) $(C_OBJECTS)

# Target executable
TARGET = sdlaudio

# Default target
all: $(TARGET)

# Build the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LIBS)

# Compile C++ source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile C source files  
%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Check if SDL2 is available
check-deps:
	@echo "Checking dependencies..."
	@pkg-config --exists sdl2 && echo "✓ SDL2 found" || echo "✗ SDL2 not found"
	@echo "OpenGL should be available through mesa or drivers"

# Print help
help:
	@echo "Available targets:"
	@echo "  all       - Build the project (default)"
	@echo "  clean     - Remove build artifacts"
	@echo "  check-deps- Check if dependencies are available"
	@echo "  help      - Show this help message"

.PHONY: all clean check-deps help
