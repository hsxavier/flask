# Variable definitions:
PROGRAM_NAME = corrlnfields
COMP = g++
BIN = ../bin
OBJ = ../obj
SRC = ../src

#FFTWFLAGS = -lfftw3

# General instructions:
all: $(BIN)/$(PROGRAM_NAME)

clean:
	rm -f $(BIN)/*
	rm -f $(OBJ)/*


# Instructions for compiling files:

$(BIN)/$(PROGRAM_NAME): $(SRC)/main.cpp $(OBJ)/ParameterList.o $(OBJ)/Utilities.o 
	$(COMP) $(SRC)/main.cpp $(OBJ)/ParameterList.o $(OBJ)/Utilities.o -o $@

$(OBJ)/ParameterList.o: $(SRC)/ParameterList.cpp $(SRC)/ParameterList.hpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/ParameterList.cpp -o $@

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/Utilities.cpp -o $@
