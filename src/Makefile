# Variable definitions:
PROGRAM_NAME = corrlnfields
COMP = g++
BIN = ../bin
OBJ = ../obj
SRC = ../src

#FFTWFLAGS = -lfftw3
GSLFLAGS = -lgsl -lgslcblas

# General instructions:
all: $(BIN)/$(PROGRAM_NAME)

clean:
	rm -f $(BIN)/*
	rm -f $(OBJ)/*


# Instructions for compiling files:

$(BIN)/$(PROGRAM_NAME): $(SRC)/main.cpp $(OBJ)/ParameterList.o $(OBJ)/Utilities.o $(OBJ)/gsl_aux.o 
	$(COMP) $(SRC)/main.cpp $(OBJ)/ParameterList.o $(OBJ)/Utilities.o $(OBJ)/gsl_aux.o -o $@ $(GSLFLAGS)

$(OBJ)/ParameterList.o: $(SRC)/ParameterList.cpp $(SRC)/ParameterList.hpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/ParameterList.cpp -o $@

$(OBJ)/gsl_aux.o: $(SRC)/gsl_aux.cpp $(SRC)/gsl_aux.hpp
	$(COMP) -c $(SRC)/gsl_aux.cpp -o $@

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/Utilities.cpp -o $@

