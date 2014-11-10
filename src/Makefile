# Variable definitions:
PROGRAM_NAME = corrlnfields
COMP = g++
BIN = ../bin
OBJ = ../obj
SRC = ../src

#FFTWFLAGS = -lfftw3
GSLFLAGS = -lgsl -lgslcblas

# General instructions:
all: $(BIN)/$(PROGRAM_NAME) $(BIN)/gencovl

clean:
	rm -f $(BIN)/*
	rm -f $(OBJ)/*


# Instructions for compiling main program:

$(BIN)/$(PROGRAM_NAME): $(SRC)/corrlnfields.cpp $(OBJ)/ParameterList.o $(OBJ)/Utilities.o $(OBJ)/gsl_aux.o 
	$(COMP) $(SRC)/corrlnfields.cpp $(OBJ)/ParameterList.o $(OBJ)/Utilities.o $(OBJ)/gsl_aux.o -o $@ $(GSLFLAGS)

$(OBJ)/ParameterList.o: $(SRC)/ParameterList.cpp $(SRC)/ParameterList.hpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/ParameterList.cpp -o $@

$(OBJ)/gsl_aux.o: $(SRC)/gsl_aux.cpp $(SRC)/gsl_aux.hpp
	$(COMP) -c $(SRC)/gsl_aux.cpp -o $@

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/Utilities.cpp -o $@

$(OBJ)/interpol.o: $(SRC)/interpol.c $(SRC)/interpol.h
	$(COMP) -c $(SRC)/interpol.c -o $@


# Instructions for auxiliary programs:

$(BIN)/gencovl: $(SRC)/gencovl.cpp $(OBJ)/Utilities.o $(OBJ)/interpol.o $(OBJ)/ParameterList.o
	$(COMP) $(SRC)/gencovl.cpp $(OBJ)/Utilities.o $(OBJ)/interpol.o $(OBJ)/ParameterList.o -o $@
