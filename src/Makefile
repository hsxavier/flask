# Variable definitions:
PROGRAM_NAME = corrlnfields
COMP = g++
BIN = ../bin
OBJ = ../obj
SRC = ../src

CXXFLAGS = -I/home/skems/prog/Healpix_3.11/src/cxx/generic_gcc/include
LDFLAGS = -L/home/skems/prog/Healpix_3.11/src/cxx/generic_gcc/lib
LIBS = -lgsl -lgslcblas -lhealpix_cxx -lcxxsupport -lcfitsio -lsharp -lfftpack -lc_utils -lgomp

# General instructions:
all: $(BIN)/$(PROGRAM_NAME)

clean:
	rm -f $(BIN)/*
	rm -f $(OBJ)/*

$(BIN)/$(PROGRAM_NAME): $(OBJ)/corrlnfields.o $(OBJ)/ParameterList.o $(OBJ)/Utilities.o \
			$(OBJ)/gsl_aux.o $(OBJ)/s2kit10_naive.o $(OBJ)/interpol.o \
			$(OBJ)/corrlnfields_aux.o $(OBJ)/Cosmology.o $(OBJ)/Integral.o
	$(COMP) $(LDFLAGS) $(OBJ)/corrlnfields.o $(OBJ)/ParameterList.o $(OBJ)/Utilities.o \
			$(OBJ)/gsl_aux.o $(OBJ)/s2kit10_naive.o $(OBJ)/interpol.o \
			$(OBJ)/corrlnfields_aux.o $(OBJ)/Cosmology.o $(OBJ)/Integral.o -o $@ $(LIBS)

$(OBJ)/corrlnfields.o: $(SRC)/corrlnfields.cpp 
	$(COMP) $(CXXFLAGS) -c $(SRC)/corrlnfields.cpp -o $@

$(OBJ)/corrlnfields_aux.o: $(SRC)/corrlnfields_aux.cpp $(SRC)/corrlnfields_aux.hpp $(SRC)/ParameterList.hpp \
	$(SRC)/gsl_aux.hpp 
	$(COMP) -c $(SRC)/corrlnfields_aux.cpp -o $@

$(OBJ)/Cosmology.o: $(SRC)/Cosmology.cpp $(SRC)/Cosmology.hpp $(SRC)/Utilities.hpp $(SRC)/Integral.hpp \
	$(SRC)/ParameterList.hpp $(SRC)/interpol.h
	$(COMP) -c $(SRC)/Cosmology.cpp -o $@

$(OBJ)/Integral.o: $(SRC)/Integral.cpp $(SRC)/Integral.hpp $(SRC)/Utilities.hpp $(SRC)/Cosmology.hpp
	$(COMP) -c $(SRC)/Integral.cpp -o $@

$(OBJ)/ParameterList.o: $(SRC)/ParameterList.cpp $(SRC)/ParameterList.hpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/ParameterList.cpp -o $@

$(OBJ)/gsl_aux.o: $(SRC)/gsl_aux.cpp $(SRC)/gsl_aux.hpp
	$(COMP) -c $(SRC)/gsl_aux.cpp -o $@

$(OBJ)/Utilities.o: $(SRC)/Utilities.cpp $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/Utilities.cpp -o $@

$(OBJ)/interpol.o: $(SRC)/interpol.c $(SRC)/interpol.h
	$(COMP) -c $(SRC)/interpol.c -o $@

$(OBJ)/s2kit10_naive.o: $(SRC)/s2kit10_naive.cpp $(SRC)/s2kit10_naive.hpp $(SRC)/interpol.h $(SRC)/Utilities.hpp
	$(COMP) -c $(SRC)/s2kit10_naive.cpp -o $@

