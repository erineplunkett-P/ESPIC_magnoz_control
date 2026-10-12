CC=g++
CFLAGS=-std=c++11 -O3 -fopenmp -pedantic -Wall# -DNOPOISSON # -g -DDEBUG
LINEAR_SOLVER_INCLUDE_DIR=./eigen-3.3.9

BIN=../bin/

#PROG=espic_debug
PROG=espic
EXEC=$(PROG).x
INCLUDES=-I$(LINEAR_SOLVER_INCLUDE_DIR)

OBJS = es_PIC.o pic_part.o inject_flow_inout.o pic_io.o \
			 solve_field.o setup_matrix_rhs.o update_matrix_rhs.o \
			 eval_hist_vel_field.o eval_azimuthal.o objects_in_domain.o beams_in_domain.o \
			 str_split.o parse.o source.o sample_field.o vdf.o edf.o comp_ellint.o

default : esPIC

all : esPIC gen_ambient

esPIC : $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(EXEC) -lm -lstdc++
	/bin/mv $(EXEC) $(BIN)

gen_ambient : gen_ambient.o
	$(CC) $(CFLAGS) gen_ambient.o -o gen_ambient.x -lm -lstdc++
	/bin/mv gen_ambient.x $(BIN)

.cpp.o :
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean :
	/bin/rm -f *.o *.x .*.un~

run :
	./$(EXEC)
