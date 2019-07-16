EXE  = experiment
MPI_EXE = mpi_experiment
SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include
LDFLAGS += -lm

SRC = $(wildcard $(SRC_DIR)/*.cc)
OBJ = $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)
INC = -I $(INC_DIR)/

CC      = g++
CFLAGS  = -Wall -std=c++14 -g

COCOFLAGS = -g -ggdb -std=c89 -pedantic -Wall -Wextra -Wstrict-prototypes -Wshadow -Wno-sign-compare -Wconversion

.PHONY: all
all: $(EXE)

.PHONY: mpi
mpi: $(MPI_EXE)

$(EXE): $(OBJ) obj/coco.o obj/experiment.o
	${CC} ${CFLAGS} -o $(EXE) $(OBJ) obj/coco.o obj/experiment.o ${LDFLAGS}

$(MPI_EXE): $(OBJ) obj/coco.o obj/mpi.o
	mpiCC ${CFLAGS} -o $(MPI_EXE) $(OBJ) obj/coco.o obj/mpi.o $(LDFLAGS)
	
obj/mpi.o: include/coco.h src/coco.c src/mpi.c
	mpiCC -c $(CFLAGS) $(INC) -o obj/mpi.o src/mpi.c

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc $(INC_DIR)/*
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

obj/coco.o: include/coco.h src/coco.c
	gcc -c $(COCOFLAGS) $(INC) -o obj/coco.o src/coco.c
obj/experiment.o: include/coco.h src/coco.c src/experiment.c
	g++ -c $(CFLAGS) $(INC) -o obj/experiment.o src/experiment.c

.PHONY:  clean
clean:
	rm -f $(OBJ_DIR)/*.o $(EXE) $(MPI_EXE)

