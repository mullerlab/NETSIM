# ================================================================================
# ================================================================================
# ================================================================================
#
#   Makefile: NETSIM
#
#   Copyright (C) 2018-2020 Lyle Muller
#	http://mullerlab.ca
#
#	COMMAND LINE OPTIONS:
#
# 	EXTERNAL_INPUT=yes		compile in support for external Poisson noise
# 	RELEASE_PROBABILITY=yes		for probabilistic release
#
# ================================================================================
# ================================================================================
# ================================================================================


SRCDIR = src
BUILDDIR = .

#
# compiler settings
#

CC = gcc

#
# complation type
#

FLAGS = -std=gnu11

# debug
ifeq ($(COMPILE_TYPE),debug)
	FLAGS += -ggdb -Wall
endif

# profiler
ifeq ($(COMPILE_TYPE),profiler) 
	FLAGS += -g
endif

# performance
ifeq ($(COMPILE_TYPE),performance)
	FLAGS += -O3 -Wall -march=native -flto
endif

#
# differential compile options
#

ifeq ($(EXTERNAL_INPUT),yes) 
	FLAGS += -DEXTERNAL_INPUT
endif

ifeq ($(RELEASE_PROBABILITY),yes) 
	FLAGS += -DRELEASE_PROBABILITY
endif

#
# library settings
#

ifeq ($(shell uname),Linux)
	LIBS += -lm  
endif

#
# files
#

OBJS = $(BUILDDIR)/netsim.o $(BUILDDIR)/rng.o $(BUILDDIR)/isaac64.o $(BUILDDIR)/gamma_dist_rng.o

TARGET = $(BUILDDIR)/netsim

#
# rules
#

$(TARGET): $(OBJS)
	$(CC) $(FLAGS) $(OBJS) $(LIBS) -o $(TARGET)

$(BUILDDIR)/netsim.o: $(SRCDIR)/netsim.c $(SRCDIR)/helper_functions/helper_functions.h $(SRCDIR)/helper_functions/helper_gaussian.h $(SRCDIR)/rng.c $(SRCDIR)/rng.h
	$(CC) $(FLAGS) $(LIBS) -c $(SRCDIR)/netsim.c -o $(BUILDDIR)/netsim.o

$(BUILDDIR)/rng.o: $(SRCDIR)/rng.c $(SRCDIR)/rng.h $(SRCDIR)/isaac64.c $(SRCDIR)/isaac64.h
	$(CC) $(FLAGS) $(LIBS) -c $(SRCDIR)/rng.c -o $(BUILDDIR)/rng.o

$(BUILDDIR)/isaac64.o: $(SRCDIR)/isaac64.c $(SRCDIR)/isaac64.h
	$(CC) $(FLAGS) $(LIBS) -c $(SRCDIR)/isaac64.c -o $(BUILDDIR)/isaac64.o

$(BUILDDIR)/gamma_dist_rng.o: $(SRCDIR)/gamma_dist_rng.c $(SRCDIR)/gamma_dist_rng.h $(SRCDIR)/rng.c $(SRCDIR)/rng.h
	$(CC) $(FLAGS) $(LIBS) -c $(SRCDIR)/gamma_dist_rng.c -o $(BUILDDIR)/gamma_dist_rng.o

clean:
	rm -f $(TARGET) $(OBJS)
