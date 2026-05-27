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
#	MULTILAYER=yes			for networks with multiple layers
# 	EXTERNAL_INPUT=yes		compile in support for external Poisson noise
#       VARIABLE_EXTERNAL_INPUT=yes	for per neuron and time external Poisson noise
# 	RELEASE_PROBABILITY=yes		for probabilistic release
#	CURRENT_INJECTION=yes		external file for current injection term I_e
#       SPARSE_CURRENT_INJECTION=yes	external sparse file for current injection term I_e
#   SINGLE_NEURON_PERTURBATION=yes  DS3
#
# ================================================================================
# ================================================================================
# ================================================================================


SRCDIR = src
BUILDDIR = ./obj

#
# compiler settings
#

CC = clang

#
# complation type
#

CFLAGS = -std=c11 -D_POSIX_C_SOURCE=200809L -D_XOPEN_SOURCE=500
LDFLAGS = 

# debug
ifeq ($(COMPILE_TYPE),debug)
    $(info Using debug flags)
    CFLAGS += -g -ggdb -gdwarf-4 -Wall -pedantic
endif

# profiler
ifeq ($(COMPILE_TYPE),profiler)
    $(info Using profiler flags)
    CFLAGS += -O3 -Wall -march=native -flto -g -ggdb -gdwarf-4
endif

# 
ifeq ($(COMPILE_TYPE),profilegen)
    $(info Using profile generate flags)
    CFLAGS += -O3 -march=native -flto -fprofile-generate -fprofile-arcs -ftest-coverage
    LDFLAGS += -fprofile-arcs
endif

ifeq ($(COMPILE_TYPE),profileuse)
    $(info Using profile use flags)
    CFLAGS += -O3 -march=native -flto -fprofile-use
endif

# performance
ifeq ($(COMPILE_TYPE),performance)
    $(info Using performance flags)
    CFLAGS += -O3 -flto
    ifeq ($(shell uname), Darwin)
        ifeq ($(shell uname -m), arm64)
            CFLAGS += -mcpu=apple-m1
        else
            CFLAGS += -march=native
        endif
    else
        CFLAGS += -march=native
    endif
endif



#
# differential compile options
#

ifeq ($(EXTERNAL_INPUT),yes)
    $(info Using external input)
    CFLAGS += -DEXTERNAL_INPUT
endif

ifeq ($(VARIABLE_EXTERNAL_INPUT),yes)
    $(info Using variable external input)
    CFLAGS += -DVARIABLE_EXTERNAL_INPUT
endif

ifeq ($(RELEASE_PROBABILITY),yes)
    $(info Using release probability)
    CFLAGS += -DRELEASE_PROBABILITY
endif

ifeq ($(MULTILAYER),yes)
    $(info Using Multilayer connectivity)
    CFLAGS += -DMULTILAYER
endif

ifeq ($(CURRENT_INJECTION),yes)
    $(info Using current injection)
    CFLAGS += -DCURRENT_INJECTION
endif

ifeq ($(SPARSE_CURRENT_INJECTION),yes)
    $(info Using sparse current injection)
    CFLAGS += -DSPARSE_CURRENT_INJECTION
endif

ifeq ($(FAST_BINOMIAL),yes)
    $(info Using fast binomial)
    CFLAGS += -DFAST_BINOMIAL
endif

ifeq ($(LONG_SIMULATION),yes)
    $(info Write out small binary files)
    CFLAGS += -DLONG_SIMULATION
endif

ifeq ($(SINGLE_NEURON_PERTURBATION),yes)
    $(info Using single neuron perturbation code (set neuron to threshold))
    CFLAGS += -DSINGLE_NEURON_PERTURBATION
endif

ifeq ($(USE_CONNECTOR_RNG),yes)
    $(info Using separate connector RNG)
    CFLAGS += -DUSE_SEPARATE_CONNECTOR_RNG
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

SRC = $(wildcard $(SRCDIR)/*.c)
INC = $(wildcard $(SRCDIR)/*.h)

OBJS = $(patsubst %.c, %.o, $(filter %.c, $(subst $(SRCDIR), $(BUILDDIR), $(SRC))))

TARGET = netsim

#
# rules
#

.PHONY: build
build: $(BUILDDIR) $(OBJS) $(INC)
	@$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)
	@echo [LD] Linked $(OBJS) into $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c $(INC)
	@$(CC) $(CFLAGS) -c $< -o $@
	@echo [CC] Compiled $< into $@

.PHONY: clean
clean:
	@rm -f $(OBJS) $(TARGET)
	@rm -d -f $(BUILDDIR)
	@echo Cleaned $(OBJS), $(BUILDDIR), and $(TARGET)

.PHONY: rebuild
rebuild: clean build

$(BUILDDIR):
	mkdir -p $@
