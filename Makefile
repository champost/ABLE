################################################################################
# Makefile for building ABLE. Edit only if you know what you are doing!
################################################################################

.PHONY: all clean deps

RM := rm -rf

GSLPATH := $(PWD)/deps/gsl
NLOPTPATH := $(PWD)/deps/nlopt

CC = gcc
CXX = g++

INCLUDES := -I deps/gsl/include/ -I deps/nlopt/include/

LDFLAGS := -L deps/gsl/lib/ -L deps/nlopt/lib/
LDFLAGS += -lgsl -lgomp -lgslcblas -lnlopt -lm

CPP_SRCS += \
./main.cpp \
./utils.cpp 

C_SRCS += \
./ms_new.c \
./streec.c 

OBJS += \
./main.o \
./ms_new.o \
./streec.o \
./utils.o 


# Building C++ objects
%.o: ./%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(INCLUDES) -O3 -c -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Building C objects
%.o: ./%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	$(CC) $(INCLUDES) -O3 -c -fopenmp -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# All Target
all: ABLE

# BuildABLE
ABLE: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -o "ABLE" $(OBJS) $(LDFLAGS)
	@echo 'Finished building target: $@'
	@echo ' '

# make deps
deps:
	@echo 'Locally installing GSL in ./deps/gsl'
	$(RM) deps/gsl && mkdir deps/gsl
	tar -xzf deps/gsl-2.1.tar.gz
	cd gsl-2.1 && ./configure --disable-shared --prefix=$(GSLPATH) && make && make install
	@echo ' '
	
	@echo 'Locally installing NLopt in ./deps/nlopt'
	$(RM) deps/nlopt && mkdir deps/nlopt
	tar -xzf deps/nlopt-2.4.2.tar.gz
	cd nlopt-2.4.2 && ./configure --disable-shared --prefix=$(NLOPTPATH) && make && make install
	@echo ' '
	
	$(RM) gsl-2.1 nlopt-2.4.2

# make clean
clean:
	$(RM) $(OBJS)
	@echo ' '
