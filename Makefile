#!/bin/bash

# Makefile for reduced gravity model.

F90=gfortran-mp-6

INCLUDES=-I/opt/local/include
LIBS=-L/opt/local/lib -lnetcdff -lnetcdf
FFLAGS=-cpp -ffree-form -fbounds-check

all:
	make -i clean
	make modules
	make main

clean:
	@rm *.mod
	@rm *.o
	@rm *.out

diagnostics:
	@$(F90) -c src/diagnostics.F90

domain:
	@$(F90) -c src/domain.F90

double:
	@$(F90) -c src/double.F90

forcing:
	@$(F90) -c src/forcing.F90

grid:
	@$(F90) -c src/grid.F90

initial:
	@$(F90) -c src/initial.F90

main:
	$(F90) $(INCLUDES) $(LIBS) $(FFLAGS) src/rgmodel_main.F90 *.o -o rgmodel.out

modules:
	@make double
	@make domain
	@make numerical
	@make physical
	@make variables
	@make overlaps
	@make grid
	@make variables
	@make initial
	@make diagnostics
	@make forcing
	@make params
	@make netcdf
	@make tendency
	@make user
	@make postprocess

netcdf:
	@$(F90) $(INCLUDES) $(LIBS) $(FFLAGS) -g -c src/netcdfinterface.F90

numerical:
	@$(F90) $(FFLAGS) -g -c src/numerical.F90
	
overlaps:
	@$(F90) -c src/overlaps.F90
	
params:
	@$(F90) -c src/params.F90
	
physical:
	@$(F90) -c src/physical.F90
	
postprocess:
	@$(F90) -c src/postprocess.F90

tendency:
	@$(F90) $(FFLAGS) -g -c src/thickness.F90
	@$(F90) $(FFLAGS) -g -c src/uvelocity.F90
	@$(F90) $(FFLAGS) -g -c src/vvelocity.F90

user:
	@$(F90) -c src/user.F90
	
variables:
	@$(F90) -c src/variables.F90
