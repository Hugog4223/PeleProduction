TOP = /home/hugog/PeleProduction/Submodules
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELE_PHYSICS_HOME ?= ${TOP}/PelePhysics
PELELM_HOME ?= ${TOP}/PeleLM

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = FALSE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

Reactions_dir   := Fuego
Transport_dir   := Simple
Eos_dir         := Fuego
Chemistry_Model := grimech30

USE_SUNDIALS_PP = TRUE


Blocs   := .

F90EXE_sources += Prob_nd.F90 user_defined_fcts_nd.F90 probdata.F90

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
