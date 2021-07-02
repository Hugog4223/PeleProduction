TOP = /home/hugog/PeleProduction/Submodules
AMREX_HOME ?= ${TOP}/amrex
IAMR_HOME ?= ${TOP}/IAMR
PELE_PHYSICS_HOME ?= ${TOP}/PelePhysics
PELELM_HOME ?= ${TOP}/PeleLM
AMREX_HYDRO_HOME	 = ${TOP}/AMReX-Hydro

DIM             = 2
COMP            = gnu
DEBUG           = FALSE
USE_MPI         = TRUE
USE_OMP         = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE

Reactions_dir   := Fuego
Transport_dir   := Simple
Eos_dir         := Fuego
Chemistry_Model := drm19

# CVODE
USE_SUNDIALS_PP = TRUE
USE_KLU_PP      = FALSE

Blocs   := .

CEXE_sources +=
F90EXE_sources +=
CEXE_headers += 
FEXE_headers += 

include $(PELELM_HOME)/Tools/Make/Make.PeleLM
