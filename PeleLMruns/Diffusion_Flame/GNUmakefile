SUBMODS ?= /Users/mday2/src/CCSE/Pele/ALT/HUGO/PeleProduction/Submodules
AMREX_HOME         ?= ${SUBMODS}/amrex
IAMR_HOME          ?= ${SUBMODS}/IAMR
PELE_PHYSICS_HOME  ?= ${SUBMODS}/PelePhysics
PELELM_HOME        ?= ${SUBMODS}/PeleLM
AMREX_HYDRO_HOME   ?= ${SUBMODS}/AMReX-Hydro

DIM             = 2
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

DEBUG           = TRUE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

USE_SUNDIALS_PP = TRUE
USE_KLU_PP      = FALSE

Chemistry_Model = drm19
Eos_Model       = Fuego
Reactor_dir     = cvode
Transport_Model = Simple

Blocs   := .

include $(PELELM_HOME)/Tools/Make/Make.PeleLM

ifeq ($(USE_CUDA),TRUE)
CXXFLAGS += -Xptxas --disable-optimizer-constants
endif
