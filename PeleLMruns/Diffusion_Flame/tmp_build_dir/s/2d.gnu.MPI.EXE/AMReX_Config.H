#ifndef AMREX_HAVE_NO_CONFIG_H
#define AMREX_HAVE_NO_CONFIG_H
#ifndef BL_USE_MPI
#define BL_USE_MPI 1
#endif
#ifndef AMREX_USE_MPI
#define AMREX_USE_MPI 1
#endif
#ifndef USE_SUNDIALS_PP
#define USE_SUNDIALS_PP 1
#endif
#ifndef AMREX_GIT_VERSION
#define AMREX_GIT_VERSION "21.06-53-g14dde9aaf741"
#endif
#ifndef BL_SPACEDIM
#define BL_SPACEDIM 2
#endif
#ifndef AMREX_SPACEDIM
#define AMREX_SPACEDIM 2
#endif
#ifndef BL_FORT_USE_UNDERSCORE
#define BL_FORT_USE_UNDERSCORE 1
#endif
#ifndef AMREX_FORT_USE_UNDERSCORE
#define AMREX_FORT_USE_UNDERSCORE 1
#endif
#ifndef BL_Linux
#define BL_Linux 1
#endif
#ifndef AMREX_Linux
#define AMREX_Linux 1
#endif
#ifndef NDEBUG
#define NDEBUG 1
#endif
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX 1
#endif
#ifndef SIMPLE_TRANSPORT
#define SIMPLE_TRANSPORT 1
#endif
#ifndef BL_USE_VELOCITY
#define BL_USE_VELOCITY 1
#endif
#ifndef BL_NOLINEVALUES
#define BL_NOLINEVALUES 1
#endif
#ifndef AMREX_NOLINEVALUES
#define AMREX_NOLINEVALUES 1
#endif
#ifndef AMREX_PARALLEL_IO
#define AMREX_PARALLEL_IO 1
#endif
#ifndef CG_USE_OLD_CONVERGENCE_CRITERIA
#define CG_USE_OLD_CONVERGENCE_CRITERIA 1
#endif
#ifdef __cplusplus
#endif
#if defined(AMREX_USE_OMP) && !defined(_OPENMP)
#error libamrex was built with OpenMP
#endif
#endif
