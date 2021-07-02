#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#define NUM_ZONES 6

using namespace amrex::literals;

struct ProbParm
   : amrex::Gpu::Managed
{
   amrex::Real P_mean = 101325.0;      /* change */
   amrex::Real fuel_ox_split = .0009;
   amrex::Real ox_air_split =  0.01;
   amrex::Real Vel_zero_th = 0.0002;
   amrex::Real pipeTh = .000;
   amrex::Real pipeBL = .001;
   amrex::Real blobx = .001;                /* initial spark positionx */
   amrex::Real bloby = .00127;
   amrex::Real blobr = .001;
   amrex::Real blobw = .000001;            /* initial spark withd */
   amrex::Real blobT = 1200;              /* initial spark Temperature */  
   amrex::Real V_fu = 0.30;     /* change */
   amrex::Real V_ox = 0.30;     /* change */  
   amrex::Real V_air = 0;
   amrex::Real T_fu = 300.0;     /* change */
   amrex::Real T_ox = 300.0;     /* change */
   amrex::Real T_air = 300.0;     /* change */
   
   const int BL_FUELPIPE = 1;
   const int BL_OUTFLOW  = 2;
   const int BL_OXIDIZER = 3;
   const int BL_AIR      = 4;
   const int BL_PIPEEND  = 5;
   const int BL_VOLUME   = 6;
   const int BL_NZONES   = 6;

   const amrex::Real BL_XLO = 0;
   const amrex::Real BL_YLO = 1;
   const amrex::Real BL_XHI = 2;
   const amrex::Real BL_YHI = 3;
   const amrex::Real BL_ZLO = 4;
   const amrex::Real BL_ZHI = 5;
   const amrex::Real BL_INTERIOR = 6;

   amrex::Array<amrex::Array<amrex::Real, NUM_ZONES>, NUM_SPECIES > Y_bc;

   amrex::Array<amrex::Real,NUM_ZONES> u_bc;
   amrex::Array<amrex::Real,NUM_ZONES> v_bc;
   amrex::Array<amrex::Real,NUM_ZONES> w_bc;
   amrex::Array<amrex::Real,NUM_ZONES> T_bc;
   amrex::Array<amrex::Real,NUM_ZONES> h_bc;
   amrex::Array<amrex::Real,NUM_ZONES> rho_bc;

   amrex::Real bathID = -1;
   amrex::Real fuelID = -1;
   amrex::Real oxidID = -1;
};
#endif
