#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <pelelm_prob_parm.H>
#include <pmf.H>
#include <pmf_data.H>
#include <PeleLM_parm.H>

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <EOS.H>
#include <Mechanism.H>
/* #include <PelePhysics.H> */
#include <cmath>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setupbc(amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 PmfData const *pmf_data) 
{
    amrex::Real fuelZone;
    amrex::Real airZone;
    amrex::Real oxZone;
    amrex::Real volZone;
    amrex::Real peZone;
    amrex::Real ofZone;
    amrex::Real region;

    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0}; 
   
    amrex::Real iN2  = -1;
    amrex::Real iO2  = -1;
    amrex::Real iH2  = -1;
    amrex::Real iCH4 = -1;
    amrex::Real iCO2 = -1;
    amrex::Real iH2O = -1;

    amrex::Vector<std::string> spec_names;
    CKSYMS_STR(spec_names);

    for (int k = 0; k < NUM_SPECIES; k++){
      if (spec_names[k]="O2") {iO2=k;
      } else if (spec_names[k]="CH4") {iCH4=k;
      } else if (spec_names[k]="N2") {iN2=k;
      } else if (spec_names[k]="H2O") {iH2O=k;
      } else if (spec_names[k]="CO2") {iCO2=k;}
    }

    /*A diffusion flame*/

    fuelZone = prob_parm.BL_FUELPIPE;
    oxZone   = prob_parm.BL_OXIDIZER;
    airZone  = prob_parm.BL_AIR;
    volZone  = prob_parm.BL_VOLUME;
    peZone   = prob_parm.BL_PIPEEND;
    ofZone   = prob_parm.BL_OUTFLOW;
    
    /*Fuel*/

    for (int n = 0; n < NUM_SPECIES; n++){
      massfrac[n] = 0.0;
    }

    massfrac[(int) iCH4] = 1.0;
    massfrac[(int) iO2] = 0.0;
    massfrac[(int) iH2O] = 0.0;
    massfrac[(int) iN2] = 1.0 - massfrac[(int) iCH4] - massfrac[(int) iO2] - massfrac[(int) iH2O];

    for (int n = 0; n < NUM_SPECIES; n++){
      Y_bc[fuelZone][n-1] = 0.0;
    }
    T_bc[fuelZone] = prob_parm.T_fu;
    u_bc[(int) fuelZone] = 0.0;
    v_bc[(int) fuelZone] = prob_parm.V_fu;
    w_bc[(int) fuelZone] = 0.0;

    /*Oxidizer*/
    for (int n = 0; n < NUM_SPECIES; n++){
      massfrac[n] = 0.0;
    }
    massfrac[(int) iO2] = 0.2395;
    massfrac[(int) iCO2] = 0.000;
    massfrac[(int) iN2] = 1.0 - massfrac[(int) iO2] - massfrac[(int) iCO2];

    for (int n = 0; n < NUM_SPECIES; n++){
      Y_bc[(int) oxZone][n-1] = massfrac[n];
    }
    T_bc[(int) oxZone] = prob_parm.T_ox;
    u_bc[(int) oxZone] = 0.0;
    v_bc[(int) oxZone] = prob_parm.V_ox;
    w_bc[(int) oxZone] = 0.0;

    /*air*/
    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = 0.0;
    }
    molefrac[(int) iO2] = 0.21;
    molefrac[(int) iN2] = 1.0 - molefrac[(int) iO2];

    EOS::X2Y(molefrac, massfrac);
    for (int n = 0; n < NUM_SPECIES; n++){
      Y_bc[(int) airZone][n-1] = molefrac[n];
    }

    T_bc[(int) airZone] = prob_parm.T_air;
    u_bc[(int) airZone] = 0.0;
    v_bc[(int) airZone] = prob_parm.V_air;
    w_bc[(int) airZone] = 0.0;

    /*pipeend*/
    for (int n = 0; n < NUM_SPECIES; n++){
      Y_bc[(int) peZone][n-1] = Y_bc[(int) oxZone][n-1];
    }
    T_bc[(int) peZone] = T_bc[(int) oxZone];
    u_bc[(int) peZone] = u_bc[(int) oxZone];
    v_bc[(int) peZone] = 0.0;
    w_bc[(int) peZone] = 0.0;

    /*volume(as air)*/
    for (int n = 0; n < NUM_SPECIES; n++){
      Y_bc[(int) volZone][n-1] = Y_bc[(int) airZone][n-1];
    }
    T_bc[(int) volZone] = T_bc[(int) airZone];
    u_bc[(int) volZone] = u_bc[(int) airZone];
    v_bc[(int) volZone] = v_bc[(int) airZone];
    w_bc[(int) volZone] = 0.0;

    /*outflow(as air)*/
    for (int n = 0; n < NUM_SPECIES; n++){
      Y_bc[(int) ofZone][n-1] = Y_bc[(int) airZone][n-1];
    }
    T_bc[(int) ofZone] = T_bc[(int) airZone];
    u_bc[(int) ofZone] = u_bc[(int) airZone];
    v_bc[(int) ofZone] = v_bc[(int) airZone];
    w_bc[(int) ofZone] = 0.0;

    for (int zone = 0; zone < prob_parm.BL_NZONES; zone++){
    EOS::PYT2R(prob_parm.Pamb, &(Y_bc[zone][-1]), T_bc[zone], rho_bc[zone]);  // CGS -> MKS conversion

    EOS::TY2H(T_bc[zone], &(Y_bc[zone][-1]),h_bc[zone]);  // CGS -> MKS conversion
    }
    
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k, 
                 amrex::Array4<amrex::Real> const& state,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 PmfData const *pmf_data)
{

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    const amrex::Real L_x = prob_hi[0] - prob_lo[0];

    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert, r;
    amrex::Real y1, y2;
    amrex::Real eta, rho, T, h, u, v;
    amrex::Real coor [3];
    amrex::Real velocity [3];

    coor[0] = x; 
    coor[1] = y;
    coor[2] = z;

    pert = 0.0;

    bcnormal(coor/*const amrex::Real* /*x[AMREX_SPACEDIM]*/*/, dx /*amrex::Real s_ext[DEF_NUM_STATE]*/, 1 /*const int idir*/,-1 /*const int sgn*/,time /*,amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  ACParm const& ac_parm,
  PmfData const *pmf_data*/)

    for (int n = 0; n < NUM_SPECIES; n++){
      state(i,j,k,DEF_first_spec+n)= massfrac[n];
    }
    EOS::X2Y(molefrac, massfrac);

    state(i,j,k,DEF_Temp) = T;

    if (prob_parm.blobr>0.0)
    {
     r = sqrt((x-prom_parm.blobx)**2 + (y-prob_parm.bloby)**2);
     eta = 0.5*(1.0 - tanh(2.0*(r-prob_parm.blobr)/prob_parm.blobw));
     state(i,j,k,DEF_Temp) = prob_parm.blobT*eta + state(i,j,k,DEF_Temp)*(1.0-eta);
    }
    state(i,j,k,Xvel) = velocity[0];
#if ( AMREX_SPACEDIM == 2 ) 
    state(i,j,k,Yvel) = velocity[1];
#elif ( AMREX_SPACEDIM == 3 ) 
    state(i,j,k,Zvel) = velocity[2];
#endif

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.Pamb * 10.0;

    EOS::PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    EOS::TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
int getZone(
  const int x,
  const int y)
{
  int getZone=prob_parm.BL_VOLUME
  if (y <= prob_lo[1]) {
      if (abs(x) <= fuel_ox_split) { 
          getZone = prob_parm.BL_FUELPIPE; 
      } else if (abs(x) <= fuel_ox_split+pipeTh) {
          getZone = prob_parm.BL_PIPEEND;
      } else if (abs(x) <= ox_air_split) {
          getZone = prob_parm.BL_OXIDIZER;
      } else { getZone = BL_AIR; }

   
  } else if (y >= prob_hi[1]) {
      getZone = prob_parm.BL_OUTFLOW;}

  return getZone;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  amrex::Real s_ext[DEF_NUM_STATE],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  ACParm const& ac_parm,
  PmfData const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real zone;
  amrex::Real rho, patm;
  amrex::Real T;
  amrex::Real h;
  amrex::Real eta;
  amrex::Real zoneL, zoneR, zoneT;

  if (sgn == 1) {
    zone = getZone(x[0],x[1]);

    if (!ac_parm.ctrl_active || time <= 0.0) {
       if (zone == BL_OUTFLOW) {
          rho= prob_parm.rho_bc[zone];
          for (int n = 0; n < NUM_SPECIES; n++) {
             massfrac[n]=prob_parm.Y_bc[zone][n];
           }
          T = T_bc[zone];
          h = h_bc[zone];
       s_ext[Xvel] = 0.0;
       s_ext[Yvel] = prob_parm.v_bc[zone];
       s_ext[Zvel] = 0.0;
    } else {
         zoneL = BL_FUELPIPE;
         zoneR = BL_OXIDIZER;
         zoneT = BL_VOLUME;
         eta = 0.50*(1.0 - tanh(2.0*(ABS(x[0])-fuel_ox_split)/blobw));
         for (int n = 0; n < NUM_SPECIES; n++) {
            massfrac[n]=prob_parm.Y_bc[zoneL][n]*eta + prob_parm.Y_bc[zoneR][n]*(1.0-eta);
         }       
         T = T_bc[zoneL]*eta + T_bc[zoneR]*(1.0-eta);
         eta = 0.50*(1.0 - tanh(2.0*(x[1]-bloby)/blobw));
         for (int n = 0; n < NUM_SPECIES; n++) {
            massfrac[n]=massfrac[n]*eta + prob_parm.Y_bc[zoneT][n]*(1.0-eta);
         }
         T = T*eta + T_bc[zoneT]*(1.0-eta);

         EOS::X2Y(molefrac, massfrac);
         patm=prob_parm.pamb*10.0;
         EOS::PYT2R(patm, massfrac, T, rho);  // CGS -> MKS conversion
         s_ext[Density] = rho * 1.0e3;

         EOS::TY2H(T,massfrac,h);  // CGS -> MKS conversion
         s_ext[Xvel] = 0.0;

         for (int n = 0; n < NUM_SPECIES; n++) {
           s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
          }

         zone=getZone(x[0],X[1]);

         if (zone==BL_FUELPIPE) {
               zoneL = BL_PIPEEND;
               zoneR = BL_FUELPIPE;
               r = fuel_ox_split - abs(x[0]);
               eta = (1.0 - tanh(2.0*r/pipeBL));
               s_ext[Yvel] = v_bc[zoneR]*2.0**(1-((abs(x[0]))**(2)/(fuel_ox_split**(2))));
               s_ext[Zvel] = 0.0;
    } else if (zone==BL_OXIDIZER) {
               zoneL = BL_PIPEEND;
               zoneR = BL_OXIDIZER;
               if (x[0]<(fuel_ox_split+pipeTh+Vel_zero_th)) {
               s_ext[Yvel] = 0.0;            
         } else {
               r = abs(x[0]) - (fuel_ox_split+pipeTh);
               eta = (1.0 - tanh(2.0*r/pipeBL));
               s_ext[Yvel] = v_bc[zoneL]*eta + v_bc[zoneR]*(1.0-eta);
         }
               s_ext[Zvel] = 0.0;
       } else {
         s_ext[Yvel] = v_bc[zone];
         s_ext[Zvel] = 0.0;
       }
   }
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int /*i*/, int /*j*/, int /*k*/,
           amrex::Array4<amrex::Real> const& /*beta*/,
           amrex::GeometryData const& /*geomdata*/,
           amrex::Box const& /*dBox*/,
           const int  /*dir*/,
           const int  /*state_comp*/,
           const int  /*nComp*/)
{
}
#endif
