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
#include <PeleLM.H>

#include <NavierStokesBase.H>
#include <IndexDefines.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setupbc(ProbParm *prob_parm) 
{
    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0}; 

    /*A diffusion flame*/

    int fuelZone = prob_parm->BL_FUELPIPE;
    int oxZone   = prob_parm->BL_OXIDIZER;
    int airZone  = prob_parm->BL_AIR;
    int volZone  = prob_parm->BL_VOLUME;
    int peZone   = prob_parm->BL_PIPEEND;
    int ofZone   = prob_parm->BL_OUTFLOW;
    
    /*Fuel*/

    for (int n = 0; n < NUM_SPECIES; n++){
      massfrac[n] = 0.0;
    }

    massfrac[CH4_ID] = 1.0;
    massfrac[O2_ID] = 0.0;
    massfrac[H2O_ID] = 0.0;
    massfrac[N2_ID] = 1.0 - massfrac[CH4_ID] - massfrac[O2_ID] - massfrac[H2O_ID];

    for (int n = 0; n < NUM_SPECIES; n++){
      prob_parm->Y_bc[fuelZone][n] = massfrac[n];
    }
    prob_parm->T_bc[fuelZone] = prob_parm->T_fu;
    prob_parm->u_bc[fuelZone] = 0.0;
    prob_parm->v_bc[fuelZone] = prob_parm->V_fu;
    prob_parm->w_bc[fuelZone] = 0.0;

    /*Oxidizer*/
    for (int n = 0; n < NUM_SPECIES; n++){
      massfrac[n] = 0.0;
    }
    massfrac[O2_ID] = 0.2395;
    massfrac[CO2_ID] = 0.000;
    massfrac[N2_ID] = 1.0 - massfrac[O2_ID] - massfrac[CO2_ID];

    for (int n = 0; n < NUM_SPECIES; n++){
      prob_parm->Y_bc[oxZone][n] = massfrac[n];
    }
    prob_parm->T_bc[oxZone] = prob_parm->T_ox;
    prob_parm->u_bc[oxZone] = 0.0;
    prob_parm->v_bc[oxZone] = prob_parm->V_ox;
    prob_parm->w_bc[oxZone] = 0.0;

    /*air*/
    for (int n = 0; n < NUM_SPECIES; n++){
      massfrac[n] = 0.0;
    }
    massfrac[O2_ID] = 0.2329;
    massfrac[N2_ID] = 1.0 - massfrac[O2_ID];

    for (int n = 0; n < NUM_SPECIES; n++){
      prob_parm->Y_bc[airZone][n] = massfrac[n];
    }
    prob_parm->T_bc[airZone] = prob_parm->T_air;
    prob_parm->u_bc[airZone] = 0.0;
    prob_parm->v_bc[airZone] = prob_parm->V_air;
    prob_parm->w_bc[airZone] = 0.0;

    /*pipeend*/
    for (int n = 0; n < NUM_SPECIES; n++){
      prob_parm->Y_bc[peZone][n] = prob_parm->Y_bc[oxZone][n];
    }
    prob_parm->T_bc[peZone] = prob_parm->T_bc[oxZone];
    prob_parm->u_bc[peZone] = prob_parm->u_bc[oxZone];
    prob_parm->v_bc[peZone] = 0.0;
    prob_parm->w_bc[peZone] = 0.0;

    /*volume(as air)*/
    for (int n = 0; n < NUM_SPECIES; n++){
      prob_parm->Y_bc[volZone][n] = prob_parm->Y_bc[airZone][n];
    }
    prob_parm->T_bc[volZone] = prob_parm->T_bc[airZone];
    prob_parm->u_bc[volZone] = prob_parm->u_bc[airZone];
    prob_parm->v_bc[volZone] = prob_parm->v_bc[airZone];
    prob_parm->w_bc[volZone] = 0.0;

    /*outflow(as air)*/
    for (int n = 0; n < NUM_SPECIES; n++){
      prob_parm->Y_bc[ofZone][n] = prob_parm->Y_bc[airZone][n];
    }
    prob_parm->T_bc[ofZone] = prob_parm->T_bc[airZone];
    prob_parm->u_bc[ofZone] = prob_parm->u_bc[airZone];
    prob_parm->v_bc[ofZone] = prob_parm->v_bc[airZone];
    prob_parm->w_bc[ofZone] = 0.0;

    amrex::Real P_cgs = prob_parm->P_mean * 10.0;
    for (int zone = 0; zone < prob_parm->BL_NZONES; zone++){
      eos.PYT2R(P_cgs, &(prob_parm->Y_bc[zone][0]), prob_parm->T_bc[zone], prob_parm->rho_bc[zone]);  // CGS -> MKS conversion
      eos.TY2H(prob_parm->T_bc[zone], &(prob_parm->Y_bc[zone][0]),prob_parm->h_bc[zone]);  // CGS -> MKS conversion
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
int getZone(amrex::GeometryData const& geomdata,
            ProbParm const& prob_parm,
            const amrex::Real x,
            const amrex::Real y)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  int getZone=prob_parm.BL_VOLUME;

  if (y <= prob_lo[1]) {
      if (abs(x) <= prob_parm.fuel_ox_split) { 
          getZone = prob_parm.BL_FUELPIPE; 
      } else if (abs(x) <= prob_parm.fuel_ox_split+prob_parm.pipeTh) {
          getZone = prob_parm.BL_PIPEEND;
      } else if (abs(x) <= prob_parm.ox_air_split) {
          getZone = prob_parm.BL_OXIDIZER;
      } else { 
         getZone = prob_parm.BL_AIR;
      }
  } else if (y >= prob_hi[1]) {
      getZone = prob_parm.BL_OUTFLOW;}

  return getZone;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
         const amrex::Real x[AMREX_SPACEDIM],
         amrex::Real s_ext[DEF_NUM_STATE],
         const int idir,
         const int sgn,
         const amrex::Real time,
         amrex::GeometryData const& geomdata,
         ProbParm const& prob_parm,
         ACParm const& ac_parm,
         PmfData const *pmf_data)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real zone;
  amrex::Real rho, patm;
  amrex::Real T,h;
  amrex::Real eta;
  amrex::Real zoneL, zoneR, zoneT;

  if (sgn == 1) {
    zone = getZone(geomdata,prob_parm,x[0],x[1]);
    if (!ac_parm.ctrl_active || time <= 0.0) {
      if (zone == prob_parm.BL_OUTFLOW) {
        rho= prob_parm.rho_bc[zone];
        for (int n = 0; n < NUM_SPECIES; n++) {
           massfrac[n]=prob_parm.Y_bc[zone][n];
        }
        s_ext[DEF_Temp] = prob_parm.T_bc[zone];
        h = prob_parm.h_bc[zone];
        s_ext[Density] = rho * 1.0e3;
        s_ext[DEF_RhoH] = h * 1.0e-4 * s_ext[Density];
        s_ext[Xvel] = 0.0;
        s_ext[Yvel] = prob_parm.v_bc[zone];
        for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
        }
#if ( AMREX_SPACEDIM == 3 )
        s_ext[Zvel] = 0.0;
#endif
    } else {
        zoneL = prob_parm.BL_FUELPIPE;
        zoneR = prob_parm.BL_OXIDIZER;
        zoneT = prob_parm.BL_VOLUME;
        eta = 0.50*(1.0 - tanh(2.0*(abs(x[0])-prob_parm.fuel_ox_split)/prob_parm.blobw));
        for (int n = 0; n < NUM_SPECIES; n++) {
            massfrac[n]=prob_parm.Y_bc[zoneL][n]*eta + prob_parm.Y_bc[zoneR][n]*(1.0-eta);
        }       
        T = prob_parm.T_bc[zoneL]*eta + prob_parm.T_bc[zoneR]*(1.0-eta);
        eta = 0.50*(1.0 - tanh(2.0*(x[1]-prob_parm.bloby)/prob_parm.blobw));
        for (int n = 0; n < NUM_SPECIES; n++) {
           massfrac[n]=massfrac[n]*eta + prob_parm.Y_bc[zoneT][n]*(1.0-eta);
        }
        amrex::Real RhoH_temp;
        s_ext[DEF_Temp] = T*eta + prob_parm.T_bc[zoneT]*(1.0-eta);

        eos.Y2X(massfrac, molefrac);
        patm = prob_parm.P_mean*10.0;
        eos.PYT2R(patm, massfrac, s_ext[DEF_Temp], rho);  // CGS -> MKS conversion
        s_ext[Density] = rho * 1.0e3;

        eos.TY2H(s_ext[DEF_Temp],massfrac,RhoH_temp);  // CGS -> MKS conversion
        s_ext[DEF_RhoH] = RhoH_temp * 1.0e-4 * s_ext[Density];   // CGS -> MKS conversion
        s_ext[Xvel] = 0.0;

        for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[DEF_first_spec+n] = massfrac[n] * s_ext[Density];
        }

        if ( zone == prob_parm.BL_FUELPIPE) {
          zoneL = prob_parm.BL_PIPEEND;
          zoneR = prob_parm.BL_FUELPIPE;
          amrex::Real r = prob_parm.fuel_ox_split - abs(x[0]);
          eta = (1.0 - tanh(2.0*r/prob_parm.pipeBL));
          s_ext[Yvel] = prob_parm.v_bc[zoneR] * std::pow(2.0, ( 1.0 - ( (abs(x[0]))*(abs(x[0]))
                                                                        / (prob_parm.fuel_ox_split*prob_parm.fuel_ox_split) ) ) );
#if ( AMREX_SPACEDIM == 3 )
          s_ext[Zvel] = 0.0;
#endif
    } else if (zone==prob_parm.BL_OXIDIZER) {
          zoneL = prob_parm.BL_PIPEEND;
          zoneR = prob_parm.BL_OXIDIZER;
          if (x[0]<(prob_parm.fuel_ox_split+prob_parm.pipeTh+prob_parm.Vel_zero_th)) {
            s_ext[Yvel] = 0.0;            
          } else {
            amrex::Real r = abs(x[0]) - (prob_parm.fuel_ox_split+prob_parm.pipeTh);
            eta = (1.0 - tanh(2.0*r/prob_parm.pipeBL));
            s_ext[Yvel] = prob_parm.v_bc[zoneL]*eta + prob_parm.v_bc[zoneR]*(1.0-eta);
          }
#if ( AMREX_SPACEDIM == 3 )
          s_ext[Zvel] = 0.0;
#endif
        } else {
          s_ext[Yvel] = prob_parm.v_bc[zone];
#if ( AMREX_SPACEDIM == 3 )
          s_ext[Zvel] = 0.0;
#endif
        }
      }  // if !prob_parm.BL_OUTFLOW
    }  // if !active control
  } // if sign
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

    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real r;
    amrex::Real eta;
    const amrex::Real coor[AMREX_SPACEDIM] = {AMREX_D_DECL(x,y,z)};
    amrex::Real s_ext[DEF_NUM_STATE] = {0.0};


    bcnormal(coor, s_ext, 0, 1, 0, geomdata, prob_parm, *(PeleLM::ac_parm.get()), pmf_data);

    for (int n = 0; n < NUM_SPECIES; n++){
      
      massfrac[n]= s_ext[DEF_first_spec+n]/s_ext[Density];
    }

    state(i,j,k,DEF_Temp) = s_ext[DEF_Temp];

    if (prob_parm.blobr>0.0)
    {
       r = sqrt((x-prob_parm.blobx)*(x-prob_parm.blobx) + (y-prob_parm.bloby)*(y-prob_parm.bloby));
       eta = 0.5*(1.0 - tanh(2.0*(r-prob_parm.blobr)/prob_parm.blobw));
       state(i,j,k,DEF_Temp) = prob_parm.blobT*eta + state(i,j,k,DEF_Temp)*(1.0-eta);
    }
    state(i,j,k,Xvel) = s_ext[Xvel];
#if ( AMREX_SPACEDIM == 2 ) 
    state(i,j,k,Yvel) = s_ext[Yvel];
#elif ( AMREX_SPACEDIM == 3 ) 
    state(i,j,k,Zvel) = s_ext[Zvel];
#endif

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i,j,k,DEF_Temp), rho_cgs);
    state(i,j,k,Density) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,DEF_Temp), massfrac, state(i,j,k,DEF_RhoH));
    state(i,j,k,DEF_RhoH) = state(i,j,k,DEF_RhoH) * 1.0e-4 * state(i,j,k,Density);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,DEF_first_spec+n) = massfrac[n] * state(i,j,k,Density);
    }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata /*probhi,problo,coor*/,
           amrex::Box const& dBox,
           const int  dir,
           const int  state_comp,
           const int  nComp)
{
      int Tid, RHid, YSid, YEid, ys, ye;
      int xl, xr, xh, y;

      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* dx      = geomdata.CellSize();

      Tid  = DEF_Temp  - state_comp + AMREX_SPACEDIM;
      RHid = DEF_RhoH    - state_comp + AMREX_SPACEDIM;
      YSid = DEF_first_spec - state_comp + AMREX_SPACEDIM;
      YEid = (DEF_first_spec+NUM_SPECIES-1)  - state_comp + AMREX_SPACEDIM;

      ys = std::max(YSid,1);
      ye = std::min(YEid,nComp);

      if ((dir==1) && (j <= geomdata.domain.smallEnd(1)/*domlo(2)*/) && (((Tid >= 1) && (Tid  <= nComp)) || ((RHid >= 1) && (RHid <= nComp)) || ((ye - ys + 1) >= 1)) ) {
  
         y = (float) j*dx[1]+prob_lo[1];
         xl = (float) i*dx[0]+prob_lo[0]; 
         xr = ((float) i +1.0)*dx[0]+prob_lo[0]; 
         xh = 0.50*(xl+xr);
         if ((getZone(geomdata,*(PeleLM::prob_parm.get()),xl,y)== PeleLM::prob_parm->BL_PIPEEND) ||
             (getZone(geomdata,*(PeleLM::prob_parm.get()),xh,y)== PeleLM::prob_parm->BL_PIPEEND) ||
             (getZone(geomdata,*(PeleLM::prob_parm.get()),xr,y)== PeleLM::prob_parm->BL_PIPEEND)  ) {
              if (((ye - ys + 1) >= 1)) {
                  for (int n =ys; ye;) {
                     beta(i,j,k,n) = 0.0;
                  }
             } // endif
           } //endif
     }//endif 

}
#endif 
