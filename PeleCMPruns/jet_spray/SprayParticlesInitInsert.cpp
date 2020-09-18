
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include "prob.H"

using namespace amrex;

int
interpolateInjectTime(const Real& time)
{
  const int nvals = ProbParm::inject_N;
  int i = 0;
  Real ctime = ProbParm::d_inject_time[i];
  while (ctime < time) {
    ctime = ProbParm::d_inject_time[++i];
  }
  return i;
}

bool
SprayParticleContainer::insertParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level)
{
  return false;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int  nstep,
                                        int  lev,
                                        int  finest_level)
{
  if (lev != 0) return false;
  if (time < ProbParm::jet_start_time
      || time > ProbParm::jet_end_time) return false;
  // Number of particles per parcel
  const Real num_ppp = m_parcelSize;
  const Geometry& geom = this->m_gdb->Geom(lev);
  const auto plo = geom.ProbLoArray();
  const auto phi = geom.ProbHiArray();
  const auto dx = geom.CellSize();
  RealVect dom_len(AMREX_D_DECL(geom.ProbLength(0),
                                geom.ProbLength(1),
                                geom.ProbLength(2)));
  Real mass_flow_rate = ProbParm::mass_flow_rate;
  Real jet_vel = ProbParm::jet_vel;
  Real jet_dia = ProbParm::jet_dia;
  Real jr2 = jet_dia*jet_dia/4.; // Jet radius squared
#if AMREX_SPACEDIM == 3
  Real jet_area = M_PI*jr2;
#else
  Real jet_area = jet_dia;
#endif
  Real part_temp = ProbParm::part_temp;
  Real part_rho = ProbParm::part_rho;
  if (ProbParm::inject_N > 0) {
    const int time_indx = interpolateInjectTime(time);
    const Real time1 = ProbParm::d_inject_time[time_indx - 1];
    const Real time2 = ProbParm::d_inject_time[time_indx];
    const Real mf1 = ProbParm::d_inject_mass[time_indx - 1];
    const Real mf2 = ProbParm::d_inject_mass[time_indx];
    const Real invt = (time - time1)/(time2 - time1);
    mass_flow_rate = mf1 + (mf2 - mf1)*invt;
    const Real jv1 = ProbParm::d_inject_vel[time_indx - 1];
    const Real jv2 = ProbParm::d_inject_vel[time_indx];
    jet_vel = jv1 + (jv2 - jv1)*invt;
  }
  // This absolutely must be included with any injection or insertion
  // function or significant issues will arise
  if (jet_vel*dt/dx[0] > 0.5) {
    Real max_vel = dx[0]*0.5/dt;
    if (ParallelDescriptor::IOProcessor()) {
      std::string warn_msg = "Injection velocity of " + std::to_string(jet_vel)
        + " is reduced to maximum " + std::to_string(max_vel);
      amrex::Warning(warn_msg);
    }
    m_injectVel = jet_vel;
    jet_vel = max_vel;
  }
  Real part_dia = ProbParm::part_mean_dia;
  Real part_stdev = ProbParm::part_stdev_dia;
  Real stdsq = part_stdev*part_stdev;
  Real meansq = part_dia*part_dia;
  Real log_mean = 2.*std::log(part_dia) - 0.5*std::log(stdsq + meansq);
  Real log_stdev = std::sqrt(amrex::max(-2.*std::log(part_dia)
                                        + std::log(stdsq + meansq), 0.));
  RealVect jet_cent(AMREX_D_DECL(dom_len[0]/2., plo[1], dom_len[2]/2.));
  Real Pi_six = M_PI/6.;
  Real spray_angle = ProbParm::spray_angle;
  Real lo_angle = -0.5*spray_angle;
  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    const RealBox& temp = RealBox(bx,geom.CellSize(),geom.ProbLo());
    const Real* xlo = temp.lo();
    const Real* xhi = temp.hi();
    RealVect box_len(AMREX_D_DECL(temp.length(0), 0., temp.length(2)));
    Gpu::HostVector<ParticleType> host_particles;
#ifdef USE_SPRAY_SOA
    std::array<Gpu::HostVector<Real>, NAR_SPR > host_real_attribs;
#endif
    if (xlo[1] == plo[1]) {
      // Box locations relative to jet center
      const RealVect xloJ(AMREX_D_DECL(xlo[0] - jet_cent[0], plo[1],
                                       xlo[2] - jet_cent[2]));
      const RealVect xhiJ(AMREX_D_DECL(xhi[0] - jet_cent[0], plo[1],
                                       xhi[2] - jet_cent[2]));
      Real cur_jet_area = 0.;
      Real testdx = dx[0]/100.;
      Real testdx2 = testdx*testdx;
      Real curx = xloJ[0];
      Real hix = xloJ[0];
      Real lox = xhiJ[0];
      // Loop over each cell and check how much overlap there is with the jet
#if AMREX_SPACEDIM == 3
      Real hiz = xloJ[2];
      Real loz = xhiJ[2];
      while (curx < xhiJ[0]) {
        Real curz = xloJ[2];
        while (curz < xhiJ[2]) {
          Real r2 = curx*curx + curz*curz;
          if (r2 <= jr2) {
            cur_jet_area += testdx2;
            lox = amrex::min(curx, lox);
            hix = amrex::max(curx, hix);
            loz = amrex::min(curz, loz);
            hiz = amrex::max(curz, hiz);
          }
          curz += testdx;
        }
        curx += testdx;
      }
      Real zlen = hiz - loz;
      loz += jet_cent[2];
#else
      while (curx < xhiJ[0]) {
        Real r2 = curx*curx;
        if (r2 <= jr2) {
          cur_jet_area += testdx;
          lox = amrex::min(curx, lox);
          hix = amrex::max(curx, hix);
        }
        curx += testdx;
      }
#endif
      Real xlen = hix - lox;
      lox += jet_cent[0];
      Real jet_perc = cur_jet_area/jet_area;
      Real perc_mass = jet_perc*mass_flow_rate*dt;
      Real total_mass = 0.;
      while (total_mass < perc_mass) {
        RealVect part_loc(AMREX_D_DECL(lox + amrex::Random()*xlen,
                                       plo[1], loz + amrex::Random()*zlen));
        Real r2 = AMREX_D_TERM(std::pow(part_loc[0] - jet_cent[0], 2),,
                               + std::pow(part_loc[2] - jet_cent[2], 2));
        if (r2 <= jr2) {
          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          Real theta = lo_angle + spray_angle*amrex::Random();
#if AMREX_SPACEDIM == 3
          Real theta2 = 2.*M_PI*amrex::Random();
#else
          Real theta2 = 0.;
#endif
          Real x_vel = jet_vel*std::sin(theta)*std::cos(theta2);
          Real y_vel = jet_vel*std::cos(theta);
          Real z_vel = jet_vel*std::sin(theta)*std::sin(theta2);
#ifdef USE_SPRAY_SOA
          AMREX_D_TERM(host_real_attribs[PeleC::pstateVel].push_back(x_vel);,
                       host_real_attribs[PeleC::pstateVel+1].push_back(y_vel);,
                       host_real_attribs[PeleC::pstateVel+2].push_back(z_vel););
#else
          AMREX_D_TERM(p.rdata(PeleC::pstateVel) = x_vel;,
                       p.rdata(PeleC::pstateVel+1) = y_vel;,
                       p.rdata(PeleC::pstateVel+2) = z_vel;);
#endif
          Real cur_dia = amrex::RandomNormal(log_mean, log_stdev);
          // Use a log normal distribution
          cur_dia = std::exp(cur_dia);
          Real part_y = plo[1] + amrex::Random()*dt*y_vel;
          AMREX_D_TERM(p.pos(0) = part_loc[0];,
                       p.pos(1) = part_y;,
                       p.pos(2) = part_loc[2];);
#ifdef USE_SPRAY_SOA
          host_real_attribs[PeleC::pstateT].push_back(part_temp);
          host_real_attribs[PeleC::pstateDia].push_back(cur_dia);
          host_real_attribs[PeleC::pstateRho].push_back(part_rho);
          host_real_attribs[PeleC::pstateY].push_back(1.);
          for (int sp = 1; sp != SPRAY_FUEL_NUM; ++sp)
            host_real_attribs[PeleC::pstateY + sp].push_back(0.);
#else
          p.rdata(PeleC::pstateT) = part_temp;
          p.rdata(PeleC::pstateDia) = cur_dia;
          p.rdata(PeleC::pstateRho) = part_rho;
          for (int sp = 0; sp != SPRAY_FUEL_NUM; ++sp)
            p.rdata(PeleC::pstateY + sp) = 0.;
          p.rdata(PeleC::pstateY) = 1.;
#endif
          host_particles.push_back(p);
          Real pmass = Pi_six*part_rho*std::pow(cur_dia, 3);
          total_mass += num_ppp*pmass;
        }
      }
    }
    auto& particle_tile = GetParticles(lev)[std::make_pair(mfi.index(),
                                                           mfi.LocalTileIndex())];
    auto old_size = particle_tile.GetArrayOfStructs().size();
    auto new_size = old_size + host_particles.size();
    particle_tile.resize(new_size);

    Gpu::copy(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
              particle_tile.GetArrayOfStructs().begin() + old_size);
#ifdef USE_SPRAY_SOA
    for (int i = 0; i != NAR_SPR; ++i) {
      Gpu::copy(Gpu::hostToDevice, host_real_attribs[i].begin(), host_real_attribs[i].end(),
                particle_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
    }
#endif
  }
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(AmrLevel* pelec, const int& lev, const int& num_ppc)
{
  // Start without any particles
  return;
}
