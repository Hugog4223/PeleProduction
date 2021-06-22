#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

#include "PeleC.H"
#include "prob.H"
#include "turbinflow.H"
#include "PelePhysics.H"
#include "IndexDefines.H"

AMREX_GPU_DEVICE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
set_loc(amrex::GeometryData const& geom,
        const amrex::IntVect&      iv)
{
  const amrex::Real* prob_lo = geom.ProbLo();
  const amrex::Real* dx = geom.CellSize();
  return {AMREX_D_DECL(
      prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
      prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
      prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};
}

struct PCHypFillExtDir
{
  ProbParmDevice const* probparmDD;

  AMREX_GPU_HOST
  constexpr explicit PCHypFillExtDir(const ProbParmDevice* d_prob_parm)
    : probparmDD(d_prob_parm)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& dest,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& geom,
    const amrex::Real time,
    const amrex::BCRec* bcr,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const auto x = set_loc(geom,iv);

    const auto* bc = bcr->data();
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real molefrac[NUM_SPECIES], massfrac[NUM_SPECIES], rho, T, e;
    constexpr int dim = AMREX_SPACEDIM;

    /*
      These are the 6 tests for whether we are on a Dirichlet boundary that needs to be filled:
 
        XLO:  if (           (bc[0]     == amrex::BCType::ext_dir) && (iv[0] < domlo[0]))
        XHI:  if (           (bc[0+dim] == amrex::BCType::ext_dir) && (iv[0] > domhi[0]))
        YLO:  if ((dim>1) && (bc[1]     == amrex::BCType::ext_dir) && (iv[1] < domlo[1]))
        YHI:  if ((dim>1) && (bc[1+dim] == amrex::BCType::ext_dir) && (iv[1] > domhi[1]))
        ZLO:  if ((dim>2) && (bc[2]     == amrex::BCType::ext_dir) && (iv[2] < domlo[2]))
        ZHI:  if ((dim>2) && (bc[2+dim] == amrex::BCType::ext_dir) && (iv[2] > domhi[2]))
    */

    // Here, set inflow on ZHI
    // This is the “U-Turn” geometry, where for large y boundary is inflow and for small y boundary is outflow
    //   (we cheat here and make outflow = FOEXTRAP)
    if ((dim > 2) && (bc[2+dim] == amrex::BCType::ext_dir) && (iv[2] > domhi[2])) {
      const auto* prob_hi = geom.ProbHi();
      int dir = 2;
      if (x[1] > 0.5*prob_hi[1]-0.1) {
        amrex::GpuArray<amrex::Real,4+NUM_SPECIES> pmf_vals;
        amrex::GpuArray<amrex::Real,dim> u = {{0.0}};
        pmf(prob_hi[2], prob_hi[2], pmf_vals, *probparmDD);
        for (int n = 0; n < NUM_SPECIES; n++) {
          molefrac[n] = pmf_vals[3 + n];
        }
        auto T = pmf_vals[0];
        auto pres = probparmDD->pamb;

        if (x[1] > 0.5*prob_hi[1]-0.1) {
          u[dim - 1] = pmf_vals[1];
          eos.X2Y(molefrac, massfrac);
          eos.PYT2RE(pres, massfrac, T, rho, e);

          if (probparmDD->turb_ok[dir+dim]) {
            for (int n=0; n<dim; ++n) {
              u[n] += probparmDD->turbarr[dir+dim](iv[0],iv[1],iv[2],n);
            }
          }
          dest(iv,URHO) = rho;
          dest(iv,UMX) = rho * u[0];
          dest(iv,UMY) = rho * u[1];
          dest(iv,UMZ) = rho * u[2];
          dest(iv,UEINT) = rho * e;
          dest(iv,UEDEN) = rho * (e + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
          dest(iv,UTEMP) = T;
          for (int n = 0; n < NUM_SPECIES; n++) {
            dest(iv,UFS+n) = rho * massfrac[n];
          }
        }
        else {
          amrex::IntVect ivi(iv[0],iv[1],iv[2]-1);
          for (int n = 0; n < NVAR; n++) {
            dest(iv,n) = dest(ivi,n); // FOEXTRAP: Copy values from just inside
          }
        }
      }
    }
  }
};

struct PCReactFillExtDir
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dest*/,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* /*bcr*/,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
  }
};

void
pc_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{
  ProbParmDevice* probparmDD = PeleC::d_prob_parm_device; // probparm data for device
  ProbParmDevice* probparmDH = PeleC::h_prob_parm_device; // host copy of probparm data for device
  ProbParmHost* probparmH = PeleC::prob_parm_host;        // probparm data for host
  constexpr int dim = AMREX_SPACEDIM;

  if (probparmH->do_turb) {

    // Copy problem parameter structs to host
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, probparmDD, probparmDD + 1, probparmDH);

    for (int dir=0; dir<dim; ++dir) {
      auto bndryBoxLO = amrex::Box(amrex::adjCellLo(geom.Domain(),dir) & bx);
      if (bcr[1].lo()[dir]==EXT_DIR && bndryBoxLO.ok()) {

        probparmH->turbfab[dir].resize(bndryBoxLO,dim);
        probparmH->turbfab[dir].setVal(0);
        add_turb(bndryBoxLO, probparmH->turbfab[dir], 0, geom, time, dir, amrex::Orientation::low, probparmDH->tp);
        probparmDH->turbarr[dir] = probparmH->turbfab[dir].array();
        probparmDH->turb_ok[dir] = true;

      }

      auto bndryBoxHI = amrex::Box(amrex::adjCellHi(geom.Domain(),dir) & bx);
      if (bcr[1].hi()[dir]==EXT_DIR && bndryBoxHI.ok()) {

        probparmH->turbfab[dir+dim].resize(bndryBoxHI,dim);
        probparmH->turbfab[dir+dim].setVal(0);
        add_turb(bndryBoxHI, probparmH->turbfab[dir+dim], 0, geom, time, dir, amrex::Orientation::high, probparmDH->tp);
        probparmDH->turbarr[dir+dim] = probparmH->turbfab[dir].array();
        probparmDH->turb_ok[dir+dim] = true;

      }
    }

    // Copy problem parameter structs back to device
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, probparmDH, probparmDH + 1, probparmDD);
  }

  amrex::GpuBndryFuncFab<PCHypFillExtDir> hyp_bndry_func(PCHypFillExtDir{probparmDD});
  hyp_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);

  if (probparmH->do_turb) {

    // Copy problem parameter structs to host
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, probparmDD, probparmDD + 1, probparmDH);

    for (int dir=0; dir<dim; ++dir) {
      if (probparmDH->turb_ok[dir]) {
        probparmH->turbfab[dir].clear();
        probparmDH->turb_ok[dir] = false;
      }
      if (probparmDH->turb_ok[dir+dim]) {
        probparmH->turbfab[dir+dim].clear();
        probparmDH->turb_ok[dir+dim] = false;
      }
    }

    // Copy problem parameter structs back to device
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, probparmDH, probparmDH + 1, probparmDD);
    
  }
}

#ifdef PELEC_USE_REACTIONS
void
pc_reactfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp)
{
  amrex::GpuBndryFuncFab<PCReactFillExtDir> react_bndry_func(
    PCReactFillExtDir{});
  react_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}
#endif

void
pc_nullfill(
  amrex::Box const& /*bx*/,
  amrex::FArrayBox& /*data*/,
  const int /*dcomp*/,
  const int /*numcomp*/,
  amrex::Geometry const& /*geom*/,
  const amrex::Real /*time*/,
  const amrex::Vector<amrex::BCRec>& /*bcr*/,
  const int /*bcomp*/,
  const int /*scomp*/)
{
}