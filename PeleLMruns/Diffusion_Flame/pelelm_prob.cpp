#include <PeleLM.H>
#include <pelelm_prob.H>
#include <pmf.H>

extern "C" {
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex_real* /*problo*/,
                         const amrex_real* /*probhi*/)
    {
        amrex::ParmParse pp("prob");

        pp.query("Pamb",   PeleLM::prob_parm->P_mean);
        pp.query("fuel_ox_split", PeleLM::prob_parm->fuel_ox_split);
        pp.query("ox_air_split",  PeleLM::prob_parm->ox_air_split);
        pp.query("pipeTh",  PeleLM::prob_parm->pipeTh);
        pp.query("pipeBL",  PeleLM::prob_parm->pipeBL);
        pp.query("blobx",  PeleLM::prob_parm->blobx);
        pp.query("bloby",  PeleLM::prob_parm->bloby);
        pp.query("blobr",  PeleLM::prob_parm->blobr);
        pp.query("blobw",  PeleLM::prob_parm->blobw);
        pp.query("blobT",  PeleLM::prob_parm->blobT);
        pp.query("V_fu",  PeleLM::prob_parm->V_fu);
        pp.query("V_ox",  PeleLM::prob_parm->V_ox);
        pp.query("V_air",  PeleLM::prob_parm->V_air);
        pp.query("T_fu",  PeleLM::prob_parm->T_fu);
        pp.query("T_ox",  PeleLM::prob_parm->T_ox);
        pp.query("T_air",  PeleLM::prob_parm->T_air);

        setupbc(PeleLM::prob_parm.get());

        PeleLM::prob_parm->bathID = N2_ID;  
        PeleLM::prob_parm->fuelID = CH4_ID;  
        PeleLM::prob_parm->oxidID = O2_ID;

        std::string pmf_datafile;
        pp.query("pmf_datafile", pmf_datafile);
    }
}
