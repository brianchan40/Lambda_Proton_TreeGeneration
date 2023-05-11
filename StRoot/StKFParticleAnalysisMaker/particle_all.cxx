#include "particle_all.h"

particle_all::particle_all(float in_px, float in_py, float in_pz, float in_charge, float in_dcaglobal, float in_nSigmaProton, float in_nSigmaPion, bool in_isTofTrack, int in_trk_id, bool in_is_pion, bool in_is_all, bool in_is_proton /*, float in_hits_ratio*/, float in_nhitsfit, float in_nhitsmax)
{
    px = in_px;
    py = in_py;
    pz = in_pz;
    Charge = in_charge;
    dcaglobal = in_dcaglobal;
    nSigmaProton = in_nSigmaProton;
    nSigmaPion = in_nSigmaPion;
    isTofTrack = in_isTofTrack;
    trk_id = in_trk_id;
    is_pion = in_is_pion;
    is_all = in_is_all;
    is_proton = in_is_proton;

    // hits_ratio = in_hits_ratio;
    nhitsfit = in_nhitsfit;
    nhitsmax = in_nhitsmax;
}
