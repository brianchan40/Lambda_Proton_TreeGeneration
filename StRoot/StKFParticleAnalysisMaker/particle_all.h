#ifndef particle_all_hh
#define particle_all_hh

class particle_all
{
public:
    particle_all() {
        px = -3.14999;
        py = -3.14999;
        pz = -3.14999;
        Charge = -3.14999;
        dcaglobal = -3.14999;
        nSigmaProton = -3.14999;
        nSigmaPion = -3.14999;
        isTofTrack = false;
        trk_id = -999;
        is_pion = false;
        is_all = false;
        is_proton = false;
        // hits_ratio = -3.149999;
        nhitsfit = -3.1499999;
        nhitsmax = -3.1499999;
    }
    particle_all(float in_px, float in_py, float in_pz, float in_charge, float in_dcaglobal, float in_nSigmaProton, float in_nSigmaPion, bool in_isTofTrack, int in_trk_id, bool in_is_pion, bool in_is_all, bool in_is_proton, /*float in_hits_ratio,*/ float in_nhitsfit, float in_nhitsmax);

    virtual       ~particle_all() { }

    float px, py, pz;
    float Charge;
    float dcaglobal;
    float nSigmaProton;
    float nSigmaPion;
    bool isTofTrack;
    int trk_id;
    bool is_pion, is_all, is_proton;
    // float hits_ratio;
    float nhitsfit;
    float nhitsmax;
};



#endif
