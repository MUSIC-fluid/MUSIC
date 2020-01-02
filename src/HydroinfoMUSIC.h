// Copyright Chun Shen @ 2018

#ifndef SRC_HYDROINFOMUSIC_H_
#define SRC_HYDROINFOMUSIC_H_

#include <vector>
#include <string>
#include "data_struct.h"
#include "data.h"
#include "pretty_ostream.h"

class HydroinfoMUSIC {
 private:
    double hydroTau0;       // tau_0 in the hydro data files
    double hydroTauMax;     // tau_max in the hydro data files
    double hydroDtau;       // step dtau in fm/c in the hydro data files
    double hydroXmax;       // maximum x in fm in the hydro data files
                            // [-xmax, +xmax] for both x and y
    double hydro_eta_max;   // maximum z in fm in the hydro data files
                            // [-zmax, +zmax] for 3D hydro
    double hydroDx;         // step dx in fm in the hydro data files
    double hydroDeta;       // step dz in fm in the hydro data files in
                            // the z-direction for 3D hydro

    int use_tau_eta_coordinate;
    bool boost_invariant;

    int itaumax, ixmax, ietamax;

    std::vector<fluidCell_ideal> lattice_ideal;
    pretty_ostream music_message;

 public:
    HydroinfoMUSIC();       // constructor
    ~HydroinfoMUSIC();      // destructor

    void clean_hydro_event();
    double get_hydro_tau_max() const {return(hydroTauMax);}
    double get_hydro_tau0() const    {return(hydroTau0);}
    double get_hydro_dtau() const    {return(hydroDtau);}
    double get_hydro_dx() const      {return(hydroDx);}
    double get_hydro_deta() const    {return(hydroDeta);}
    double get_hydro_eta_max() const {return(hydro_eta_max);}
    double get_hydro_x_max() const   {return(hydroXmax);}
    int get_ntau() const {return(itaumax);}
    int get_nx()   const {return(ixmax  );}
    int get_neta() const {return(ietamax);}
    bool is_boost_invariant() const {return(boost_invariant);}

    void getHydroValues(const double x, const double y,
                        const double z, const double t,
                        fluidCell *info);
    void set_grid_infomatioin(const InitData &DATA);
    void print_grid_information();
    void dump_ideal_info_to_memory(double tau, float eta, float epsilon,
                                   float pressure, float entropy, float T,
                                   float ux, float uy, float ueta);

    int get_number_of_fluid_cells() const {return(lattice_ideal.size());}
    void get_fluid_cell_with_index(const int idx, fluidCell *info) const;
};

#endif  // SRC_HYDROINFO_MUSIC_H_

