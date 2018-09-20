// Copyright Chun Shen @ 2018

#ifndef SRC_HYDROINFOMUSIC_H_
#define SRC_HYDROINFOMUSIC_H_

#include <vector>
#include <string>
#include "data_struct.h"
#include "data.h"

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

 public:
    HydroinfoMUSIC();       // constructor
    ~HydroinfoMUSIC();      // destructor

    void clean_hydro_event();
    double get_hydro_tau_max() {return(hydroTauMax);}
    double get_hydro_tau0()    {return(hydroTau0);}
    double get_hydro_dtau()    {return(hydroDtau);}
    double get_hydro_dx()      {return(hydroDx);}
    double get_hydro_deta()    {return(hydroDeta);}
    double get_hydro_eta_max() {return(hydro_eta_max);}
    double get_hydro_x_max()   {return(hydroXmax);}

    void readHydroData(int whichHydro, int nskip_tau_in,
            std::string input_filename_in, std::string hydro_ideal_filename,
            std::string hydro_shear_filename, std::string hydro_bulk_filename);

    void getHydroValues(const double x, const double y,
                        const double z, const double t,
                        fluidCell *info);
    void set_grid_infomatioin(const InitData &DATA);
    void print_grid_information();
    void dump_ideal_info_to_memory(double tau, float epsilon, float pressure,
                                   float entropy, float T,
                                   float ux, float uy, float ueta);
    int get_number_of_fluid_cells() {return(lattice_ideal.size());};
};

#endif  // SRC_HYDROINFO_MUSIC_H_

