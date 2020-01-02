// Copyright Chun Shen @ 2018

#include "util.h"
#include "HydroinfoMUSIC.h"

HydroinfoMUSIC::HydroinfoMUSIC() {
    hydroTauMax = 0.0;
    itaumax = 0;
}

HydroinfoMUSIC::~HydroinfoMUSIC() {
    clean_hydro_event();
}

void HydroinfoMUSIC::clean_hydro_event() {
    lattice_ideal.clear();
    hydroTauMax = 0.;
    itaumax = 0;
}

void HydroinfoMUSIC::getHydroValues(
        const double x, const double y, const double z, const double t,
        fluidCell* info) {
// For interpolation of evolution files in tau-eta coordinates. Only the
// reading of MUSIC's evolution_xyeta.dat file is implemented here.
// For simplicity, hydro_eta_max refers to MUSIC's eta_size, and similarly for
// hydroDeta; however, x, y, z, and t are as usual to stay compatible with
// MARTINI.
    double tau, eta;
    if (use_tau_eta_coordinate == 1) {
        if (t*t > z*z) {
            tau = sqrt(t*t - z*z);
            eta = 0.5*log((t + z)/(t - z));
        } else {
            tau = 0.;
            eta = 0.;
        }
    } else {
        // if the medium is given in cartesian coordinates
        // set tau and eta to t and z
        tau = t;
        eta = z;
    }

    int ieta = static_cast<int>((hydro_eta_max + eta)/hydroDeta + 0.0001);
    if (boost_invariant == 1) {
        ieta = 0;
    }

    const int itau = static_cast<int>((tau - hydroTau0)/hydroDtau + 0.0001);
    const int ix   = static_cast<int>((hydroXmax + x)/hydroDx + 0.0001);
    const int iy   = static_cast<int>((hydroXmax + y)/hydroDx + 0.0001);

    double xfrac = (x - (static_cast<double>(ix)*hydroDx - hydroXmax))/hydroDx;
    double yfrac = (y - (static_cast<double>(iy)*hydroDx - hydroXmax))/hydroDx;
    double etafrac = (eta/hydroDeta - static_cast<double>(ieta)
                      + 0.5*static_cast<double>(ietamax));
    double taufrac = (tau - hydroTau0)/hydroDtau - static_cast<double>(itau);

    if (ix < 0 || ix >= ixmax) {
        music_message << "[HydroinfoMUSIC::getHydroValues]: "
                      << "WARNING - x out of range x=" << x
                      << ", ix=" << ix << ", ixmax=" << ixmax;
        music_message.flush("warning");
        music_message << "x=" << x << " y=" << y << " eta=" << eta
                      << " ix=" << ix << " iy=" << iy << " ieta=" << ieta;
        music_message.flush("warning");
        music_message << "t=" << t << " tau=" << tau
                      << " itau=" << itau << " itaumax=" << itaumax;
        music_message.flush("warning");
        music_message.flush("warning");

        info->temperature = 0.0;
        info->ed = 0.0;
        info->sd = 0.0;
        info->pressure = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }

    if (iy < 0 || iy >= ixmax) {
        music_message << "[HydroinfoMUSIC::getHydroValues]: "
                      << "WARNING - y out of range, y=" << y << ", iy="  << iy
                      << ", iymax=" << ixmax;
        music_message.flush("warning");
        music_message << "x=" << x << " y=" << y << " eta=" << eta
                      << " ix=" << ix << " iy=" << iy << " ieta=" << ieta;
        music_message.flush("warning");
        music_message << "t=" << t << " tau=" << tau
                      << " itau=" << itau << " itaumax=" << itaumax;
        music_message.flush("warning");
        music_message.flush("warning");

        info->temperature = 0.0;
        info->ed = 0.0;
        info->sd = 0.0;
        info->pressure = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (itau < 0 || itau > itaumax) {
        music_message << "[HydroinfoMUSIC::getHydroValues]: WARNING - "
                      << "tau out of range, itau=" << itau
                      << ", itaumax=" << itaumax;
        music_message.flush("warning");
        music_message << "[HydroinfoMUSIC::getHydroValues]: tau= " << tau
                      << ", hydroTau0 = " << hydroTau0
                      << ", hydroTauMax = " << hydroTauMax
                      << ", hydroDtau = " << hydroDtau;
        music_message.flush("warning");

        info->temperature = 0.0;
        info->ed = 0.0;
        info->sd = 0.0;
        info->pressure = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (ieta < 0 || ieta >= ietamax) {
        music_message << "[HydroinfoMUSIC::getHydroValues]: WARNING - "
                      << "eta out of range, ieta=" << ieta
                      << ", ietamax=" << ietamax;
        music_message.flush("warning");
        info->temperature = 0.0;
        info->ed = 0.0;
        info->sd = 0.0;
        info->pressure = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }

    // The array of positions on the 4-dimensional rectangle:
    int position[2][2][2][2];
    for (int ipx = 0; ipx < 2; ipx++) {
        int px;
        if (ipx == 0 || ix == ixmax-1) {
            px = ix;
        } else {
            px = ix + 1;
        }
        for (int ipy = 0; ipy < 2; ipy++) {
            int py;
            if (ipy == 0 || iy == ixmax-1) {
                py = iy;
            } else {
                py = iy + 1;
            }
            for (int ipeta = 0; ipeta < 2; ipeta++) {
                int peta;
                if (ipeta == 0 || ieta == ietamax-1) {
                    peta = ieta;
                } else {
                    peta = ieta + 1;
                }
                for (int iptau = 0; iptau < 2; iptau++) {
                    int ptau;
                    if (iptau == 0 || itau == itaumax-1) {
                        ptau = itau;
                    } else {
                        ptau = itau + 1;
                    }
                    position[ipx][ipy][ipeta][iptau] = (
                                px + ixmax*(py + ixmax*(peta + ietamax*ptau)));
                }
            }
        }
    }

    // And now, the interpolation:
    double T = 0.0;
    double ed = 0.0;
    double sd = 0.0;
    double pressure = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double ux = 0.0;
    double uy = 0.0;
    double ueta = 0.0;
    double pi00 = 0.0;
    double pi01 = 0.0;
    double pi02 = 0.0;
    double pi03 = 0.0;
    double pi11 = 0.0;
    double pi12 = 0.0;
    double pi13 = 0.0;
    double pi22 = 0.0;
    double pi23 = 0.0;
    double pi33 = 0.0;
    double bulkPi = 0.0;

    fluidCell_ideal *HydroCell_ptr1, *HydroCell_ptr2;
    for (int iptau = 0; iptau < 2; iptau++) {
        double taufactor;
        if (iptau == 0)
            taufactor = 1. - taufrac;
        else
            taufactor = taufrac;
        for (int ipeta = 0; ipeta < 2; ipeta++) {
            double etafactor;
            if (ipeta == 0)
                etafactor = 1. - etafrac;
            else
                etafactor = etafrac;
            for (int ipy = 0; ipy < 2; ipy++) {
                double yfactor;
                if (ipy == 0)
                    yfactor = 1. - yfrac;
                else
                    yfactor = yfrac;

                double prefrac = yfactor*etafactor*taufactor;

                HydroCell_ptr1 = (
                        &lattice_ideal[position[0][ipy][ipeta][iptau]]);
                HydroCell_ptr2 = (
                        &lattice_ideal[position[1][ipy][ipeta][iptau]]);
                ed += prefrac*((1. - xfrac)*HydroCell_ptr1->ed
                              + xfrac*HydroCell_ptr2->ed);
                sd += prefrac*((1. - xfrac)*HydroCell_ptr1->sd
                              + xfrac*HydroCell_ptr2->sd);
                pressure += prefrac*((1. - xfrac)*HydroCell_ptr1->pressure
                                     + xfrac*HydroCell_ptr2->pressure);
                T += prefrac*((1. - xfrac)*HydroCell_ptr1->temperature
                              + xfrac*HydroCell_ptr2->temperature);
                ux += prefrac*((1. - xfrac)*HydroCell_ptr1->ux
                                + xfrac*HydroCell_ptr2->ux);
                uy += prefrac*((1. - xfrac)*HydroCell_ptr1->uy
                                + xfrac*HydroCell_ptr2->uy);
                ueta += prefrac*((1. - xfrac)*HydroCell_ptr1->ueta
                                 + xfrac*HydroCell_ptr2->ueta);
            }
        }
    }

    double eta_local = 0.5*log((t + z)/(t - z));
    double sinh_eta, cosh_eta;
    if (fabs(eta_local) < 1e-6) {
        // use Taylor expansion for small eta_s to speed up
        // avoiding to evaluate sinh and cosh
        sinh_eta = eta_local;
        cosh_eta = 1.0 + 0.5*eta_local*eta_local;
    } else {
        sinh_eta = sinh(eta_local);
        cosh_eta = cosh(eta_local);
    }
    double utau = sqrt(1. + ux*ux + uy*uy + ueta*ueta);
    double uz = utau*sinh_eta + ueta*cosh_eta;
    double ut = utau*cosh_eta + ueta*sinh_eta;
    vx = ux/ut;
    vy = uy/ut;
    vz = uz/ut;

    info->temperature = static_cast<float>(T);
    info->vx = static_cast<float>(vx);
    info->vy = static_cast<float>(vy);
    info->vz = static_cast<float>(vz);

    info->ed = static_cast<float>(ed);
    info->sd = static_cast<float>(sd);
    info->pressure = static_cast<float>(pressure);

    info->pi[0][0] = static_cast<float>(pi00);
    info->pi[0][1] = static_cast<float>(pi01);
    info->pi[0][2] = static_cast<float>(pi02);
    info->pi[0][3] = static_cast<float>(pi03);
    info->pi[1][0] = static_cast<float>(pi01);
    info->pi[1][1] = static_cast<float>(pi11);
    info->pi[1][2] = static_cast<float>(pi12);
    info->pi[1][3] = static_cast<float>(pi13);
    info->pi[2][0] = static_cast<float>(pi02);
    info->pi[2][1] = static_cast<float>(pi12);
    info->pi[2][2] = static_cast<float>(pi22);
    info->pi[2][3] = static_cast<float>(pi23);
    info->pi[3][0] = static_cast<float>(pi03);
    info->pi[3][1] = static_cast<float>(pi13);
    info->pi[3][2] = static_cast<float>(pi23);
    info->pi[3][3] = static_cast<float>(pi33);

    info->bulkPi = static_cast<float>(bulkPi);
}


void HydroinfoMUSIC::get_fluid_cell_with_index(const int idx,
                                               fluidCell *info) const {
    info->temperature = static_cast<float>(lattice_ideal[idx].temperature);
    
    double ux   = lattice_ideal[idx].ux;
    double uy   = lattice_ideal[idx].uy;
    double ueta = lattice_ideal[idx].ueta;
    double utau = sqrt(1. + ux*ux + uy*uy + ueta*ueta);
    double sinh_eta = sinh(lattice_ideal[idx].eta);
    double cosh_eta = cosh(lattice_ideal[idx].eta);
    double uz = utau*sinh_eta + ueta*cosh_eta;
    double ut = utau*cosh_eta + ueta*sinh_eta;
    info->vx = static_cast<float>(ux/ut);
    info->vy = static_cast<float>(uy/ut);
    info->vz = static_cast<float>(uz/ut);

    info->ed = static_cast<float>(lattice_ideal[idx].ed);
    info->sd = static_cast<float>(lattice_ideal[idx].sd);
    info->pressure = static_cast<float>(lattice_ideal[idx].pressure);

    info->pi[0][0] = static_cast<float>(0.0);
    info->pi[0][1] = static_cast<float>(0.0);
    info->pi[0][2] = static_cast<float>(0.0);
    info->pi[0][3] = static_cast<float>(0.0);
    info->pi[1][0] = static_cast<float>(0.0);
    info->pi[1][1] = static_cast<float>(0.0);
    info->pi[1][2] = static_cast<float>(0.0);
    info->pi[1][3] = static_cast<float>(0.0);
    info->pi[2][0] = static_cast<float>(0.0);
    info->pi[2][1] = static_cast<float>(0.0);
    info->pi[2][2] = static_cast<float>(0.0);
    info->pi[2][3] = static_cast<float>(0.0);
    info->pi[3][0] = static_cast<float>(0.0);
    info->pi[3][1] = static_cast<float>(0.0);
    info->pi[3][2] = static_cast<float>(0.0);
    info->pi[3][3] = static_cast<float>(0.0);

    info->bulkPi = static_cast<float>(0.0);
}


void HydroinfoMUSIC::set_grid_infomatioin(const InitData &DATA) {
    use_tau_eta_coordinate = 1;
    boost_invariant = DATA.boost_invariant;

    hydroTau0 = DATA.tau0;
    hydroDtau = DATA.delta_tau*DATA.output_evolution_every_N_timesteps;
    hydroDx   = DATA.delta_x*DATA.output_evolution_every_N_x;
    hydroDeta = DATA.delta_eta*DATA.output_evolution_every_N_eta;

    hydroXmax     = DATA.x_size/2.;
    hydro_eta_max = DATA.eta_size/2.;

    ixmax   = (static_cast<int>((DATA.nx - 1)
                                /DATA.output_evolution_every_N_x) + 1);
    ietamax = (static_cast<int>((DATA.neta - 1)
                                /DATA.output_evolution_every_N_eta) + 1);
}

void HydroinfoMUSIC::print_grid_information() {
    music_message << "boost_invariant = " << boost_invariant;
    music_message << "hydro_tau0 = " << hydroTau0 << " fm";
    music_message << "hydro_tau_max = " << hydroTauMax << " fm";
    music_message << "hydry_dtau = " << hydroDtau << " fm";
    music_message << "itaumax = " << itaumax;
    music_message << "hydro_Xmax = " << hydroXmax << " fm";
    music_message << "hydro_dx = " << hydroDx << " fm";
    music_message << "ixmax = " << ixmax;
    music_message << "hydro_eta_max = " << hydro_eta_max;
    music_message << "hydro_deta = " << hydroDeta;
    music_message << "ietamax = " << ietamax;
    music_message.flush("info");
}

void HydroinfoMUSIC::dump_ideal_info_to_memory(double tau, float eta,
        float epsilon, float pressure, float entropy, float T,
        float ux, float uy, float ueta) {
    if (tau > hydroTauMax) {
        hydroTauMax = tau;
        itaumax++;
    }
    fluidCell_ideal new_cell;
    new_cell.eta = eta;
    new_cell.temperature = T*Util::hbarc;
    new_cell.ed = epsilon*Util::hbarc;
    new_cell.sd = entropy;
    new_cell.pressure = pressure*Util::hbarc;
    new_cell.ux = ux;
    new_cell.uy = uy;
    new_cell.ueta = ueta;
    lattice_ideal.push_back(new_cell);
}
