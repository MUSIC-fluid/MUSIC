// Copyright 2019 Chun Shen

#include "hydro_source_bullet.h"
//#include "data_struct.h"

HydroSourceBullet::HydroSourceBullet(const InitData &DATA_in) :
                      tau_(DATA_in.bullet_initial_time), DATA_(DATA_in),
                      point_flagged(DATA_in.nx, DATA_in.ny, DATA_in.neta){


    r_ = {DATA_.bullet_position_x, DATA_.bullet_position_y, DATA_.bullet_position_eta_s};
    pmu_ = {DATA_.total_current_energy,
            DATA_.total_current_momentum_x,
            DATA_.total_current_momentum_y,
            DATA_.total_current_momentum_eta_s};

    

    set_source_tau_max(DATA_.bullet_final_time);
    set_source_tau_min(DATA_.bullet_initial_time);
    set_sigma_x(  DATA_.bullet_size_x);
    set_sigma_y(  DATA_.bullet_size_y);
    set_sigma_eta(DATA_.bullet_size_eta_s);
    
    point_flagged = TagPointsForFillEllipse(r_[0], r_[1], r_[2]);

    //Loop over grid totalizing number of points where we will add energy
    int npoints = 0;
    for (int ieta=0; ieta<DATA_.neta; ++ieta)
        for (int iy=0; iy<DATA_.ny; ++iy)
            for (int ix=0; ix<DATA_.nx; ++ix)
                if( point_flagged(ix,iy,ieta) )
                    npoints++;

    //Divide energy equally between all points
    pmu_frac_ = {
        pmu_[0]/((double) npoints),
        pmu_[1]/((double) npoints),
        pmu_[2]/((double) npoints),
        pmu_[3]/((double) npoints)
    };

    music_message << "Using bullet source term.";
    music_message.flush("info");
    music_message << "Energy per cell: " << pmu_frac_[0];
    music_message.flush("info");
    music_message << "px per cell: " << pmu_frac_[1];
    music_message.flush("info");
    music_message << "py per cell: " << pmu_frac_[2];
    music_message.flush("info");
    music_message << "peta per cell: " << pmu_frac_[3];
    music_message.flush("info");

};

void HydroSourceBullet::get_hydro_energy_source(
    const double tau, const double x, const double y, const double eta_s,
    const FlowVec &u_mu, EnergyFlowVec &j_mu) {
    const double dtau = DATA_.delta_tau;
    j_mu = {0};

    if ((tau_-tau < dtau) && (tau_-tau > 0)) {
        //std::stringstream buf;
    //    std::cout << "Dumping bullet into position x = " << x << ", y = "
    //        << y << ", eta_s = " << eta_s << std::endl;
        //music_message.info(buf.str());
        //music_message.flush("info");
        int  index[3] = {get_ix(x), get_iy(y), get_ieta(eta_s)};
        if (index[0] >= DATA_.nx) {
            music_message << "Out of bounds access in direction x";
            music_message.flush("info");
        }
        if (index[1] >= DATA_.ny) {
            music_message << "Out of bounds access in direction y";
            music_message.flush("info");
        }
        if (index[2] >= DATA_.neta) {
            music_message << "Out of bounds access in direction eta_s";
            music_message.flush("info");
        }
        //std::cout << index[0] << ", " << index[1] << ", " << index[2] << std::endl << std::flush;
        if( point_flagged(index[0],index[1],index[2]) ){
            //Compute point volume
            double dz = get_z(tau,eta_s+DATA_.delta_eta/2) - get_z(tau,eta_s-DATA_.delta_eta/2);
            double vol = DATA_.delta_x*DATA_.delta_y*dz;
            for (int idir=0; idir<4;++idir)
                j_mu[idir] = (pmu_frac_[idir]/vol/DATA_.delta_tau)/Util::hbarc;
        }
    }
    return;
}

//_________________________________________________________________________
//Algorithm to delimit the kick ellipse
GridT<int> HydroSourceBullet::TagPointsForFillEllipse(double x0, double y0, double eta0){

    //Stores which points on the grid will be filled
    GridT<int> points_flagged(DATA_.nx, DATA_.ny, DATA_.neta);

    //Get the index of the center points
    int ix0 = get_ix(x0);
    int iy0 = get_iy(y0);
    int ieta0 = get_ieta( eta0 );

    //Redefine x0, y0, z0 to match the grid
    x0 = get_x(ix0);
    y0 = get_y(iy0);
    eta0 = get_eta(ieta0);
    double z0 = get_z(tau_,get_eta(ieta0));

    //Get the index of the points which are at x0+sigma in each direction
    int ix_r_plus   = get_ix(x0+get_sigma_x());
    int iy_r_plus   = get_iy(y0+get_sigma_y());
    int ieta_r_plus = get_ieta(eta0+get_sigma_eta());
    double z_max = get_z(tau_,get_eta(ieta_r_plus));
    //Get the index of the points which are at x0-r in each direction
    int ix_r_minus   = get_ix(x0-get_sigma_x());
    int iy_r_minus   = get_iy(y0-get_sigma_x());
    int ieta_r_minus = get_ieta( eta0-get_sigma_eta() );
    double z_min = get_z(tau_,get_eta(ieta_r_minus));
    
    double rz = (z_max-z_min)/2.;

    //Perform loop tagging
    for (int ieta = ieta_r_minus; ieta <= ieta_r_plus; ++ieta){
        for (int iy = iy_r_minus; iy <= iy_r_plus; ++iy){
            for (int ix = ix_r_minus; ix <= ix_r_plus; ++ix){
                //Now we decide if we will fill it or not
                double x = get_x(ix);
                double y = get_y(iy);
                double z = get_z(tau_,  get_eta(ieta) );
                if ( pow((x-x0)/get_sigma_x(),2) + pow((y-y0)/get_sigma_x(),2)
                     + pow((z-z0)/rz,2) < 1.+1.e-7){
                    points_flagged(ix,iy,ieta) = 1;
                }
            }
        }
    }
    return points_flagged;
}

//_________________________________________________________________________
// Given a position x, computes the grid position ix
int HydroSourceBullet::get_ix(double x) const {
        return (int) std::round(  (DATA_.x_size/2 + x + 1e-6)/DATA_.delta_x ); }

//_________________________________________________________________________
// Given a position y, computes the grid position iy
int HydroSourceBullet::get_iy(double y) const {
        return (int) std::round(  (DATA_.y_size/2 + y + 1.e-6)/DATA_.delta_y ); }

//_________________________________________________________________________
// Given a position eta, computes the grid position ieta
int HydroSourceBullet::get_ieta(double eta) const{
        return (int) std::round(  (DATA_.eta_size/2 + eta + 1.e-6)/DATA_.delta_eta ); }

//_________________________________________________________________________
//Given a tau and z value, returns the spatial rapidity
double HydroSourceBullet::get_spatial_rapidity(double tau, double z) const {
    double t = sqrt( tau*tau + z*z );
    return .5*log( (t+z)/(t-z) );
}

//_________________________________________________________________________
//Get the position x of the site ix
double HydroSourceBullet::get_x (int ix) const {
        return ix*DATA_.delta_x-DATA_.x_size/2;}

//_________________________________________________________________________
//Get the position x of the site iy
double HydroSourceBullet::get_y (int iy) const {
        return iy*DATA_.delta_y-DATA_.y_size/2;}

//_________________________________________________________________________
//Get the position eta of the site ieta
double HydroSourceBullet::get_eta (int ieta) const {
        return ieta*DATA_.delta_eta-DATA_.eta_size/2;}

//_________________________________________________________________________
//Given a tau value and a a space-time rapidty, returns the z value
double HydroSourceBullet::get_z (double tau, double etaS) const{
    return tau*sinh(etaS);
}
