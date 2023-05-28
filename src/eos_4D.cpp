// Copyright 2018 @ Chun Shen

#include "eos_4D.h"
#include "util.h"

#include <sstream>
#include <fstream>
#include <cmath>

using std::stringstream;
using std::string;

EOS_4D::EOS_4D(){
    set_EOS_id(20);
    set_number_of_tables(0);
    resize_table_info_arrays();
    set_eps_max(1e5);
    set_flag_muB(true);
    set_flag_muS(false);
    set_flag_muQ(false);

    pi = 3.141592653589793;
    Nf = 3;
    alphaNf = (8/45.0 + 7/60.0*3)*3.141592653589793*3.141592653589793;
    OneoveralphaNf = 1/((8/45.0 + 7/60.0*3)*3.141592653589793*3.141592653589793);
}


EOS_4D::~EOS_4D() {}

void EOS_4D::read_header(std::string filepath){
    std::ifstream ifs(filepath);
    if (!ifs.is_open()) {
        music_message << "Can not open EOS files: "<< filepath;
        music_message.flush("error");
        exit(1);
    }
    ifs >> mubtilde0 >> muqtilde0 >> mustilde0 >> Ttilde0;
    ifs >> dmubtilde >> dmuqtilde >> dmustilde >> dTtilde;
    ifs >> N_mub >> N_muq >> N_mus >> N_T;

    mubtilde0 /= Util::hbarc;
    muqtilde0 /= Util::hbarc;
    mustilde0 /= Util::hbarc;
    Ttilde0 /= Util::hbarc;

    dmubtilde /= Util::hbarc;
    dmuqtilde /= Util::hbarc;
    dmustilde /= Util::hbarc;
    dTtilde /= Util::hbarc;

    N_T += 1;
    N_mub += 1;
    N_muq += 1;
    N_mus += 1;

    ifs.close();
}

std::vector<double> EOS_4D::read_vector(std::string filepath, int header_size){
    std::vector<double> out;
    std::ifstream eos(filepath);
    if (!eos.is_open()) {
        music_message << "Can not open EOS file: "<< filepath;
        music_message.flush("error");
        exit(1);
    }

    // skip header.
    std::string dummy;
    for(int i=0; i<header_size; i++){std::getline(eos, dummy);}

    // read file
    double dum1;
    std::string line;
        while(std::getline(eos, line)){
        for(int i=0; i<5;i++){
            std::stringstream ss(line);
            ss >> dum1;out.push_back(dum1/Util::hbarc);}
                }
        eos.close();
        out.resize(out.size());
    return out;
}

int EOS_4D::index(int i_T, int i_mub, int i_muq, int i_mus) const {
    //return i_T*N_mub*N_muq*N_mus + i_mub*N_muq*N_mus + i_muq*N_mus + i_mus;
    return i_T*N_mub*N_muq*N_mus + i_mus*N_muq*N_mub + i_muq*N_mub + i_mub;
}

int EOS_4D::Cfloor(double val) const {
    double eps = 0.00000000000001;
    int fl = std::floor(val);
    int cfl = std::floor(val+eps);
    if(cfl != fl){return cfl;}
    else{return fl;}
}

int EOS_4D::check_index(int ind, int lim) const {
    if(ind < 0){return 0;}
    else if(ind > lim-1){return lim-1;}
    else{return ind;}
}


std::vector<double> EOS_4D::FourDLInterp(std::vector<double>* data, std::vector<double> TildeVar, bool compute_derivatives) const {
    double T = TildeVar[0];
    double mub = TildeVar[1];
    double muq = TildeVar[2];
    double mus = TildeVar[3];

    if(T<0.0){T = 0.0;}
    if(T>5.07614213197969550){T = 5.0761421319796955;}

    if(mub<-1.2164216928535225){mub = -1.2164216928535225;}
    if(mub>3.0410542321338063){mub = 3.0410542321338063;}

    if(muq<-0.25342118601115055){muq = -0.25342118601115055;}
    if(muq>0.10136847440446022){muq = 0.10136847440446022;}

    if(mus<-0.5068423720223011){mus = -0.5068423720223011;}
    if(mus>1.2671059300557526){mus = 1.2671059300557526;}

    // Calculate the weights associated to the sixteen surrounding point 
    
    double indmuT = (T - Ttilde0)/dTtilde;
    double indmub = (mub - mubtilde0)/dmubtilde;
    double indmuq = (muq - muqtilde0)/dmuqtilde;
    double indmus = (mus - mustilde0)/dmustilde;

    //int iT = check_index(static_cast<int>(indmuT), N_T);int iT1 = check_index(iT + 1, N_T);
    //int ib = check_index(static_cast<int>(indmub), N_mub);int ib1 = check_index(ib + 1, N_mub);
    //int iq = check_index(static_cast<int>(indmuq), N_muq);int iq1 = check_index(iq + 1, N_muq);
    //int is = check_index(static_cast<int>(indmus), N_mus);int is1 = check_index(is + 1, N_mus);

    int iT = static_cast<int>(indmuT);int iT1 = iT + 1;
    int ib = static_cast<int>(indmub);int ib1 = ib + 1;
    int iq = static_cast<int>(indmuq);int iq1 = iq + 1;
    int is = static_cast<int>(indmus);int is1 = is + 1;
    //int iT = check_index(Cfloor(indmuT), N_T);int iT1 = check_index(iT + 1, N_T);
    //int ib = check_index(Cfloor(indmub), N_mub);int ib1 = check_index(ib + 1, N_mub);
    //int iq = check_index(Cfloor(indmuq), N_muq);int iq1 = check_index(iq + 1, N_muq);
    //int is = check_index(Cfloor(indmus), N_mus);int is1 = check_index(is + 1, N_mus);

    double dx = indmuT - (double)iT;
    double dy = indmub - (double)ib;
    double dz = indmuq - (double)iq;
    double dt = indmus - (double)is;

    double w0000 = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dt);
    double w1111 = dx * dy * dz * dt;

    double w1000 = dx * (1 - dy) * (1 - dz) * (1 - dt);
    double w0100 = (1 - dx) * dy * (1 - dz) * (1 - dt);
    double w0010 = (1 - dx) * (1 - dy) * dz * (1 - dt);
    double w0001 = (1 - dx) * (1 - dy) * (1 - dz) * dt;

    double w1001 = dx * (1 - dy) * (1 - dz) * dt;
    double w0101 = (1 - dx) * dy * (1 - dz) * dt;
    double w0011 = (1 - dx) * (1 - dy) * dz * dt;
    double w1100 = dx * dy * (1 - dz) * (1 - dt);
    double w1010 = dx * (1 - dy) * dz * (1 - dt);
    double w0110 = (1 - dx) * dy * dz * (1 - dt);

    double w0111 = (1 - dx) * dy * dz * dt;
    double w1011 = dx * (1 - dy) * dz * dt;
    double w1101 = dx * dy * (1 - dz) * dt;
    double w1110 = dx * dy * dz * (1 - dt);
    
    // store values at surrounding data points on the grid.
    
    double data_0000 = (*data)[index(iT,  ib, iq, is)];
    double data_1111 = (*data)[index(iT1,  ib1, iq1, is1)];

    double data_1000 = (*data)[index(iT1,  ib, iq, is)];
    double data_0100 = (*data)[index(iT,  ib1, iq, is)];
    double data_0010 = (*data)[index(iT,  ib, iq1, is)];
    double data_0001 = (*data)[index(iT,  ib, iq, is1)];

    double data_1001 = (*data)[index(iT1,  ib, iq, is1)];
    double data_0101 = (*data)[index(iT,  ib1, iq, is1)];
    double data_0011 = (*data)[index(iT,  ib, iq1, is1)];
    double data_1100 = (*data)[index(iT1,  ib1, iq, is)];
    double data_1010 = (*data)[index(iT1,  ib, iq1, is)];
    double data_0110 = (*data)[index(iT,  ib1, iq1, is)];

    double data_0111 = (*data)[index(iT,  ib1, iq1, is1)];
    double data_1011 = (*data)[index(iT1,  ib, iq1, is1)];
    double data_1101 = (*data)[index(iT1,  ib1, iq, is1)];
    double data_1110 = (*data)[index(iT1,  ib1, iq1, is)];

    // Interpolate the value of the target point using the weights and the values of the sixteen surrounding points
    double interpolated_value = w0000 * data_0000 + w1111 * data_1111  + w1000 * data_1000 
        + w0100 * data_0100 + w0010 * data_0010 + w0001 * data_0001 
        + w1001 * data_1001 + w0101 * data_0101 + w0011 * data_0011 
        + w1100 * data_1100 + w1010 * data_1010 + w0110 * data_0110 
        + w0111 * data_0111 + w1011 * data_1011 + w1101 * data_1101 
        + w1110 * data_1110;

    double dXoverde = 0.0; 
    double dXoverdrhob = 0.0;
    double dXoverdrhoq = 0.0;
    double dXoverdrhos = 0.0;
    if(compute_derivatives){
        // Calculate derivatives.
        
        // Ttilde direction
        double wT000 = (1-dy)*(1-dz)*(1-dt);double wT111 = dy*dz*dt;
        double wT100 = dy*(1-dz)*(1-dt);double wT110 = dy*dz*(1-dt);
        double wT010 = (1-dy)*dz*(1-dt);double wT011 = (1-dy)*dz*dt;
        double wT001 = (1-dy)*(1-dz)*dt;double wT101 = dy*(1-dz)*dt;

        double tempT1 = wT000 * data_0000 + wT100 * data_0100 + wT010 * data_0010 + wT001 * data_0001 
                + wT101 * data_0101 + wT011 * data_0011 + wT110 * data_0110 + wT111 * data_0111; 

        double tempT2 = wT111 * data_1111 + wT000 * data_1000 + wT001 * data_1001 + wT100 * data_1100 
                + wT010 * data_1010 + wT011 * data_1011 + wT101 * data_1101 + wT110 * data_1110;

        double dXdTtilde = (tempT2 - tempT1)/dTtilde; 
        // to test for speed
        //double dXdTtilde = (tempT2 - interpolated_value)/((1-dx) * dT);

        // mubtilde direction
        double wb000 = (1-dx)*(1-dz)*(1-dt);double wb111 = dx*dz*dt;
        double wb100 = dx*(1-dz)*(1-dt);double wb110 = dx*dz*(1-dt);
        double wb010 = (1-dx)*dz*(1-dt);double wb011 = (1-dx)*dz*dt;
        double wb001 = (1-dx)*(1-dz)*dt;double wb101 = dx*(1-dz)*dt;

        double tempb1 = wb000 * data_0000 + wb100 * data_1000 + wb010 * data_0010 + wb001 * data_0001 
                + wb101 * data_1001 + wb011 * data_0011 + wb110 * data_1010 + wb111 * data_1011; 

        double tempb2 =  wb111 * data_1111 + wb000 * data_0100 + wb001 * data_0101 + wb100 * data_1100 
                + wb010 * data_0110 + wb011 * data_0111 + wb101 * data_1101 + wb110 * data_1110;

        double dXdmubtilde = (tempb2 - tempb1)/dmubtilde; 
        // to test for speed
        //double dXdmubtilde = (tempb2 - interpolated_value)/((1-dy) * dmub);

        // muqtilde direction
        double wq000 = (1-dx)*(1-dy)*(1-dt);double wq111 = dx*dy*dt;
        double wq100 = dx*(1-dy)*(1-dt);double wq110 = dx*dy*(1-dt);
        double wq010 = (1-dx)*dy*(1-dt);double wq011 = (1-dx)*dy*dt;
        double wq001 = (1-dx)*(1-dy)*dt;double wq101 = dx*(1-dy)*dt;

        double tempq1 = wq000 * data_0000 + wq100 * data_1000 + wq010 * data_0100 + wq001 * data_0001 
                + wq101 * data_1001 + wq011 * data_0101 + wq110 * data_1100 + wq111 * data_1101; 

        double tempq2 =  wq111 * data_1111 + wq000 * data_0010 + wq001 * data_0011 + wq100 * data_1010 
                + wq010 * data_0110 + wq011 * data_0111 + wq101 * data_1011 + wq110 * data_1110;

        double dXdmuqtilde = (tempq2 - tempq1)/dmuqtilde; 
        // to test for speed
        //double dXdmuqtilde = (tempq2 - interpolated_value)/((1-dz) * dmuq);

        // mustilde direction
        double ws000 = (1-dx)*(1-dy)*(1-dz);double ws111 = dx*dy*dz;
        double ws100 = dx*(1-dy)*(1-dz);double ws110 = dx*dy*(1-dz);
        double ws010 = (1-dx)*dy*(1-dz);double ws011 = (1-dx)*dy*dz;
        double ws001 = (1-dx)*(1-dy)*dz;double ws101 = dx*(1-dy)*dz;

        double temps1 = ws000 * data_0000 + ws100 * data_1000 + ws010 * data_0100 + ws001 * data_0010 
                + ws110 * data_1100 + ws101 * data_1010 + ws011 * data_0110 + ws111 * data_1110;

        double temps2 =  ws111 * data_1111 + ws000 * data_0001 + ws100 * data_1001 + ws010 * data_0101 
                + ws001 * data_0011 + ws011 * data_0111 + ws101 * data_1011 + ws110 * data_1101; 

        double dXdmustilde = (temps2 - temps1)/dmustilde; 
        // to test for speed
        //double dXdmustilde = (temps2 - interpolated_value)/((1-dt) * dmus);

        dXoverde = 1/(12 * alphaNf * T * T * T) * dXdTtilde; 
        dXoverdrhob = (5 * dXdmubtilde - dXdmuqtilde + 2 * dXdmustilde)/(T * T);
        dXoverdrhoq = (- 1 * dXdmubtilde + 2 * dXdmuqtilde - dXdmustilde)/(T * T);
        dXoverdrhos = (2 * dXdmubtilde - dXdmuqtilde + 2 * dXdmustilde)/(T * T);
    }
    
    std::vector<double> out;
    out.push_back(interpolated_value);
    out.push_back(dXoverde);
    out.push_back(dXoverdrhob);
    out.push_back(dXoverdrhoq);
    out.push_back(dXoverdrhos);
    out.resize(out.size());
    return out;

}

std::vector<double> EOS_4D::get_tilde_variables(double e, double rhob, double rhoq, double rhos) const {
    double Ttilde = sqrt(sqrt(e/3.0 * OneoveralphaNf)); // fm-1 
    double mubtilde = (5.0 * rhob - rhoq + 2.0*rhos)/(Ttilde*Ttilde); // fm-1 
    double muqtilde = (2.0 * rhoq - rhob -   rhos)/(Ttilde*Ttilde); // fm-1 
    double mustilde = (2.0 * rhob - rhoq + 2.0*rhos)/(Ttilde*Ttilde); // fm-1 
                                                                      //
    
    //std::cout << mubtilde << " nb " << rhob << " nq " << rhoq << " ns " << rhos << " T " << Ttilde << std::endl;

    // convert to GeV for 4D EoS tables. 
    
    //Ttilde *= Util::hbarc; // GeV
    //mubtilde *= Util::hbarc; // GeV
    //muqtilde *= Util::hbarc; // GeV
    //mustilde *= Util::hbarc; // GeV

    std::vector<double> out;
    out.push_back(Ttilde);
    out.push_back(mubtilde);
    out.push_back(muqtilde);
    out.push_back(mustilde);
    return out;
} 

void EOS_4D::initialize_eos() {
    music_message.info("Using 4D EOS");

    auto envPath = get_hydro_env_path();
    std::stringstream spath;
    spath << envPath;
    spath << "/EOS/neos4D/";
    std::string path = spath.str();
    music_message << "from path " << path;
    music_message.flush("info");

    // read header info
    read_header(path + "neos4d_t.dat");

    // read vectors
    pressure_vec = read_vector(path + "neos4d_p.dat");
    temp_vec = read_vector(path + "neos4d_t.dat");
    mub_vec = read_vector(path + "neos4d_mub.dat");
    muq_vec = read_vector(path + "neos4d_muq.dat");
    mus_vec = read_vector(path + "neos4d_mus.dat");
    music_message.info("Done reading EOS.");
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_4D::get_temperature(double e, double rhob, double rhoq, double rhos) const {
    std::vector<double> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    std::vector<double> interp_output = FourDLInterp(t_, TildeVar);  // 1/fm^5
    double interp_T = std::max(Util::small_eps, interp_output[0]);
    if(interp_T > 9){std::cout << interp_T << std::endl;}
    return(interp_T);
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_4D::get_pressure(double e, double rhob, double rhoq, double rhos) const {
    std::vector<double> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    std::vector<double> interp_output = FourDLInterp(p_, TildeVar);  
    double interp_p = std::max(Util::small_eps, interp_output[0]);
    return(interp_p);
}

void EOS_4D::get_pressure_with_gradients(double e, double rhob, double rhoq, double rhos, double &p, double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos, double &cs2) const {

    std::vector<double> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    //std::vector<double> interp_output = FourDLInterp(pressure_vec, TildeVar, true);  // 1/fm^5
    std::vector<double> interp_output = FourDLInterp(p_, TildeVar, true);  // 1/fm^5
    p = std::max(Util::small_eps, interp_output[0]);

    dpde = interp_output[1]; 
    dpdrhob = interp_output[2]; 
    dpdrhoq = interp_output[3]; 
    dpdrhos = interp_output[4]; 

    cs2 = dpde + rhob/(e + p + Util::small_eps)*dpdrhob + 
        rhoq/(e + p + Util::small_eps)*dpdrhoq + 
        rhos/(e + p + Util::small_eps)*dpdrhos;
    cs2 = std::max(0.01, std::min(1./3, cs2));
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_4D::get_muB(double e, double rhob, double rhoq, double rhos) const {
    std::vector<double> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    //std::vector<double> interp_output = FourDLInterp(mub_vec, TildeVar);  // 1/fm^5
    std::vector<double> interp_output = FourDLInterp(mub_, TildeVar);  // 1/fm^5
    return(interp_output[0]);
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_4D::get_muS(double e, double rhob, double rhoq, double rhos) const {
    std::vector<double> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    //std::vector<double> interp_output = FourDLInterp(mus_vec, TildeVar);  // 1/fm^5
    std::vector<double> interp_output = FourDLInterp(mus_, TildeVar);  // 1/fm^5
    return(interp_output[0]);
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_4D::get_muQ(double e, double rhob, double rhoq, double rhos) const {
    std::vector<double> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    //std::vector<double> interp_output = FourDLInterp(muq_vec, TildeVar);  // 1/fm^5
    std::vector<double> interp_output = FourDLInterp(muq_, TildeVar);  // 1/fm^5
    return(interp_output[0]);
}
