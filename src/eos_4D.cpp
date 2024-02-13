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
    alphaNf = (8/45.0 + 7/60.0*3.0)*3.141592653589793*3.141592653589793;
    OneoveralphaNf = 1/((8/45.0 + 7/60.0*3.0)*3.141592653589793*3.141592653589793);
}


EOS_4D::~EOS_4D() {}

void EOS_4D::get_eos_max_values(){
    T_tilde_max = Ttilde0 + (N_T - 1)*dTtilde;
    mub_tilde_max = mubtilde0 + (N_mub - 1)*dmubtilde;
    muq_tilde_max = muqtilde0 + (N_muq - 1)*dmuqtilde;
    mus_tilde_max = mustilde0 + (N_mus - 1)*dmustilde;
}

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

    // Get variables in fm-1 divide by hbarc
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


void EOS_4D::read_header_binary(std::string filepath, int header_size){
    std::ifstream ifs(filepath, std::ios::in | std::ios::binary);
    if (!ifs.is_open()) {
        music_message << "Can not open EOS files: "<< filepath;
        music_message.flush("error");
        exit(1);
    }
    std::vector<float> hd(header_size);
    ifs.read(reinterpret_cast<char*>(hd.data()), sizeof(float) * header_size);

    Ttilde0  = hd[3];dTtilde = hd[7];
    mubtilde0 = hd[0];muqtilde0 = hd[1];mustilde0 = hd[2];
    dmubtilde = hd[4];dmuqtilde = hd[5];dmustilde = hd[6];
    N_mub = hd[8];N_muq = hd[9];N_mus = hd[10];N_T = hd[11];


    // Get variables in fm-1 divide by hbarc
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
    std::cout << N_T << " " << N_mus <<std::endl;

}



std::vector<float> EOS_4D::read_eos(std::string filepath, bool is_cs, int header_size){
    std::vector<float> out;
    std::ifstream eos(filepath);

    // EoS is in GeV, MUSIC in fm-1
    // cs is dimensionless
    float dimension = Util::hbarc;
    if(is_cs){dimension = 1.0;}


    if (!eos.is_open()) {
        music_message << "Can not open EOS file: "<< filepath;
        music_message.flush("error");
        exit(1);
    }

    // skip header.
    std::string dummy;
    for(int i = 0; i < header_size; i++) {
        std::getline(eos, dummy);
    }

    // read file
    float dum1;
    std::string line;
    while (std::getline(eos, line)) {
        std::stringstream ss(line);
        for(int i = 0; i < 5; i++) {
            ss >> dum1;
            out.push_back(dum1/dimension);
        }
    }
    eos.close();
    return out;
}

std::vector<float> EOS_4D::read_eos_binary(std::string filepath, bool is_cs, int header_size){
    std::ifstream eos_binary_file(filepath, std::ios::in | std::ios::binary);
    std::vector<float> out;

    // EoS is in GeV, MUSIC in fm-1
    // cs is dimensionless
    float dimension = Util::hbarc;
    if(is_cs){dimension = 1.0;}

    if (!eos_binary_file.is_open()) {
        music_message << "Can not open EOS file: "<< filepath;
        music_message.flush("error");
        exit(1);
    }
    float number; 
    while (eos_binary_file.read(reinterpret_cast<char*>(&number), sizeof(float))) {
        out.push_back(number/dimension);
    }
    eos_binary_file.close();

    // remove header
    out.erase(out.begin(), out.begin() + header_size);
    out.resize(out.size());
    return out;
}

int EOS_4D::index(int i_T, int i_mub, int i_muq, int i_mus) const {
    int idx = ((i_T*N_mus + i_mus)*N_muq + i_muq)*N_mub + i_mub;
    return(idx);
}

std::vector<float> EOS_4D::FourDLInterp(const std::vector<float> &data,
                                         const std::vector<float> &TildeVar,
                                         bool compute_derivatives) const {
    float T = TildeVar.at(0);
    float mub = TildeVar.at(1);
    float muq = TildeVar.at(2);
    float mus = TildeVar.at(3);

    // Constrain the input tilde variables values to the table boundaries.
    
    if(T<Ttilde0){T = Ttilde0;}
    if(T>T_tilde_max){T = T_tilde_max;}

    if(mub<mubtilde0){mub = mubtilde0;}
    if(mub>mub_tilde_max){mub = mub_tilde_max;}

    if(muq<muqtilde0){muq = muqtilde0;}
    if(muq>muq_tilde_max){muq = muq_tilde_max;}

    if(mus<mustilde0){mus = mustilde0;}
    if(mus>mus_tilde_max){mus =mus_tilde_max;}

    // Calculate the weights associated to the sixteen surrounding point 

    float indmuT = (T - Ttilde0)/dTtilde;
    float indmub = (mub - mubtilde0)/dmubtilde;
    float indmuq = (muq - muqtilde0)/dmuqtilde;
    float indmus = (mus - mustilde0)/dmustilde;

    int iT = static_cast<int>(indmuT);int iT1 = iT + 1;
    int ib = static_cast<int>(indmub);int ib1 = ib + 1;
    int iq = static_cast<int>(indmuq);int iq1 = iq + 1;
    int is = static_cast<int>(indmus);int is1 = is + 1;

    float dx = indmuT - iT;
    float dy = indmub - ib;
    float dz = indmuq - iq;
    float dt = indmus - is;

    float w0000 = (1 - dx) * (1 - dy) * (1 - dz) * (1 - dt);
    float w1111 = dx * dy * dz * dt;

    float w1000 = dx * (1 - dy) * (1 - dz) * (1 - dt);
    float w0100 = (1 - dx) * dy * (1 - dz) * (1 - dt);
    float w0010 = (1 - dx) * (1 - dy) * dz * (1 - dt);
    float w0001 = (1 - dx) * (1 - dy) * (1 - dz) * dt;

    float w1001 = dx * (1 - dy) * (1 - dz) * dt;
    float w0101 = (1 - dx) * dy * (1 - dz) * dt;
    float w0011 = (1 - dx) * (1 - dy) * dz * dt;
    float w1100 = dx * dy * (1 - dz) * (1 - dt);
    float w1010 = dx * (1 - dy) * dz * (1 - dt);
    float w0110 = (1 - dx) * dy * dz * (1 - dt);

    float w0111 = (1 - dx) * dy * dz * dt;
    float w1011 = dx * (1 - dy) * dz * dt;
    float w1101 = dx * dy * (1 - dz) * dt;
    float w1110 = dx * dy * dz * (1 - dt);

    // store values at surrounding data points on the grid.

    float data_0000 = data.at(index(iT,  ib, iq, is));
    float data_1111 = data.at(index(iT1,  ib1, iq1, is1));

    float data_1000 = data.at(index(iT1,  ib, iq, is));
    float data_0100 = data.at(index(iT,  ib1, iq, is));
    float data_0010 = data.at(index(iT,  ib, iq1, is));
    float data_0001 = data.at(index(iT,  ib, iq, is1));

    float data_1001 = data.at(index(iT1,  ib, iq, is1));
    float data_0101 = data.at(index(iT,  ib1, iq, is1));
    float data_0011 = data.at(index(iT,  ib, iq1, is1));
    float data_1100 = data.at(index(iT1,  ib1, iq, is));
    float data_1010 = data.at(index(iT1,  ib, iq1, is));
    float data_0110 = data.at(index(iT,  ib1, iq1, is));

    float data_0111 = data.at(index(iT,  ib1, iq1, is1));
    float data_1011 = data.at(index(iT1,  ib, iq1, is1));
    float data_1101 = data.at(index(iT1,  ib1, iq, is1));
    float data_1110 = data.at(index(iT1,  ib1, iq1, is));

    // Interpolate the value of the target point using the weights and the values of the sixteen surrounding points
    float interpolated_value = w0000 * data_0000 + w1111 * data_1111  + w1000 * data_1000 
        + w0100 * data_0100 + w0010 * data_0010 + w0001 * data_0001 
        + w1001 * data_1001 + w0101 * data_0101 + w0011 * data_0011 
        + w1100 * data_1100 + w1010 * data_1010 + w0110 * data_0110 
        + w0111 * data_0111 + w1011 * data_1011 + w1101 * data_1101 
        + w1110 * data_1110;

    float dXoverde = 0.0; 
    float dXoverdrhob = 0.0;
    float dXoverdrhoq = 0.0;
    float dXoverdrhos = 0.0;
    if(compute_derivatives){
        // Calculate derivatives.

        // Ttilde direction
        float wT000 = (1-dy)*(1-dz)*(1-dt);float wT111 = dy*dz*dt;
        float wT100 = dy*(1-dz)*(1-dt);float wT110 = dy*dz*(1-dt);
        float wT010 = (1-dy)*dz*(1-dt);float wT011 = (1-dy)*dz*dt;
        float wT001 = (1-dy)*(1-dz)*dt;float wT101 = dy*(1-dz)*dt;

        float tempT1 = wT000 * data_0000 + wT100 * data_0100 + wT010 * data_0010 + wT001 * data_0001 
                + wT101 * data_0101 + wT011 * data_0011 + wT110 * data_0110 + wT111 * data_0111; 

        float tempT2 = wT111 * data_1111 + wT000 * data_1000 + wT001 * data_1001 + wT100 * data_1100 
                + wT010 * data_1010 + wT011 * data_1011 + wT101 * data_1101 + wT110 * data_1110;

        float dXdTtilde = (tempT2 - tempT1)/dTtilde; 
        // to test for speed
        //float dXdTtilde = (tempT2 - interpolated_value)/((1-dx) * dTtilde);

        // mubtilde direction
        float wb000 = (1-dx)*(1-dz)*(1-dt);float wb111 = dx*dz*dt;
        float wb100 = dx*(1-dz)*(1-dt);float wb110 = dx*dz*(1-dt);
        float wb010 = (1-dx)*dz*(1-dt);float wb011 = (1-dx)*dz*dt;
        float wb001 = (1-dx)*(1-dz)*dt;float wb101 = dx*(1-dz)*dt;

        float tempb1 = wb000 * data_0000 + wb100 * data_1000 + wb010 * data_0010 + wb001 * data_0001 
                + wb101 * data_1001 + wb011 * data_0011 + wb110 * data_1010 + wb111 * data_1011; 

        float tempb2 =  wb111 * data_1111 + wb000 * data_0100 + wb001 * data_0101 + wb100 * data_1100 
                + wb010 * data_0110 + wb011 * data_0111 + wb101 * data_1101 + wb110 * data_1110;

        float dXdmubtilde = (tempb2 - tempb1)/dmubtilde; 
        // to test for speed
        //float dXdmubtilde = (tempb2 - interpolated_value)/((1-dy) * dmubtilde);

        // muqtilde direction
        float wq000 = (1-dx)*(1-dy)*(1-dt);float wq111 = dx*dy*dt;
        float wq100 = dx*(1-dy)*(1-dt);float wq110 = dx*dy*(1-dt);
        float wq010 = (1-dx)*dy*(1-dt);float wq011 = (1-dx)*dy*dt;
        float wq001 = (1-dx)*(1-dy)*dt;float wq101 = dx*(1-dy)*dt;

        float tempq1 = wq000 * data_0000 + wq100 * data_1000 + wq010 * data_0100 + wq001 * data_0001 
                + wq101 * data_1001 + wq011 * data_0101 + wq110 * data_1100 + wq111 * data_1101; 

        float tempq2 =  wq111 * data_1111 + wq000 * data_0010 + wq001 * data_0011 + wq100 * data_1010 
                + wq010 * data_0110 + wq011 * data_0111 + wq101 * data_1011 + wq110 * data_1110;

        float dXdmuqtilde = (tempq2 - tempq1)/dmuqtilde; 
        // to test for speed
        //float dXdmuqtilde = (tempq2 - interpolated_value)/((1-dz) * dmuq);

        // mustilde direction
        float ws000 = (1-dx)*(1-dy)*(1-dz);float ws111 = dx*dy*dz;
        float ws100 = dx*(1-dy)*(1-dz);float ws110 = dx*dy*(1-dz);
        float ws010 = (1-dx)*dy*(1-dz);float ws011 = (1-dx)*dy*dz;
        float ws001 = (1-dx)*(1-dy)*dz;float ws101 = dx*(1-dy)*dz;

        float temps1 = ws000 * data_0000 + ws100 * data_1000 + ws010 * data_0100 + ws001 * data_0010 
                + ws110 * data_1100 + ws101 * data_1010 + ws011 * data_0110 + ws111 * data_1110;

        float temps2 =  ws111 * data_1111 + ws000 * data_0001 + ws100 * data_1001 + ws010 * data_0101 
                + ws001 * data_0011 + ws011 * data_0111 + ws101 * data_1011 + ws110 * data_1101; 

        float dXdmustilde = (temps2 - temps1)/dmustilde; 
        // to test for speed
        //float dXdmustilde = (temps2 - interpolated_value)/((1-dt) * dmus);

        dXoverde = 3.0/(19.0 * pi *  pi * T * T * T) * dXdTtilde; 
        dXoverdrhob = (5.0 * dXdmubtilde - dXdmuqtilde + 2.0 * dXdmustilde)/(T * T);
        dXoverdrhoq = (- 1.0 * dXdmubtilde + 2.0 * dXdmuqtilde - dXdmustilde)/(T * T);
        dXoverdrhos = (2.0 * dXdmubtilde - dXdmuqtilde + 2.0 * dXdmustilde)/(T * T);
    }

    std::vector<float> out;
    out.push_back(interpolated_value);
    out.push_back(dXoverde);
    out.push_back(dXoverdrhob);
    out.push_back(dXoverdrhoq);
    out.push_back(dXoverdrhos);
    return out;
}


std::vector<float> EOS_4D::get_tilde_variables(double e, double rhob, double rhoq, double rhos) const {
    double Ttilde = sqrt(sqrt(e/3.0 * OneoveralphaNf)); // fm-1
    double mubtilde = (5.0 * rhob - rhoq + 2.0*rhos)/(Ttilde*Ttilde); // fm-1
    double muqtilde = (2.0 * rhoq - rhob -   rhos)/(Ttilde*Ttilde); // fm-1
    double mustilde = (2.0 * rhob - rhoq + 2.0*rhos)/(Ttilde*Ttilde); // fm-1

    std::vector<float> out;
    out.push_back((float) Ttilde);
    out.push_back((float) mubtilde);
    out.push_back((float) muqtilde);
    out.push_back((float) mustilde);
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
    
    /////////////////////////////////////////////////////////////////////////////
    // Set the type of EoS. Should be connected to parameters.
    set_eos_file_in_binary(true);
    // Set if the cs2 file is here.
    file_for_cs = true;
    /////////////////////////////////////////////////////////////////////////////

    if(EoS_file_in_binary){
        // Read EoS in binary
        // 1D Tables
        pressure_vec = read_eos_binary(path + "neos4d_urqmd_p_b.dat", false);
        temp_vec = read_eos_binary(path + "neos4d_urqmd_t_b.dat", false);
        mub_vec = read_eos_binary(path  + "neos4d_urqmd_mub_b.dat", false);
        muq_vec = read_eos_binary(path  + "neos4d_urqmd_muq_b.dat", false);
        mus_vec = read_eos_binary(path  + "neos4d_urqmd_mus_b.dat", false);
        if(file_for_cs){
            cs_vec = read_eos_binary(path + "neos4d_urqmd_cs_b.dat", true);
        }

        // Header info 
        read_header_binary(path + "neos4d_urqmd_t_b.dat");
    }
    else{
        // read 4D eos in txt format. 
        // read header info
        read_header(path + "neos4d_t.dat");

        // read 1D data into 1D tables
        pressure_vec = read_eos(path + "neos4d_p.dat", false);
        temp_vec = read_eos(path + "neos4d_t.dat", false);
        mub_vec = read_eos(path + "neos4d_mub.dat", false);
        muq_vec = read_eos(path + "neos4d_muq.dat", false);
        mus_vec = read_eos(path + "neos4d_mus.dat", false);
        if(file_for_cs){
            cs_vec = read_eos(path + "neos4d_cs.dat", true);
        }
    }

    // Get max values of EoS table for interpolation boundary.
    get_eos_max_values();
    music_message.info("Done reading EOS.");
}


//! This function returns the local temperature in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_4D::get_temperature(double e, double rhob, double rhoq, double rhos) const {
    std::vector<float> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    std::vector<float> interp_output = FourDLInterp(temp_vec, TildeVar);  // 1/fm
    double interp_T = std::max(Util::small_eps, (double) interp_output[0]);
    return(interp_T);
}


//! This function returns the local pressure in [1/fm^4]
//! the input local energy density [1/fm^4], rhob [1/fm^3]
double EOS_4D::get_pressure(double e, double rhob, double rhoq, double rhos) const {
    std::vector<float> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    std::vector<float> interp_output = FourDLInterp(pressure_vec, TildeVar);
    double interp_p = std::max(Util::small_eps, (double) interp_output[0]);
    return(interp_p);
}

//! This function returns the speed of sound. 
//! the input local energy density [1/fm^4], rhob [1/fm^3], rhoq [1/fm^3], rhos [1/fm^3]
double EOS_4D::get_cs2(double e, double rhob, double rhoq, double rhos) const {
    std::vector<float> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    double cs2;
    double interp_cs;
    if(file_for_cs){
        // read cs from file
        std::vector<float> interp_output = FourDLInterp(cs_vec, TildeVar);
        interp_cs = std::max(Util::small_eps, (double) interp_output[0]);
        cs2 = interp_cs*interp_cs;
    }
    else{
        // compute cs2 
        std::vector<float> interp_output = FourDLInterp(pressure_vec, TildeVar, true);
        float p = std::max(Util::small_eps, (double) interp_output[0]);
        float dpde = interp_output[1];
        float dpdrhob = interp_output[2];
        float dpdrhoq = interp_output[3];
        float dpdrhos = interp_output[4];
        cs2 = dpde + rhob/(e + p + Util::small_eps)*dpdrhob +
            rhoq/(e + p + Util::small_eps)*dpdrhoq +
            rhos/(e + p + Util::small_eps)*dpdrhos;
    }
    cs2 = std::max(0.01, std::min(1./3, cs2));
    return cs2;

}

void EOS_4D::get_pressure_with_gradients(double e, double rhob, double rhoq, double rhos, double &p, double &dpde, double &dpdrhob, double &dpdrhoq, double &dpdrhos, double &cs2) const {

    std::vector<float> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    std::vector<float> interp_output = FourDLInterp(pressure_vec, TildeVar, true);
    p = std::max(Util::small_eps, (double) interp_output[0]);

    dpde = (double) interp_output[1];
    dpdrhob = (double) interp_output[2];
    dpdrhoq = (double) interp_output[3];
    dpdrhos = (double) interp_output[4];

    if(file_for_cs){
        std::vector<float> interp_cs = FourDLInterp(cs_vec, TildeVar);
        cs2 = (double) interp_cs[0]*interp_cs[0];
    }
    else{

        cs2 = dpde + rhob/(e + p + Util::small_eps)*dpdrhob +
            rhoq/(e + p + Util::small_eps)*dpdrhoq +
            rhos/(e + p + Util::small_eps)*dpdrhos;
    }
    cs2 = std::max(0.01, std::min(1./3, cs2));
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_4D::get_muB(double e, double rhob, double rhoq, double rhos) const {
    std::vector<float> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    std::vector<float> interp_output = FourDLInterp(mub_vec, TildeVar);  // 1/fm
    return((double) interp_output[0]);
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_4D::get_muS(double e, double rhob, double rhoq, double rhos) const {
    std::vector<float> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    std::vector<float> interp_output = FourDLInterp(mus_vec, TildeVar);  // 1/fm
    return((double) interp_output[0]);
}


//! This function returns the local baryon chemical potential  mu_B in [1/fm]
//! input local energy density eps [1/fm^4] and rhob [1/fm^3]
double EOS_4D::get_muQ(double e, double rhob, double rhoq, double rhos) const {
    std::vector<float> TildeVar = get_tilde_variables(e, rhob, rhoq, rhos);
    std::vector<float> interp_output = FourDLInterp(muq_vec, TildeVar);  // 1/fm
    return((double) interp_output[0]);
}
