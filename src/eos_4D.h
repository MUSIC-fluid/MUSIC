// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_4D_H_
#define SRC_EOS_4D_H_

#include "eos_base.h"

class EOS_4D : public EOS_base {
 private:
    // variables for header infos.
    float mubtilde0, muqtilde0, mustilde0, Ttilde0;
    float dmubtilde, dmuqtilde, dmustilde, dTtilde;
    int N_mub, N_muq, N_mus, N_T;

    // max values for the tilde variables
    float T_tilde_max;
    float mub_tilde_max;
    float muq_tilde_max;
    float mus_tilde_max;
    //double Ttilde, mubtilde, muqtilde, mustilde;

    // useful constants 
    double pi;
    int Nf;
    double alphaNf;
    double OneoveralphaNf;

    bool EoS_file_in_binary;
    bool file_for_cs;


    // 1D tables.
    std::vector<float> pressure_vec;
    std::vector<float> temp_vec;
    std::vector<float> mub_vec;
    std::vector<float> muq_vec;
    std::vector<float> mus_vec;
    std::vector<float> cs_vec;

    std::vector<float>* p_ = &pressure_vec;
    std::vector<float>* t_ = &temp_vec;
    std::vector<float>* mub_ = &mub_vec;
    std::vector<float>* muq_ = &muq_vec;
    std::vector<float>* mus_ = &mus_vec;
    std::vector<float>* cs_ = &cs_vec;



    // method to read/mainupalate header info and data
    std::vector<float> read_eos(std::string filepath, bool is_cs, int header_size=2);
    std::vector<float> read_eos_binary(std::string filepath, bool is_cs, int header_size=12);

    void get_eos_max_values();

    void read_header(std::string filepath);
    void read_header_binary(std::string filepath, int header_size=12);

    // Shift in the index corresponds to the header size.
    int index(int i_T, int i_mub, int i_muq, int i_mus) const;

    std::vector<float> FourDLInterp(
        const std::vector<float> &data, const std::vector<float> &TildeVar,
        bool compute_derivatives=false) const;
    std::vector<float> get_tilde_variables(
            double e, double nb, double nq, double ns) const;

 public:
    EOS_4D();
    ~EOS_4D();

    void initialize_eos();

    void set_eos_file_in_binary(bool is_binary){EoS_file_in_binary=is_binary;}
    bool get_eos_file_in_binary(){return EoS_file_in_binary;}

    double get_temperature(double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muB        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muS        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muQ        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_pressure   (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_cs2   (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;

    void get_pressure_with_gradients(
            double epsilon, double rhob, double rhoq, double rhos,
            double &p, double &dpde, double &dpdrhob,
            double &dpdrhoq, double &dpdrhos, double &cs2) const;

    void check_eos() const {
        check_4D_eos();
    }
};

#endif  // SRC_EOS_4D_H_
