// Copyright 2018 @ Chun Shen

#ifndef SRC_EOS_4D_H_
#define SRC_EOS_4D_H_

#include "eos_base.h"

class EOS_4D : public EOS_base {
 private:
    // variables for header infos.
    double mubtilde0, muqtilde0, mustilde0, Ttilde0;
    double dmubtilde, dmuqtilde, dmustilde, dTtilde;
    int N_mub, N_muq, N_mus, N_T;
    //double Ttilde, mubtilde, muqtilde, mustilde;

    // useful constants 
    double pi;
    int Nf;
    double alphaNf;
    double OneoveralphaNf;


    // 1D tables.
    std::vector<double> pressure_vec;
    std::vector<double> temp_vec;
    std::vector<double> mub_vec;
    std::vector<double> muq_vec;
    std::vector<double> mus_vec;

    std::vector<double>* p_ = &pressure_vec;
    std::vector<double>* t_ = &temp_vec;
    std::vector<double>* mub_ = &mub_vec;
    std::vector<double>* muq_ = &muq_vec;
    std::vector<double>* mus_ = &mus_vec;



    // method to read/mainupalate header info and data
    std::vector<double> read_vector(std::string filepath, int header_size=2);
    void read_header(std::string filepath);
    int index(int i_T, int i_mub, int i_muq, int i_mus) const;
    //int Cfloor(double val);
    int check_index(int ind, int lim) const;
    std::vector<double> FourDLInterp(std::vector<double>* data, 
		    std::vector<double> TildeVar, bool compute_derivatives=false) const;
    std::vector<double> get_tilde_variables(double e, double nb, double nq, double ns) const;

 public:
    EOS_4D();
    ~EOS_4D();

    void initialize_eos();
    double get_temperature(double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muB        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muS        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_muQ        (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;
    double get_pressure   (double e, double rhob, double rhoq=0.0, double rhos=0.0) const;

    //void get_pressure_with_gradients(double epsilon, double rhob, double &p, 
	//	    double &dpde, double &dpdrhob, double &cs2, 
	//	    double rhoq=0.0, double rhos=0.0) const;

    void get_pressure_with_gradients(double epsilon, double rhob, double rhoq, double rhos, 
            double &p, 
		    double &dpde, double &dpdrhob, 
		    double &dpdrhoq, double &dpdrhos, 
		    double &cs2) const;

    void check_eos() const {
        check_eos_with_finite_muB();
        outputMutable();
    }
};

#endif  // SRC_EOS_4D_H_
