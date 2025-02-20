#ifndef SURFACECELL_H
#define SURFACECELL_H

class SurfaceCell {
 public:
    float xmu[4];               //!< Surface position [fm]
    float d3sigma_mu[4];        //!< Surface vector
    float energy_density;       //!< Local energy density [GeV/fm^3]
    float temperature;          //!< Local temperature [GeV]
    float pressure;             //!< Thermal pressure [GeV/fm^3]
    float umu[4];               //!< Flow velocity
    float rho_b;                //!< baryon density [1/fm^3]
    float rho_q;                //!< electric charge density [1/fm^3]
    float rho_s;                //!< strangeness density [1/fm^3]
    float mu_B;                 //!< Net baryon chemical potential [GeV]
    float mu_Q;                 //!< Net charge chemical potential [GeV]
    float mu_S;                 //!< Net strangeness chemical potential [GeV]
    float shear_pi[10];         //!< Shear stress tensor [GeV/fm^3]
    float bulk_Pi;              //!< Bulk viscous pressure [GeV/fm^3]

    SurfaceCell() = default;
    ~SurfaceCell(){};
};

#endif  // SURFACECELL_H
