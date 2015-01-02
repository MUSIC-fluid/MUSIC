#ifndef MINMOD_H
#define MINMOD_H

#include "data.h"

class Minmod
{
    private:
        double theta_flux;

    public:
        Minmod(InitData* DATA);
        ~Minmod();
        double minmod_dx(double up1, double u, double um1);
        double minmod_theta_dx(double up1, double u, double um1, double theta);
};
#endif
