#include "minmod.h"

Minmod::Minmod(const InitData &DATA) : theta_flux(DATA.minmod_theta) {}
Minmod::Minmod(double theta_in) : theta_flux(theta_in) {}
