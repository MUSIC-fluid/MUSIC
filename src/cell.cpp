// Copyright 2012 Bjoern Schenke, Sangyong Jeon, and Charles Gale

#include "./cell.h"

Cell::Cell() {
    prev_epsilon = 0.0;
    epsilon_t    = 0.0;
    epsilon      = 0.0;

    prev_rhob = 0.0;
    rhob_t    = 0.0;
    rhob      = 0.0;
}
