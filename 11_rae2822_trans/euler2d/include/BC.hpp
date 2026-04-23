#pragma once
#include "Solver.hpp"

void bc_wall    (Solver& s);   // j=-1 ghost: reflect normal velocity
void bc_farfield(Solver& s);   // j=ncj ghost: Riemann characteristic
void bc_wake_cut(Solver& s);   // j=-1 wake ghosts: mirror across cut
void bc_i_exit  (Solver& s);   // i=-1, i=nci ghosts: zeroth-order extrapolation
