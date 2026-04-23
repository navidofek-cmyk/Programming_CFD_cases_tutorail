#pragma once

// Roe numerical flux with Harten entropy fix.
// Inputs: left/right primitive states, face unit normal (nx,ny), gamma.
// Output: F[4] — physical normal flux (F·n scaled by face length NOT included).
void roe_flux(double rL, double uL, double vL, double pL,
              double rR, double uR, double vR, double pR,
              double nx, double ny, double gamma,
              double* F);
