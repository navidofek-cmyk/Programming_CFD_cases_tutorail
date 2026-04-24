#include "Flux.hpp"
#include <cmath>

// ---- Roe flux ---------------------------------------------------------------
void roe_flux(double rL, double uL, double vL, double pL,
              double rR, double uR, double vR, double pR,
              double nx, double ny, double gamma,
              double* F) {
    double gm1 = gamma - 1.0;

    // Total energies (per unit mass) and enthalpies
    double EL = pL / (gm1 * rL) + 0.5*(uL*uL + vL*vL);
    double ER = pR / (gm1 * rR) + 0.5*(uR*uR + vR*vR);
    double HL = EL + pL / rL;
    double HR = ER + pR / rR;

    // Physical normal fluxes
    double qnL = uL*nx + vL*ny;
    double qnR = uR*nx + vR*ny;
    double FL[4] = { rL*qnL, rL*uL*qnL + pL*nx, rL*vL*qnL + pL*ny, rL*HL*qnL };
    double FR[4] = { rR*qnR, rR*uR*qnR + pR*nx, rR*vR*qnR + pR*ny, rR*HR*qnR };

    // Roe-averaged state
    double sL    = std::sqrt(rL);
    double sR    = std::sqrt(rR);
    double denom = sL + sR;
    double u_h   = (sL*uL + sR*uR) / denom;
    double v_h   = (sL*vL + sR*vR) / denom;
    double H_h   = (sL*HL + sR*HR) / denom;
    double qn_h  = u_h*nx + v_h*ny;
    double a2    = gm1 * (H_h - 0.5*(u_h*u_h + v_h*v_h));
    if (a2 < 1e-10) a2 = 1e-10;
    double a_h   = std::sqrt(a2);

    // Tangential unit vector: t = (-ny, nx)
    double tx = -ny, ty = nx;
    double qt_h = u_h*tx + v_h*ty;

    // Eigenvalues
    double lam1 = qn_h - a_h;
    double lam4 = qn_h + a_h;
    double lam2 = qn_h;   // contact + shear (same eigenvalue)

    // Harten entropy fix on acoustic eigenvalues only
    double delta = 0.1 * a_h;
    auto efix = [delta](double lam) -> double {
        double al = std::abs(lam);
        return al < delta ? (lam*lam + delta*delta) / (2.0 * delta) : al;
    };
    double al1 = efix(lam1);
    double al4 = efix(lam4);
    double al2 = std::abs(lam2);

    // Jump in pressure, normal velocity, tangential velocity
    double dp  = pR - pL;
    double dqn = qnR - qnL;
    double qt_L = uL*tx + vL*ty;
    double qt_R = uR*tx + vR*ty;
    double dqt = qt_R - qt_L;

    // Roe-mean density for wave-strength computation
    double r_h = sL * sR;  // = sqrt(rL*rR)

    // Wave strengths (Steger-Warming / Roe decomposition)
    double alpha1 = (dp - r_h*a_h*dqn) / (2.0*a2);
    double alpha4 = (dp + r_h*a_h*dqn) / (2.0*a2);
    double alpha2 = (rR - rL) - dp / a2;     // entropy wave amplitude
    double alpha3 = r_h * dqt;               // shear wave amplitude

    // Dissipation: |A_Roe|*(UR-UL) = Σ |λ_k| α_k r_k
    // Right eigenvectors:
    //   r1 = [1,  u_h-a_h*nx,  v_h-a_h*ny,  H_h-qn_h*a_h ]
    //   r2 = [1,  u_h,         v_h,          0.5(u_h²+v_h²)]  (entropy)
    //   r3 = [0,  tx,          ty,            qt_h          ]  (shear)
    //   r4 = [1,  u_h+a_h*nx,  v_h+a_h*ny,  H_h+qn_h*a_h ]
    double ek = 0.5*(u_h*u_h + v_h*v_h);
    double d[4];
    d[0] = al1*alpha1
         + al2*alpha2
       /*+ al2*0.0 (r3[0]=0)*/
         + al4*alpha4;

    d[1] = al1*alpha1*(u_h - a_h*nx)
         + al2*alpha2*u_h
         + al2*alpha3*tx
         + al4*alpha4*(u_h + a_h*nx);

    d[2] = al1*alpha1*(v_h - a_h*ny)
         + al2*alpha2*v_h
         + al2*alpha3*ty
         + al4*alpha4*(v_h + a_h*ny);

    d[3] = al1*alpha1*(H_h - qn_h*a_h)
         + al2*alpha2*ek
         + al2*alpha3*qt_h
         + al4*alpha4*(H_h + qn_h*a_h);

    for (int k = 0; k < 4; ++k)
        F[k] = 0.5*(FL[k] + FR[k]) - 0.5*d[k];
}
