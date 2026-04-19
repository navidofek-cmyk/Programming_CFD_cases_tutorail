#pragma once

// Staggered (MAC) grid layout:
//
//   p(i,j)   : cell centre at ((i+0.5)*dx, (j+0.5)*dy),  size  Nx  x  Ny
//   u(i,j)   : face centre at ( i      *dx, (j+0.5)*dy),  size (Nx+1) x  Ny
//   v(i,j)   : face centre at ((i+0.5)*dx,  j      *dy),  size  Nx  x (Ny+1)
//
// Boundary indices:
//   u(0,j)   = inlet,  u(Nx,j)  = outlet
//   v(i,0)   = bottom wall,  v(i,Ny) = top wall

struct Mesh {
    int    Nx, Ny;
    double Lx, Ly, dx, dy;

    Mesh(int Nx, int Ny, double Lx, double Ly)
        : Nx(Nx), Ny(Ny), Lx(Lx), Ly(Ly), dx(Lx / Nx), dy(Ly / Ny) {}

    int u_nx() const { return Nx + 1; }
    int u_ny() const { return Ny; }
    int v_nx() const { return Nx; }
    int v_ny() const { return Ny + 1; }
    int p_nx() const { return Nx; }
    int p_ny() const { return Ny; }
};
