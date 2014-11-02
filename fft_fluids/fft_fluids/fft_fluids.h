/*
An updated version of Jos Stams "A Simple Fluid Solver based on the FFT"
code built against FFTW version 3.
*/

#pragma once
#ifndef _FFT_FLUIDS_H_
#define _FFT_FLUIDS_H_

#include <fftw3.h>

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include <cmath>

static const std::string stable_solve_str("stable_solve");
static const std::string fftw_solve_str("fftw_solve");
static const std::string density_advect_str("density_advect");
static const std::string velocity_advect_str("velocity_advect");

double* allocate_inplace_memory_2d(int n)
{
  return fftw_alloc_real(n * 2 * (n / 2 + 1));
}

void free_2d(double* x)
{
  fftw_free(x);
}

void init_FFT(int n, double* x, fftw_plan& plan_r2c_x, fftw_plan& plan_c2r_x)
{
  plan_r2c_x = fftw_plan_dft_r2c_2d(n, n, x, (fftw_complex*)x, FFTW_MEASURE);
  plan_c2r_x = fftw_plan_dft_c2r_2d(n, n, (fftw_complex*)x, x, FFTW_MEASURE);
}

#define floor(x) ((x) >= 0.0 ? ((int)(x)) : (-((int)(1-(x)))))

void stable_solve(int n, std::vector<double>& u,
                  std::vector<double>& v, std::vector<double>& temp,
                  double* u0, fftw_plan plan_r2c_u0, fftw_plan plan_c2r_u0, 
                  double* v0, fftw_plan plan_r2c_v0, fftw_plan plan_c2r_v0,
                  double visc, double dt)
{
  double x, y, f, r, U[2], V[2];
  int i, j;
  int ps = 2 * (n / 2 + 1);

  // Copying from u to u0 and from v to v0.  Padding is important
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      u0[i + ps*j] = u[i + n*j];
      v0[i + ps*j] = v[i + n*j];
    }
  }

  fftw_execute(plan_r2c_u0);
  fftw_execute(plan_r2c_v0);

  for (j = 0; j < n; j++) {
    y = (j <= n / 2) ? j : j - n;
    for (i = 0; i <= n; i += 2) {
      x = 0.5 * i;
      r = x*x + y*y;

      if (r == 0.0)
        continue;

      f = exp(-r*dt*visc);
      
      U[0] = u0[i + ps*j];
      V[0] = v0[i + ps*j];
      U[1] = u0[i + 1 + ps*j];
      V[1] = v0[i + 1 + ps*j];

      u0[i     + ps*j] = f*( (1-x*x/r)*U[0] -    x*y/r *V[0] );
      u0[i + 1 + ps*j] = f*( (1-x*x/r)*U[1] -    x*y/r *V[1] );
      v0[i     + ps*j] = f*(   -y*x/r *U[0] + (1-y*y/r)*V[0] );
      v0[i + 1 + ps*j] = f*(   -y*x/r *U[1] + (1-y*y/r)*V[1] );
    }
  }

  fftw_execute(plan_c2r_u0);
  fftw_execute(plan_c2r_v0);

  // Pack u0 and v0 ready for advection
  f = 1.0 / (n*n);
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      u0[i + n*j] = f * u0[i + ps*j];
      v0[i + n*j] = f * v0[i + ps*j];
    }
  }
}

#endif
