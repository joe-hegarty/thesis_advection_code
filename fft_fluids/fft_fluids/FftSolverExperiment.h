#ifndef _FFT_SOLVER_EXPERIMENT_
#define _FFT_SOLVER_EXPERIMENT_

#include "Experiment.h"

struct FftSolverExperiment : public Experiment
{
  FftSolverExperiment()
  {}

  virtual double getDensity(size_t i, size_t j)
  {
    return d[i + j*N];
  }

  virtual Vector2<double> getVelocity(size_t i, size_t j)
  {
    return Vector2<double>(u[i + j*N], v[i + j*N]);
  }

  virtual bool allocate_data()
  {
    Experiment::allocate_data();

    int size = N * N;

    t.resize(size);
    t0.resize(size);

    u.resize(size);
    v.resize(size);

    u0 = allocate_inplace_memory_2d(N);
    v0 = allocate_inplace_memory_2d(N);

    if ( !u0 || !v0 ) {
      fprintf ( stderr, "cannot allocate data\n" );
      return false;
    }

    init_FFT(N, u0, plan_r2c_u0, plan_c2r_u0);
    init_FFT(N, v0, plan_r2c_v0, plan_c2r_v0);

    return true;
  }

  virtual void free_data()
  {
    Experiment::free_data();

    fftw_destroy_plan(plan_r2c_u0);
    fftw_destroy_plan(plan_c2r_u0);
    fftw_destroy_plan(plan_r2c_v0);
    fftw_destroy_plan(plan_c2r_v0);

    free_2d(u0);
    free_2d(v0);
  }

  void clear_data()
  {
    Experiment::clear_data();

    int size = N * N;

    for (int i = 0 ; i < size ; i++ ) {
      t[i] = t0[i] = 0.0;
    }

    size = N * 2 * (N / 2 + 1);

    for (int i = 0 ; i < size ; i++ ) {
      u0[i] = v0[i] = 0.0;
    }
  }

  virtual void addForces()
  {
    // Note that even though we are actually dealing with u0 and v0 here we still operate
    // over the first N * N elements, we do not need to worry about the buffer zone
    int size = N * N;
    for (int i = 0; i < size; i++ ) {
      t0[i] = 0.0;
      u0[i] = 0.0;
      v0[i] = 0.0;
    }

    buoyancy(t, buoyancy_multiplyer);

    // Initially u0 and v0 stores the velocity forces.  We add these to u and v scaled by dt
    for (int i = 0; i<N*N; i++) {
      t[i] += dt * t0[i];
      u[i] += dt * u0[i];
      v[i] += dt * v0[i];
    }
  }

  virtual void project()
  {
    stable_solve(N, u, v, temp, u0, plan_r2c_u0,
                 plan_c2r_u0, v0, plan_r2c_v0,
                 plan_c2r_v0, visc, dt);
    
    // Copy new velocity into u and v so it's available in getVelocity function immediately
    int size = N * N;
    for (int i = 0; i < size; i++) {
      u[i] = u0[i];
      v[i] = v0[i];
    }
  }

  virtual void advectScalars()
  {
    Experiment::advectScalars();

    std::swap(t, t0);

    switch (scalarAdvectMethod) {
      case SemiLagragian:
        advect(N, &t[0], &t0[0], dt);
        break;
      case MacCormack:
        advect_MacCormack(N, &t[0], &t0[0], &temp[0], &posX[0], &posY[0], dt);
        break;
      case BFECC:
        advect_BFECC(N, &t[0], &t0[0], &temp[0], &posX[0], &posY[0], dt);
        break;
    }
  }

  virtual void advectVelocity()
  {
    advect(N, &temp[0], u0, dt);
    advect(N, &temp2[0], v0, dt);

    std::swap(u, temp);
    std::swap(v, temp2);
  }

  virtual void emit()
  {
    if (emission_length-- > 0)
      emit(t, dt * source, 0.0675, 0.0675, 0.5, 0.125);
  }

private:
  void buoyancy(std::vector<double>& s, double factor)
  {
    double ambient = 0;
    int size = N * N;
    for (int i = 0; i < size; i++) {
      ambient += s[i];
    }
    ambient /= size;

    // for each cell compute buoyancy force
    for ( int j = 0; j < N; j++ ) {
      for ( int i = 0; i < N; i++ ) {
        v0[i + j*N] += dt * (s[i + j*N] - ambient) * factor;
      }
    }
  }

  void emit(std::vector<double>& s, double value, double halfWidth, double halfHeight, double x, double y)
  {
    x *= N;
    y *= N;
    halfWidth *= N;
    halfHeight *= N;

    int lowerLeft_i = x - halfWidth;
    int lowerLeft_j = y - halfHeight;
    int upperRight_i = x + halfWidth;
    int upperRight_j = y + halfHeight;

    if ( lowerLeft_i < 0) lowerLeft_i = 0;
    if ( lowerLeft_j < 0) lowerLeft_j = 0;
    if (upperRight_i < 0) upperRight_i = 0;
    if (upperRight_j < 0) upperRight_j = 0;
    if (lowerLeft_i > N - 1) lowerLeft_i = N - 1;
    if (lowerLeft_j > N - 1) lowerLeft_j = N - 1;
    if (upperRight_i > N - 1) upperRight_i = N - 1;
    if (upperRight_j > N - 1) upperRight_j = N - 1;

    for (int j = lowerLeft_j; j <= upperRight_j; j++) {
      for (int i = lowerLeft_i; i <= upperRight_i; i++) {
        s[i + j*N] += value;
      }
    }
  }

  double source = 1.0;
  double visc = 0.0;

  double buoyancy_multiplyer = 0.01;

  int emission_length = 2 / dt;

  fftw_plan plan_r2c_u0, plan_c2r_u0;
  fftw_plan plan_r2c_v0, plan_c2r_v0;

  std::vector<double> u;
  std::vector<double> v;

  double* u0;
  double* v0;

  std::vector<double> t;
  std::vector<double> t0;
};

#endif