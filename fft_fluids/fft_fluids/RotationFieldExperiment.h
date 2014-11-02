#ifndef _ROTATION_FIELD_EXPERIMENT_
#define _ROTATION_FIELD_EXPERIMENT_

#include "Experiment.h"

struct RotationFieldExperiment : public Experiment
{
  RotationFieldExperiment(double magnitude)
    : magnitude(magnitude)
  {
    limited = true;

  	N = 1024;

    dt = 1.0 / 1000.0f;

    simulation_length = 1.0 / dt;
  }

  virtual void initialise_data()
  {
    Experiment::initialise_data();

    // Black out values outside the radius of the circle
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < N; i++) {
        double i_norm = (double(i) / N) - 0.5;
        double j_norm = (double(j) / N) - 0.5;
        double dist = sqrt(i_norm*i_norm + j_norm*j_norm);

        if (dist > 0.40) {
          d[i + j*N] *= std::max(0.0, 20.0 * (0.45 - dist));
        }
      }
    }
  }

  virtual Vector2<double> getVelocity(size_t i, size_t j)
  {
    double i_norm = magnitude * ((double(i) / N) - 0.5);
    double j_norm = magnitude * ((double(j) / N) - 0.5);

    return Vector2<double>(-j_norm, i_norm);
  }

private:
  double magnitude;
};

#endif
