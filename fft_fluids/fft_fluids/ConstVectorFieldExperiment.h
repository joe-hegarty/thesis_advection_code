#ifndef _CONST_VECTOR_FIELD_EXPERIMENT_
#define _CONST_VECTOR_FIELD_EXPERIMENT_

#include "Experiment.h"

struct ConstVectorFieldExperiment : public Experiment
{
  ConstVectorFieldExperiment(Vector2<double> value)
    : vectorField(value)
  {
    limited = true;

  	N = 128;

    dt = 1.0 / 1000.0f;

    simulation_length = 1.0 / dt;
  }

  virtual Vector2<double> getVelocity(size_t i, size_t j)
  {
    return vectorField;
  }

private:
  Vector2<double> vectorField;
};

#endif
