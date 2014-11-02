#pragma once
#ifndef _VECTOR_FIELD_H_
#define _VECTOR_FIELD_H_

#include "vector2.h"

struct VectorField
{
  virtual Vector2<double> operator()(double x, double y, double t) = 0;
};

struct ConstVectorField : public VectorField
{
  ConstVectorField(Vector2<double> value)
    : value(value)
  {}
  
  virtual Vector2<double> operator()(double x, double y, double t)
  {
    return value;
  }

private:
  Vector2<double> value;
};

#endif
