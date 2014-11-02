#pragma once
#ifndef _VECTOR_2_H_
#define _VECTOR_2_H_

template <class T>
class Vector2
{
public:
  // Member not initialised in default constructor
  Vector2() {}

  Vector2(T x, T y) : x_(x), y_(y) {}

  T x() const { return x_; }
  T y() const { return y_; }

  void x(T x) { x_ = x; }
  void y(T y) { y_ = y; }

  T& operator[](size_t i)
  {
    if (i == 0)
      return x_;
    else if (i == 1)
      return y_;
    else
      throw std::exception("Invalid 2d index");
  }

  const T& operator[](size_t i) const
  {
    if (i == 0)
      return x_;
    else if (i == 1)
      return y_;
    else
      throw std::exception("Invalid 2d index");
  }

private:
  T x_;
  T y_;
};

#endif
