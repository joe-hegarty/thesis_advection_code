#pragma once
#ifndef _EXPERIMENT_H_
#define _EXPERIMENT_H_

#include "vector_field.h"

#include <GL/glut.h>

#include <memory>

enum AdvectMethod {SemiLagragian, MacCormack, BFECC};

AdvectMethod scalarAdvectMethod = MacCormack;

struct Experiment
{
  Experiment()
    : mouse_down_(std::vector<int>(3))
  {}

  virtual void advectScalars()
  {
    std::swap(d, d0);

    switch (scalarAdvectMethod) {
      case SemiLagragian:
        advect(N, &d[0], &d0[0], dt);
        break;
      case MacCormack:
        advect_MacCormack(N, &d[0], &d0[0], &temp[0], &posX[0], &posY[0], dt);
        break;
      case BFECC:
        advect_BFECC(N, &d[0], &d0[0], &temp[0], &posX[0], &posY[0], dt);
        break;
    }
  }

  virtual Vector2<double> getVelocity(size_t i, size_t j) = 0;

  virtual double getDensity(size_t i, size_t j)
  {
    return d[i + j*N];
  }

  virtual void idle_func()
  {
    if ( isActive() ) {
      addForces();

      project();

      advectScalars();

      advectVelocity();

      emit();
    }
  }

  virtual bool isActive()
  {
    return simulation_length-- > 0;
  }

  virtual void advectVelocity()
  {}

  virtual void addForces()
  {}

  virtual void project()
  {}

  virtual void emit()
  {}

  virtual void initialise_data()
  {
    clear_data();

    int size = N * N;
    int period = 8;
    int spacing = (N / period < 1) ? 1 : N / period;

    double value = 0;

    for (int j = 0; j < N; j++) {
      if (j % spacing == 0) {
        value = 1.0 - value;
      }

      for (int i = 0; i < N; i++) {
        if (i % spacing == 0) {
          value = 1.0 - value;
        }

        d[i + j*N] = value;
      }
    }
  }

  virtual bool allocate_data()
  {
    int size = N * N;

    d.resize(size);
    d0.resize(size);

    temp.resize(size);
    temp2.resize(size);
    posX.resize(size);
    posY.resize(size);

    return true;
  }

  virtual void free_data()
  {}

  virtual void clear_data()
  {
    int size = N * N;

    for (int i = 0; i < size; i++) {
      d[i] = d0[i] = temp[i] = temp2[i] = posX[i] = posY[i] = 0.0;
    }
  }

  virtual void draw_density()
  {
    int i, j;
    double x, y, h;

    h = 1.0 / N;

    glBegin ( GL_QUADS );

    for (j = 0; j < N; j++) {
      y = j * h;
      for (i = 0; i < N; i++) {
        x = i * h;

        double v = getDensity(i, j);

        glColor3f ( v, v, v ); glVertex2f ( x,   y );
        glColor3f ( v, v, v ); glVertex2f ( x+h, y );
        glColor3f ( v, v, v ); glVertex2f ( x+h, y+h );
        glColor3f ( v, v, v ); glVertex2f ( x,   y+h );
      }
    }

    glEnd ();
  }

  virtual void draw_velocity()
  {
    int i, j;
    double x, y, h;

    h = 1.0f / N;

    glColor3f(1.0, 1.0, 0.0);
    glLineWidth(1.0);

    glBegin(GL_LINES);

    int step = (N < 256) ? 1 : (N / 256);
    int offset = step / 2;
    for (j = offset; j < N; j += step) {
      y = (j + 0.5) * h;
      for (i = offset; i < N; i += step) {
        x = (i + 0.5) * h;

        const Vector2<double>& vel = getVelocity(i, j);

        glVertex2f(x, y);
        glVertex2f(x + vel.x(), y + vel.y());
      }
    }

    glEnd();
  }

  int win_x() { return win_x_; }
  int win_y() { return win_y_; }
  void win_x(int value) { win_x_ = value; }
  void win_y(int value) { win_y_ = value; }

  void advect(int n, double* c, const double* c0, double dt)
  {
    double x, y, x0, y0, s, t;
    int i0, j0, i1, j1;

    for ( int j = 0; j < n; j++ ) {
      y = (0.5 + j) / n;
      for ( int i = 0; i < n; i++ ) {
        x = (0.5 + i) / n;

        const Vector2<double> vel = getVelocity(i, j);

        x0 = n * (x - dt*vel.x()) - 0.5;
        y0 = n * (y - dt*vel.y()) - 0.5;

        i0 = floor(x0);
        s = x0 - i0;
        i0 = (n + (i0%n)) % n;
        i1 = (i0 + 1) % n;

        j0 = floor(y0);
        t = y0 - j0;
        j0 = (n + (j0%n)) % n;
        j1 = (j0 + 1) % n;

        c[i + n*j] = (1 - s)*((1 - t)*c0[i0 + n*j0] + t*c0[i0 + n*j1]) +
                          s *((1 - t)*c0[i1 + n*j0] + t*c0[i1 + n*j1]);
      }
    }
  }

  void advectWithPos(int n, double* c, const double* c0, double* xPos, double* yPos, double dt)
  {
    double x, y, x0, y0, s, t;
    int i0, j0, i1, j1;

    for ( int j = 0; j < n; j++ ) {
      y = (0.5 + j) / n;
      for ( int i = 0; i < n; i++ ) {
        x = (0.5 + i) / n;

        const Vector2<double> vel = getVelocity(i, j);

        x0 = n * (x - dt*vel.x()) - 0.5;
        y0 = n * (y - dt*vel.y()) - 0.5;

        i0 = floor(x0);
        s = x0 - i0;
        i0 = (n + (i0%n)) % n;
        i1 = (i0 + 1) % n;

        j0 = floor(y0);
        t = y0 - j0;
        j0 = (n + (j0%n)) % n;
        j1 = (j0 + 1) % n;

        c[i + n*j] = (1 - s)*((1 - t)*c0[i0 + n*j0] + t*c0[i0 + n*j1]) +
                          s *((1 - t)*c0[i1 + n*j0] + t*c0[i1 + n*j1]);

        xPos[i + n*j] = x0;
        yPos[i + n*j] = y0;
      }
    }
  }

  void clamp_limiter(int n, double* c, const double* c0, double* xPos, double* yPos, double dt)
  {
    double x0, y0, s, t;
    int i0, j0, i1, j1;

    for ( int j = 0; j < n; j++ ) {
      for ( int i = 0; i < n; i++ ) {
        x0 = xPos[i + n*j];
        y0 = yPos[i + n*j];

        i0 = floor(x0);
        s = x0 - i0;
        i0 = (n + (i0%n)) % n;
        i1 = (i0 + 1) % n;

        j0 = floor(y0);
        t = y0 - j0;
        j0 = (n + (j0%n)) % n;
        j1 = (j0 + 1) % n;

        double value = c[i + n*j];

        double min = std::min(c0[i0 + n*j0], std::min(c0[i0 + n*j1], std::min(c0[i1 + n*j0], c0[i1 + n*j1])));
        double max = std::max(c0[i0 + n*j0], std::max(c0[i0 + n*j1], std::max(c0[i1 + n*j0], c0[i1 + n*j1])));

        value = std::min(value, max);
        value = std::max(value, min);

        c[i + n*j] = value;
      }
    }
  }

  void advect_BFECC(int n, double* c, const double* c0, double* temp, double* xPos, double* yPos, double dt)
  {
    advectWithPos(n, c, c0, xPos, yPos, dt);
    advect(n, temp, c, -dt);

    int size = n * n;
    for ( int i = 0; i < size; i++ )
      temp[i] = c0[i] + 0.5 * (c0[i] - temp[i]);

    advect(n, c, temp, dt);

    if (limited)
      clamp_limiter(n, c, c0, xPos, yPos, dt);
  }

  void advect_MacCormack(int n, double* c, const double* c0, double* temp, double* xPos, double* yPos, double dt)
  {
    advectWithPos(n, c, c0, xPos, yPos, dt);
    advect(n, temp, c, -dt);

    int size = n * n;
    for ( int i = 0; i < size; i++ )
      c[i] += 0.5 * (c0[i] - temp[i]);

    if (limited)
      clamp_limiter(n, c, c0, xPos, yPos, dt);
  }

  std::vector<int> mouse_down_;

  int omx_ = 0;
  int omy_ = 0;
  int mx_ = 0;
  int my_ = 0;
  int win_x_ = 0;
  int win_y_ = 0;

  std::vector<double> d;
  std::vector<double> d0;

  // Needed for BFECC / MacCormack
  std::vector<double> temp;
  std::vector<double> temp2;
  std::vector<double> posX;
  std::vector<double> posY;

  bool limited = true;

	int N = 128;
  
  double 	dt = 1.0 / 20.0;
  int simulation_length = 400 / dt;
};

#endif
