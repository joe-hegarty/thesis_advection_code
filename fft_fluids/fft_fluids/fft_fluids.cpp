/*
A heavily modified version of Jos Stams GDC2003 "Real-Time Fluid Dynamics
for Games" code.  Original commment from the code below (note that the
file name has changed).
*/

/*
  ======================================================================
   demo.c --- protoype to show off the simple solver
  ----------------------------------------------------------------------
   Author : Jos Stam (jstam@aw.sgi.com)
   Creation Date : Jan 9 2003

   Description:

	This code is a simple prototype that demonstrates how to use the
	code provided in my GDC2003 paper entitles "Real-Time Fluid Dynamics
	for Games". This code uses OpenGL and GLUT for graphics and interface

  =======================================================================
*/
#include "stdafx.h"

#include "fft_fluids.h"
#include "FftSolverExperiment.h"
#include "RotationFieldExperiment.h"
#include "ConstVectorFieldExperiment.h"

#include <stdlib.h>
#include <GL/glut.h>

#include <vector>

enum ExperimentType {fftSim, rotate, constVec};

ExperimentType experimentType = fftSim;

static int win_id;
static int dvel;

std::unique_ptr<Experiment> experiment;

static void pre_display ( void )
{
	glViewport ( 0, 0, experiment->win_x(), experiment->win_y() );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ();
	gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
	glClearColor ( 0.0, 0.0, 0.0, 1.0 );
	glClear ( GL_COLOR_BUFFER_BIT );
}

static void post_display ( void )
{
	glutSwapBuffers ();
}

static void draw_velocity(void)
{
  experiment->draw_velocity();
}

static void draw_density(void)
{
  experiment->draw_density();
}

static void key_func ( unsigned char key, int x, int y )
{
	switch ( key )
	{
		case 'c':
		case 'C':
			experiment->clear_data();
			break;

		case 'q':
		case 'Q':
			experiment->free_data();
			exit ( 0 );
			break;

		case 'v':
		case 'V':
			dvel = ! dvel;
			break;
	}
}

static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	experiment->win_x(width);
	experiment->win_y(height);
}


static void idle_func(void)
{
  experiment->idle_func();

  glutSetWindow(win_id);
  glutPostRedisplay();
}

static void display_func(void)
{
  pre_display();

  if ( dvel )
    draw_velocity();
  else
    draw_density();

  post_display();
}

static void open_glut_window ( void )
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(experiment->win_x(), experiment->win_y());
	win_id = glutCreateWindow("2D Fluids");

	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	pre_display();

	glutKeyboardFunc(key_func);
	glutReshapeFunc(reshape_func);
  glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
}

int main ( int argc, char ** argv )
{
	glutInit ( &argc, argv );

  switch (experimentType) {
  case fftSim:
      experiment.reset(new FftSolverExperiment);
      break;
  case rotate:
      experiment.reset(new RotationFieldExperiment(2.0 * acos(-1)));
      break;
  case constVec:
    double u = 1.0;
    double v = (1.0 / 8.0) * u;
    experiment.reset(new ConstVectorFieldExperiment(Vector2<double>(u, v)));
    break;
  }

	experiment->win_x(512);
	experiment->win_y(512);

  dvel = 0;

	if (! experiment->allocate_data()) {
    exit(1);
  }

	experiment->initialise_data();

	open_glut_window();

	glutMainLoop();

	exit(0);
}