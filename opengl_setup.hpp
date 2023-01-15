#ifndef OPENGL_SETUP_HPP_
#define OPENGL_SETUP_HPP_

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <string>

extern float* image;
extern int width;
extern int height; 

void runSimulation(char*);
void keyboard(unsigned char, int, int);
void mouse(int, int, int, int);
void render(void);
void loop(void);
void loop_code(void);

#endif // OPENGL_SETUP_HPP_
