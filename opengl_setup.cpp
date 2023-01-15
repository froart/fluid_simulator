#include "opengl_setup.hpp"

void runSimulation(char* name) {
	int argc = 1;
  char* argv[1] = {(char*) "OpenGL"};
  glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(width, height);
	glutCreateWindow(name);
	glutDisplayFunc(render);
	glutIdleFunc(loop);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMainLoop();	
  return;
}

void render(void) {
  glClearColor(0, 0, 0, 1);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-1, 1, 1, -1, 1, 0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glEnable(GL_TEXTURE_2D); 
  glBindTexture(GL_TEXTURE_2D, 1);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glClear(GL_COLOR_BUFFER_BIT);
  glTexImage2D(GL_TEXTURE_2D, 0, 1, width, height, 0, GL_RED, GL_FLOAT, image);
  glBegin(GL_QUADS);
    glTexCoord2f(0.f,0.f);  glVertex3f(-1, 1, -1);
    glTexCoord2f(1.f,0.f);  glVertex3f(1, 1, -1);
    glTexCoord2f(1.f,1.f);  glVertex3f(1, -1, -1);
    glTexCoord2f(0.f,1.f);  glVertex3f(-1, -1, -1);
  glEnd();
  glFlush();
  glutSwapBuffers();
  return;
}

void loop(void) {
	loop_code();
	glutPostRedisplay();
}
