#include "opengl_setup.hpp"
#include "fluid_dynamics.hpp"
#include <iostream>
#include <omp.h>

using namespace std;

int width = 300;
int height = 300;
float* image = new float[width * height];
Fluid f(image, width, height, 0.1, 0.0, 0.01, 10);			
int brush_size = 100;

void loop_code() {
	// post fluid code here
	float speed = 100;
 	for(int k = -5; k < 5; ++k) {
//	f.addVelocity(width/2+k, height/2+5, 0, speed); 
//	f.addVelocity(width/2+k, height/2-5, 0, -speed); 
//	f.addVelocity(width/2-5, height/2+k, -speed, 0); 
//	f.addVelocity(width/2+5, height/2+k, speed, 0); 
  	f.addVelocity(width/3+k, height/3, 0, speed); 
  	f.addVelocity(2*width/3, height/3+k, -speed, 0); 
  	f.addVelocity(width/3, 2*height/3+k, speed, 0); 
  	f.addVelocity(2*width/3+k, 2*height/3, 0, -speed); 
  }
	f.evaluate();
	return;
}

int main() {
	runSimulation("Fluid Simulator");
	return 0;
}

void mouse(int button, int state, int x, int y){
	if(button == GLUT_RIGHT_BUTTON)
		std::fill_n(image, width*height, 0);
	if(button == GLUT_LEFT_BUTTON) {
		for(int i = -brush_size/2; i <= brush_size/2; ++i)
			for(int j = -brush_size/2; j <= brush_size/2; ++j)
				f.addSource(x+i, height-y+j, 1.0);
	}
}

void keyboard(unsigned char c, int x, int y) {
  if(c == 27) {
	float sum = 0;
	float d;
	/*
	for(int i = 0; i < width; ++i)
		for(int j = 0; j < height; ++j) {
  		d = image[i+j*width];
			sum += d;	
		}
	cout << endl << "sum of all densities field-wide: " << sum/(width*height);
	*/
  exit(0);
	
}
}


