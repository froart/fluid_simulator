#include "opengl_setup.hpp"
#include "fluid_dynamics.hpp"
#include <iostream>

using namespace std;

int width = 200;
int height = 200;
float* image = new float[width * height];
Fluid f(image, width, height, 0.1, 1000.0, 0.00001);			

void loop_code() {
	// post fluid code here
	f.evaluate(5);	
	return;
}

int main() {
	windowInit("Fluid Simulator");
	return 0;
}

void mouse(int button, int state, int x, int y){
	if(button == GLUT_RIGHT_BUTTON)
		exit(0);
	if(button == GLUT_LEFT_BUTTON) {
		for(int i = -3; i <= 3; ++i)
			for(int j = -3; j <= 3; ++j)
				f.addSource(width/2+i, height/2+j,  10.5);
	 	f.addVelocity(width/2, height/2, 10*(x-width/2), 10*(height-y-height/2)); 
	}
}
