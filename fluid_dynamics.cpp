#include "fluid_dynamics.hpp"
#include "opengl_setup.hpp"
#include <math.h>
#include <iostream>

using namespace std;

#define I(x,y) ((x)+(width_)*(y))
#define SWAP_ARRAYS(a,b) {float* tmp=a;a=b;b=tmp;}

Fluid::Fluid(float* image = NULL,
						 int width = 100, 
						 int height = 100, 
						 float dt = 0.1,
						 float viscosity = 0.1, 
						 float diffusion_rate = 0.1) 
						 : dens_(image),
							 width_(width), 	
							 height_(height),
							 dt_(dt),
						 	 viscosity_(viscosity), 
						 	 diffusion_rate_(diffusion_rate) {
	vx_ = new float[width * height];
	vy_ = new float[width * height];
	vx0_ = new float[width * height];
	vy0_ = new float[width * height];
	s_ = new float[width * height];
};

void Fluid::addSource(int x, int y, float amount) {
	dens_[I(x,y)] += amount;
}

void Fluid::addVelocity(int x, int y, float x_amount, float y_amount) {
	vx_[I(x,y)] += x_amount;
	vy_[I(x,y)] += y_amount;
}

void Fluid::evaluate(int iterations) {
	SWAP_ARRAYS(vx0_, vx_); diffuse(1, vx_, vx0_, iterations);
	SWAP_ARRAYS(vy0_, vy_); diffuse(2, vy_, vy0_, iterations);
	project(vx_, vy_, vx0_, vy0_, iterations);

	SWAP_ARRAYS(vx0_, vx_); SWAP_ARRAYS(vy0_, vy_);
	advect(1, vx_, vx0_, vx0_, vy0_);
	advect(2, vy_, vy0_, vx0_, vy0_);
	project(vx_, vy_, vx0_, vy0_, iterations);
	
	SWAP_ARRAYS(dens_, s_); diffuse(0, dens_, s_, iterations);
  SWAP_ARRAYS(dens_, s_); advect(0, dens_, s_, vx_, vy_);
}

void Fluid::diffuse(int b, float* x, float* x0, int iter) {
	float a = dt_ * diffusion_rate_ * (width_) * (height_);
	solveLinear(b, x, x0, a, 1+4*a, iter);
}

void Fluid::project(float* vx, float* vy, float* p, float* div, int iter) {
	for(int j = 1; j < height_ - 1; ++j)
		for(int i = 1; i < width_ - 1; ++i) {
			div[I(i, j)] = -0.5f * (vx[I(i+1, j)] - vx[I(i-1, j)] 
														+ vy[I(i, j+1)] - vy[I(i, j-1)]) / width_;
			p[I(i,j)] = 0;
		}
	setBoundary(0, div);
	setBoundary(0, p);	
	
	solveLinear(0, p, div, 1, 6, iter);

	for(int j = 1; j < height_ - 1; ++j)
		for(int i = 1; i < width_ - 1; ++i) {
			vx[I(i, j)] -= 0.5f * (p[I(i+1, j)] - p[I(i-1, j)]) * width_; 
			vy[I(i, j)] -= 0.5f * (p[I(i, j+1)] - p[I(i, j-1)]) * height_; 
		}
	setBoundary(1, vx);
	setBoundary(2, vy);
}

void Fluid::advect(int b, float* d, float* d0, float* vx, float* vy) {
	float x, y;
	for(int j = 1; j < height_ - 1; ++j)
		for(int i = 1; i < width_ - 1; ++i) {
			x = ((float) i) - dt_ * (width_-2) * vx[I(i,j)];
			y = ((float) j) - dt_ * (height_-2) * vy[I(i,j)];
			if(x < 0.5f) x = 0.5f;
			if(x > (((float) width_) + 0.5f)) x = ((float) width_) + 0.5f;
			float i0 = floorf(x);
			float i1 = i0 + 1.0f;
			if(y < 0.5f) y = 0.5f;
			if(y > (((float) height_) + 0.5f)) y = ((float) height_) + 0.5f;
			float j0 = floorf(y);
			float j1 = j0 + 1.0f;
			
			float s1 = x - i0;
			float s0 = 1.0f - s1;
			float t1 = y - j0;
			float t0 = 1.0f - t1;

			int i0i = (int) i0;
			int i1i = (int) i1;
			int j0i = (int) j0;
			int j1i = (int) j1;
			d[I(i, j)] = s0 * (t0 * d0[I(i0i, j0i)] + t1 * d0[I(i0i, j1i)])
							 	 + s1 * (t0 * d0[I(i1i, j0i)] + t1 * d0[I(i1i, j1i)]);
		}
	setBoundary(b, d);
}

void Fluid::setBoundary(int b, float* x) {
		// mirror values at top and bottom borders
		for(int i = 1; i < width_-1; i++) {
			x[I(i, 0			 )] = b == 2 ? -x[I(i,        1)] : x[I(i,        1)];
			x[I(i, height_-1)] = b == 2 ? -x[I(i, height_-2)] : x[I(i, height_-2)];
		}
		// mirror values at left and right borders
		for(int j = 1; j < height_-1; j++) {
			x[I(0, j			  )] = b == 1 ? -x[I(1,         j)] : x[I(1        , j)];
			x[I(width_ - 1, j)] = b == 1 ? -x[I(width_ - 2, j)] : x[I(width_ - 2, j)];
		}
		// set values at the vortexes
		x[I(0, 0)] = 0.33f * (x[I(1,0)] + x[I(0,1)]);
		x[I(width_-1, 0)] = 0.33f * (x[I(width_-2,0)] + x[I(width_-1,1)]);
		x[I(width_-1, height_-1)] = 0.33f * (x[I(width_-2,height_-1)] + x[I(width_-1,height_-2)]);
		x[I(0, height_-1)] = 0.33f * (x[I(0,height_-2)] + x[I(1,height_-1)]);
}

void Fluid::solveLinear(int b, float* x, float* x0, float a, float c, int iter) {
	float cRecip = 1.0 / c;
	for(int k = 0; k < iter; k++) {
		for(int j = 1; j < height_-1; j++) {
			for(int i = 1; i < width_-1; i++) {
				x[I(i, j)] = (x0[I(i, j)] + a * (x[I(i+1, j)] + x[I(i-1, j)]
																			 + x[I(i, j+1)] + x[I(i, j-1)])) * cRecip;
			}
		}
		setBoundary(b, x);
	}
} 

Fluid::~Fluid() {
//delete[] vx_; 
//delete[] vy_;
//delete[] vx0_;
//delete[] vx0_;
//delete[] s_;
}
