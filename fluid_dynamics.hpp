#ifndef FLUID_DYNAMICS_HPP_
#define FLUID_DYNAMICS_HPP_

class Fluid {
	private:
		int width_; 
		int height_;	
		float dt_;
		float viscosity_;
		float diffusion_rate_;
		float* vx_;
		float* vy_;
		float* vx0_;
		float* vy0_;
		float* dens_;
		float* s_;
		void setBoundary(int, float*);
		void solveLinear(int, float*, float*, float, float, int);
		void diffuse(int, float*, float*, int);
		void project(float*, float*, float*, float*, int);
		void advect(int, float*, float*, float*, float*);
	public:
		Fluid(float*, int, int, float, float, float);
		~Fluid();
		void addSource(int, int, float);
		void addVelocity(int, int, float, float);
		void evaluate(int);
};

#endif // FLUID_DYNAMICS_HPP_
