/**
 * @file fluidsolver.h
 * @author Caleb Johnston
 * @version 0.1
 *
 * @section LICENSE
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef __FDL_FLUIDSOLVER_H
#define __FDL_FLUIDSOLVER_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>

#include "core/common.h"
#include "core/vector.hpp"
#include "core/grid.hpp"

namespace fdl {

class FluidSolver {
public:
	FluidSolver(fdl::Grid* _grid);
	~FluidSolver();
	
	void step(float dt=0);
	void project(float dt);
	void applyForces(float dt);
	void advect(float dt);
	void diffuse(float dt, float rate=0.25f);
	void diffuseHeat(float dt, float rate);
	void applyViscosity(float dt);
	void setCGTolerance(float tol); 
	void setCGMaxIter(unsigned N);
	
	const Vector& getCurlMagnitude() const { return m_curlMagnitude; }
	const Vector& getDivergence() const { return m_divergence; }
	const Vector& getPressure() const { return m_pressure; }
	
protected:
	float computeMaxTimeStep() const;
	void axpy_prod(const Vector& x, Vector& y) const;
	void solvePreconditioner(const Vector& b, Vector& x);
	void constructMatrix(float dx, float dt, float rho=0.25f);//, bool variable_density=false);
	//void constructMatrix(float dx, float dt, float rho=0.25f, bool variable_density=false);
	void constructPreconditioner(float rho=0.25f, float tau=0.97f);
	float cgSolve(const Vector& b, Vector& x);
	float pcgSolve(const Vector& b, Vector& x);
	
private:
	/* Vorticity confinement vectors */
	std::vector<fdl::Vector3> m_curl;
	std::vector<fdl::Vector3> m_vorticityConfinementForce;
	Vector m_curlMagnitude;
	
	/* Grid discretized domain */
	fdl::Grid* grid;

	/* Vectors for the Hodge decomposition */
	Vector m_divergence;
	Vector m_pressure;

	/* Temporary vectors for the PCG iteration */
	Vector m_tempP, m_tempW, m_tempZ, m_tempR, m_tempQ;

	/* Pressure matrix */
	Vector m_ADiag, m_APlusX, m_APlusY, m_APlusZ;
	
	/* Modified incomplete cholesky preconditioner */
	Vector m_precond;
	
	/* Linear algebra stopping conditions */
	float tol_cg; 
	unsigned maxiter_cg;
	
	/* Gravity vector */
	Vector3 m_gravity;

	/* Grid resolution */
	int m_gridX;
	int m_gridY;
	int m_gridZ;
	int m_numPoints;
	int m_slice;
	int m_velSlice;
	float m_time;
	float m_dx;
	
	// performs 1D convolution in place
	template<class T>
	static void convolve(const std::vector<double>& kernel, const std::vector<T>& signal, std::vector<T>& output)
	{
		// allocate the output
		std::vector<double>* temp = new std::vector<double>();
		temp->resize(signal.size());
		temp->clear();
	
		// check for degenerate input
		if(kernel.empty()){
			ERROR() << "convolve - input kernel is empty!";
			return;
		}

		// check for degenerate input
		if(signal.empty()){
			ERROR() << "convolve - input signal is empty!";
			return;
		}
	
		int kernel_offset = (int) kernel.size()/2;
	
		// convolve
		int i, j;
		for(i=0; i<signal.size(); ++i){
			for(j=0; j<kernel.size(); ++j){
				int index = i + (j - kernel_offset);
				if(index > 0 && index < signal.size()){
					temp->at(i) += signal.at(index) * kernel.at(j);
				}
			}
		}

		// Copy results to the output just in case the output was the same as the signal...
		if(output.size() != signal.size()){
			output.resize(signal.size());
		}
		
		output = *temp;
	}
	
	// this function really doens't work
	static std::vector<double>* gaussianHeatKernel(const int size, double diffusion, float dt)
	{
		// allocate
		std::vector<double>* gauss = new std::vector<double>();
		gauss->resize(size);

		// OLD
		/*
		double N = (double)size - 1;
		double alpha = max(2.0, std);
		double exponent, sample, t;
		
		// compute
		for (int i=0; i<size; i++) {
			double n = (double)i - N/2;
			double t = alpha * (n/(N/2));
			exponent = -0.5*(t*t);
			sample = exp(exponent);
			gauss->at(i) = sample;
			//cout << sample << endl;
		}
		*/
		
		double x;
		int i = 0;
		double j = size - 1;
		std::vector<double>::iterator it;
		for(it = gauss->begin(); it < gauss->end(); it++){
			double denom = std::pow(4.0 * PI * diffusion * dt, 3.0/2.0);
			
			x = i - j/2.0;
			*it = (1.0/denom) * std::exp(-(x*x)/(4.0 * diffusion * dt));

			// test
			// std::cout << "gauss[" << i << "] = " << *it << std::endl;
			i++;
		}

		return gauss;
	}
	
};

}	// namespace fdl

#endif // __FDL_FLUIDSOLVER_H
