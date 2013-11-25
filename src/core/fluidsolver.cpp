#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <cfloat>

#include <png.h>

#include "core/common.h"
#include "core/fluidsolver.h"
#include "logger/logger.h"


namespace fdl {

/**
 * FluidSolver constructor collect information about the grid in order to allocate all
 * the data for the linear solve(s) and forces. 
 *
 * @param _grid domain to be solved
 *
 */
FluidSolver::FluidSolver(Grid* _grid)
{	
	// copy some grid params
	this->grid = _grid;
	m_numPoints = grid->getNumberOfGridCells();
	m_gridX = grid->getGridSizeX();
	m_gridY = grid->getGridSizeY();
	m_gridZ = grid->getGridSizeZ();
	m_slice = grid->getDensityGridSlice();
	m_velSlice = grid->getVelocityGridSlice();
	m_dx = grid->getVoxelSize();
	m_time = 0.0f;
	
	// define global forces
	m_gravity = Vector3(0.0, -9.8, 0.0);
	
	// allocate memory
	m_tempW.resize(m_numPoints);
	m_tempP.resize(m_numPoints);
	m_tempZ.resize(m_numPoints);
	m_tempR.resize(m_numPoints);
	m_tempQ.resize(m_numPoints);
	m_curl.resize(m_numPoints);
	m_vorticityConfinementForce.resize(m_numPoints);
	m_curlMagnitude.resize(m_numPoints);
	m_divergence.resize(m_numPoints);
	m_pressure.resize(m_numPoints);
	m_precond.resize(m_numPoints);
	m_ADiag.resize(m_numPoints);
	m_APlusX.resize(m_numPoints);
	m_APlusY.resize(m_numPoints);
	m_APlusZ.resize(m_numPoints);
	m_time = 0.0f;
	
	// initialize matrix/preconditioner
	constructMatrix(m_dx, 0.1f);
	constructPreconditioner();
	
	// initialize CG stopping criterion
	setCGTolerance( sqrt( FLT_EPSILON ) );
	
}


FluidSolver::~FluidSolver()
{
}

/**
 * Sets the tolerance for the residual at step n, \f$|r_n|_2 \leq \mathrm{tol}\f$,
 * in the Euclidean 2-norm.
 */
void FluidSolver::setCGTolerance(float tol) { tol_cg = tol; }

/**
 * Sets the maximum number of iterations for the conjugate gradients solver(s).
 */
void FluidSolver::setCGMaxIter(unsigned N) { maxiter_cg = N; }
	
/**
 * Uses the formulation of the maximum timestep from Foster and Fedkiw '01.<br/>
 * See: <a href="http://physbam.stanford.edu/~fedkiw/papers/stanford2001-02.pdf">Nick Foster Ronald Fedkiw. Practical Animation of Liquids. SIGGRAPH, pages 23-30, 2001.</a>
 *
 * @return the max timestep float for deltaT.
 */
float FluidSolver::computeMaxTimeStep() const
{
	// TODO: Consider all body forces, including gravity
	// See: Foster and Fedkiw 01
	fdl::Vector3 g = m_gravity;
	g.x *= (5.0 * m_dx);
	g.y *= (5.0 * m_dx);
	g.z *= (5.0 * m_dx);
	
	float maxVelocity = grid->getMaximumVelocity();
	if (maxVelocity < 0.2) maxVelocity = 0.2f;
	
	maxVelocity += m_gravity.getSquareRoot().length();
	return (5.0*m_dx) / maxVelocity;
}


/**
 * This is the main function for updating the fluid system. Input of zero or less 
 * (which is the default) will trigger a call to compute the maximum time step 
 * using computeMaxTimeStep.
 *
 * @param dt delta time value to step forward
 */
void FluidSolver::step(float dt)
{	
	if(dt <= 0){
		dt = computeMaxTimeStep();
	}else{
		dt = std::max(dt, computeMaxTimeStep());
	}
	INFO() << "FluidSolver: step ("  << dt << ") ..";
	
	INFO() << "  + Performing advection step ..";
	advect(dt);
	
	INFO() << "  + Apply forces ..";
	applyForces(dt);

	INFO() << "  + Performing projection step ..";
	project(dt);


	// add density
	float targetTemp = 28.5;
	float rateDensity = 100.0f;
	INFO() << "  + Updating densities ..";
	int rad = 4;
	Sample cell(dt*rateDensity,targetTemp,0);
//	for (int z=0; z<m_gridZ; ++z) {
//		int pos = z * m_slice;
//		int velIdx = z * m_velSlice;
//		for (int y=0; y<m_gridY; ++y, ++velIdx) {
//			for (int x=0; x<m_gridX; ++x, ++velIdx, ++pos) {
	for (int i=-rad; i<=rad; ++i) {
		for (int j=-rad; j<=rad; ++j) {
			int x = (int) m_gridX/2 + j;
			int y = (int) m_gridY/2 + i;
			int z = (int) m_gridZ/2;
			float tmp = rad*rad - i*i - j*j;
			if (tmp < 0){
				continue;
			}
			int dist = -std::max((int) std::sqrt(tmp)-1, 0);
			grid->setDensity(dist + x + y * m_gridX + z * m_slice, cell);
			grid->getForce(1)[dist + x + y * (m_gridX+1) + z * m_velSlice] = 0.3 * m_gravity.y;
		}
	}

	// apply forces to velocity field
	INFO() << "  + Adding forces ..";
	for (int z=0; z<m_gridZ; ++z) {
		int velIdx = z * m_velSlice;
		for (int y=0; y<m_gridY; ++y, ++velIdx) {
			for (int x=0; x<m_gridX; ++x, ++velIdx) {
				grid->getVelocity(0)[velIdx] += dt * grid->getForce(0)[velIdx];
				grid->getVelocity(1)[velIdx] += dt * grid->getForce(1)[velIdx];
				grid->getVelocity(2)[velIdx] += dt * grid->getForce(2)[velIdx];
			}
		}
	}
	
	// diffuse(dt, 1.5);
	
	m_time += dt;
}

/**
 * Computes and applies all body forces to the system. Currently this includes vorticity
 * confinement and bouyancy. More to come...
 *
 * @param dt delta time value to step forward
 *
 */
void FluidSolver::applyForces(float dt)
{
	grid->clearForces();

	// Calculate the magnitude of the curl
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z * m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (grid->isSolid(pos))
					continue;
					
				fdl::Vector3 velTop     = grid->getVelocity(fdl::Point3((x+0.5f)*m_dx, (y+0.0f)*m_dx, (z+0.5f)*m_dx));
				fdl::Vector3 velBot     = grid->getVelocity(fdl::Point3((x+0.5f)*m_dx, (y+1.0f)*m_dx, (z+0.5f)*m_dx));
				fdl::Vector3 velLeft    = grid->getVelocity(fdl::Point3((x+0.0f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx));
				fdl::Vector3 velRight   = grid->getVelocity(fdl::Point3((x+1.0f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx));
				fdl::Vector3 velFront   = grid->getVelocity(fdl::Point3((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+0.0f)*m_dx));
				fdl::Vector3 velBack    = grid->getVelocity(fdl::Point3((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+1.0f)*m_dx));

				fdl::Vector3 dudx = (velRight - velLeft) / m_dx;
				fdl::Vector3 dudy = (velBot - velTop) / m_dx;
				fdl::Vector3 dudz = (velBack - velFront) / m_dx;
				fdl::Vector3 curl(dudy.z - dudz.y, dudz.x - dudx.z, dudx.y-dudy.x);
				m_curl[pos] = curl;
				m_curlMagnitude[pos] = curl.length();
			}
		}
	}

	// Add the vorticity confinement force
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z * m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (y==0 || x==0 || z == 0 || x==m_gridX-1 
					|| y==m_gridY-1 || z == m_gridZ-1 || grid->isSolid(pos)) {
					m_curl[pos] = fdl::Vector3(0,0,0);
					continue;
				}
				fdl::Vector3 N(
					m_curlMagnitude[pos+1] - m_curlMagnitude[pos-1],
					m_curlMagnitude[pos+m_gridX] - m_curlMagnitude[pos-m_gridX],
					m_curlMagnitude[pos+m_slice] - m_curlMagnitude[pos-m_slice]
				);
				float length = N.length();
				if (length < EPSILON) continue;
				
				m_vorticityConfinementForce[pos] = cross(N / length, m_curl[pos]) * (m_dx * 0.8f);
			}
		}
	}

	// Add the buoyancy force 
	float ambient = 0;
	float a = 0.0625f*0.5f;
	float b = 0.025f;
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z * m_slice;
		int velIdx = z * m_velSlice;
		for (int y=0; y<m_gridY; ++y, ++velIdx) {
			for (int x=0; x<m_gridX; ++x, ++velIdx, ++pos) {
				if (grid->isSolid(pos))
					continue;

				const Sample top = grid->getDensity((x+0.5f)*m_dx, y*m_dx, (z+0.5f)*m_dx);

				if (x != 0 && !grid->isSolid(pos-1)) 
					grid->getForce(0)[velIdx] += (m_vorticityConfinementForce[pos].x + m_vorticityConfinementForce[pos-1].x) * 0.5f;

				if (y!= 0 && !grid->isSolid(pos-m_gridX)) {
					grid->getForce(1)[velIdx] -= -a*top.density + b*(top.temperature - ambient);
					grid->getForce(1)[velIdx] += (m_vorticityConfinementForce[pos].y + m_vorticityConfinementForce[pos-m_gridX].y) * 0.5f;
				}

				if (z != 0 && !grid->isSolid(pos-m_slice)) 
					grid->getForce(2)[velIdx] += (m_vorticityConfinementForce[pos].z + m_vorticityConfinementForce[pos-m_slice].z) * 0.5f;
			}
		}
	}
}


/**
 * Takes the fluid field update and projects it down to a divergence-free field using
 * a linear solve on the new pressure values. Afterwards, the solved pressures will be
 * applied to the grid.
 *
 * @param dt delta time value to step forward
 *
 */
void FluidSolver::project(float dt)
{
	// Calculate the divergence 
	// Here we use rho if the entire medium has constant density.
	// Otherwise we must add it to the pressure update.
	float rho = 0.25f;
	float inv_dx = 1.0 / m_dx;
	
	// CONSTANT DENSITY
	for (int z=0; z<m_gridZ; ++z) {
		int velIdx = z * m_velSlice, pos = z * m_slice;
		for (int y=0; y<m_gridY; ++y, ++velIdx) {
			for (int x=0; x<m_gridX; ++x, ++pos, ++velIdx) {
				if (grid->isSolid(pos)) {
					m_divergence[pos] = 0;
					continue;
				}
				m_divergence[pos] = 
					( grid->getVelocity(0)[velIdx+1]			- grid->getVelocity(0)[velIdx]
					+ grid->getVelocity(1)[velIdx+m_gridX+1]	- grid->getVelocity(1)[velIdx]
					+ grid->getVelocity(2)[velIdx+m_velSlice]	- grid->getVelocity(2)[velIdx])
					* inv_dx;
			}
		}
	}
	
	// VARIABLE DENSITY
	/*
	float t1,t2;
	float v1,v2,v3,v4,v5,v6;
	for (int z=0; z<m_gridZ; ++z) {
		int velIdx = z * m_velSlice, pos = z * m_slice;
		for (int y=0; y<m_gridY; ++y, ++velIdx) {
			for (int x=0; x<m_gridX; ++x, ++pos, ++velIdx) {
				if (grid->isSolid(pos)) {
					m_divergence[pos] = 0;
					continue;
				}
				
				t1 = t2 = 0;
				float scale = 1.0 / m_dx;
				
				v1 = grid->isSolid(pos+1)? 1: 0;
				v2 = grid->isSolid(pos-1)? 1: 0;
				v3 = grid->isSolid(pos+m_gridX)? 1: 0;
				v4 = grid->isSolid(pos-m_gridX)? 1: 0;
				v5 = grid->isSolid(pos+m_slice)? 1: 0;
				v6 = grid->isSolid(pos-m_slice)? 1: 0;
				
				t1 =	(v1*grid->getVelocity(0)[velIdx+1]			- v2*grid->getVelocity(0)[velIdx]) +
						(v3*grid->getVelocity(1)[velIdx+m_gridX+1]	- v4*grid->getVelocity(1)[velIdx]) +
						(v5*grid->getVelocity(2)[velIdx+m_velSlice]	- v6*grid->getVelocity(2)[velIdx]);

				// t2 =	(v1*grid->getVelocity(0)[velIdx+1]			- v2*grid->getVelocity(0)[velIdx]) +
				// 		(v3*grid->getVelocity(1)[velIdx+m_gridX+1]	- v4*grid->getVelocity(1)[velIdx]) +
				// 		(v5*grid->getVelocity(2)[velIdx+m_velSlice]	- v6*grid->getVelocity(2)[velIdx]);
				m_divergence[pos] = -scale*t1 + scale*t2;
				
				if(m_divergence[pos] != m_divergence[pos]){
					std::cerr << "m_divergence[" << pos << "] = nan" << std::endl;
					exit(1);
				}
			}
		}
	}
	*/
	
	// TO-DO: Include smoke divergence control
	
	// Construct new coefficient matrix
	//constructMatrix(m_dx, dt, rho, true);
	constructMatrix(m_dx, dt, rho);
	
	// Re-construct the preconditioner 
	constructPreconditioner(rho);

	// Perform linear solve
	pcgSolve(m_divergence, m_pressure);

	// Apply the computed pressure gradients [CONSTANT DENSITY]
	float scale = dt / (rho * m_dx);
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z * m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (grid->isSolid(pos)) 
					continue;
				int velIdx = x + y * (m_gridX+1) + z * m_velSlice;
				if (x < m_gridX-1 && !grid->isSolid(pos + 1)) {
					grid->getLastVelocity(0)[velIdx+1] = 
						grid->getVelocity(0)[velIdx+1] + (m_pressure[pos+1] - m_pressure[pos]) * scale;
				}
				if (y < m_gridY-1 && !grid->isSolid(pos + m_gridX)) {
					grid->getLastVelocity(1)[velIdx+m_gridX+1] = 
						grid->getVelocity(1)[velIdx+m_gridX+1] + (m_pressure[pos+m_gridX] - m_pressure[pos]) * scale;
				}
				if (z < m_gridZ-1 && !grid->isSolid(pos + m_slice)) {
					grid->getLastVelocity(2)[velIdx+m_velSlice] = 
						grid->getVelocity(2)[velIdx+m_velSlice] + (m_pressure[pos+m_slice] - m_pressure[pos]) * scale;
				}
			}
		}
	}
	
	// Apply the computed gradients [VARIABLE DENSITY], not ready yet. Must update A first..
	// for (int z=0; z<m_gridZ; ++z) {
	// 	int pos = z * m_slice;
	// 	for (int y=0; y<m_gridY; ++y) {
	// 		for (int x=0; x<m_gridX; ++x, ++pos) {
	// 			int velIdx = x + y * (m_gridX+1) + z * m_velSlice;
	// 			grid->setDensity(pos);
	// 			
	// 		}
	// 	}
	// }
	/*
	float scale = 0;
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z * m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (grid->isSolid(pos)) 
					continue;
					
				int velIdx = x + y * (m_gridX+1) + z * m_velSlice;
				if (x < m_gridX-1 && !grid->isSolid(pos + 1)) {
					// define scale... (with interpolated denisty value)
					scale = dt / (std::max(grid->getDensity((float)x + 0.5,y,z).density, 0.05f) * m_dx);
					grid->getLastVelocity(0)[velIdx+1] = 
						grid->getVelocity(0)[velIdx+1] + (m_pressure[pos+1] - m_pressure[pos]) * scale;
				}
				if (y < m_gridY-1 && !grid->isSolid(pos + m_gridX)) {
					scale = dt / (std::max(grid->getDensity(x,(float)y + 0.5,z).density, 0.05f) * m_dx);
					grid->getLastVelocity(1)[velIdx+m_gridX+1] = 
						grid->getVelocity(1)[velIdx+m_gridX+1] + (m_pressure[pos+m_gridX] - m_pressure[pos]) * scale;
				}
				if (z < m_gridZ-1 && !grid->isSolid(pos + m_slice)) {
					scale = dt / (std::max(grid->getDensity(x,y,(float)z + 0.5).density, 0.05f) * m_dx);
					grid->getLastVelocity(2)[velIdx+m_velSlice] = 
						grid->getVelocity(2)[velIdx+m_velSlice] + (m_pressure[pos+m_slice] - m_pressure[pos]) * scale;
				}
			}
		}
	}
	*/
	
	grid->swapVelocities();
}

/**
 * Appy diffusion by scaling the density values by an arbitrary constant rate. 
 * The calculation used here is d = e^(-dt*rate).
 *
 * @param dt delta time value to step forward
 * @param rate the rate by which we will scale down the density
 *
 */
void FluidSolver::diffuse(float dt, float rate)
{
	double exponent = -dt*rate;	// rate of density decay = e^-dt*[desired_rate]
	float scale = (float)std::pow(E,exponent);
	Sample cell;
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z * m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (grid->isSolid(pos)) 
					continue;
				
				cell.density = grid->getDensity(pos).density * scale;
				cell.smoke = grid->getDensity(pos).smoke;// * scale;	// just density for now...
				cell.temperature = grid->getDensity(pos).temperature;// * scale;	// just density for now...
				grid->setDensity(pos, cell);
			}
		}
	}
}


void FluidSolver::diffuseHeat(float dt, float rate)
{
	/* SOLVE HEAT EQUATION WITH CONVOLUTION KERNEL. */

	// Separable gaussian convolution.
	// This could be computed globally and re-used forever instead...
	std::vector<double>* kernel = gaussianHeatKernel(9, (double) rate, dt);
	
	int gauss_size = kernel->size();
	int gauss_half = (int)gauss_size/2;
	int offset, index;
	Sample cell;
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (pos >= m_numPoints) {
					break;
				}
				if (grid->isSolid(pos)) {
					continue;
				}
				
				// convolve on X
				cell = 0;
				for(int i=0; i<gauss_size; i++){
					offset = x - gauss_half + i;
					if(x >= 0 && x < m_gridX){
						index = (z * m_gridY * m_gridX) + (y * m_gridX) + offset;
						cell += grid->getDensity(index) * kernel->at(i);
					}
				}
				
				// convolve on Y
				for(int i=0; i<gauss_size; i++){
					offset = y - gauss_half + i;
					if(y >= 0 && y < m_gridY){
						index = (z * m_gridY * m_gridX) + (offset * m_gridX) + x;
						cell += grid->getDensity(index) * kernel->at(i);
					}
				}
				
				// convolve on Z
				for(int i=0; i<gauss_size; i++){
					offset = z - gauss_half + i;
					if(z >= 0 && z < m_gridZ){
						index = (offset * m_gridY * m_gridX) + (y * m_gridX) + x;
						cell += grid->getDensity(index) * kernel->at(i);
					}
				}
				
				grid->setLastDensity(pos, cell);
			}
		}
	}
	
	grid->swapDensities();
}


/**
 * NOT IMPLEMENTED YET
 *
 * @param dt delta time value to step forward
 *
 */
void FluidSolver::applyViscosity(float dt)
{
	// compute viscosity term for each fluid cell
	// for (z=0; z<m_gridZ; ++z) {
	// 	int pos = z*m_slice;
	// 	for (int y=0; y<m_gridY; ++y) {
	// 		for (int x=0; x<m_gridX; ++x, ++pos) {
	// 			int velIdx = x + y * (m_gridX+1) + z * m_velSlice;
	// 			if (grid->isSolid(pos)) 
	// 				continue;
	// 
	// 			// Advect X velocities
	// 			if (x < m_gridX-1 && !grid->isSolid(pos + 1)) {
	// 				//dt * V * laplacian
	// 				float U = (grid->getVelocity(0)[velIdx+1]			- grid->getVelocity(0)[velIdx]
	// 				+ grid->getVelocity(1)[velIdx+m_gridX+1]	- grid->getVelocity(1)[velIdx]
	// 				+ grid->getVelocity(2)[velIdx+m_velSlice]	- grid->getVelocity(2)[velIdx])
	// 				grid->getVelocity(0)[velIdx+1] = grid->getLastVelocity(p).x;
	// 			}
	// 
	// 			// Advect Y velocities
	// 			if (y < m_gridY-1 && !grid->isSolid(pos + m_gridX)) {
	// 				
	// 				grid->getVelocity(1)[velIdx+m_gridX+1] = grid->getLastVelocity(p).y;
	// 			}
	// 
	// 			// Advect Z velocities
	// 			if (z < m_gridZ-1 && !grid->isSolid(pos + m_slice)) {
	// 				
	// 				grid->getLastVelocity(2)[velIdx+m_velSlice] = grid->getLastVelocity(p).z;
	// 			}
	// 		}
	// 	}
	// }
	// 
	// grid->swapVelocities();
}

/**
 * Performs a semi-lagrangian particle back-trace from each cell center. 
 * 
 * See: <a href="http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf">Jos Stam. Stable Fluids. SIGGRAPH, pages 121–128, 1999.</a>
 *
 * @param dt delta time value to step forward
 *
 */
void FluidSolver::advect(float dt)
{
	// Advect the density field 	
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z*m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (pos >= m_numPoints) {
					break;
				}
				if (grid->isSolid(pos)) {
					continue;
				}
				fdl::Point3 p((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx);
				fdl::Vector3 vel = grid->getVelocity(p + grid->getVelocity(p) * (-dt * 0.5f)) * -dt;
				p += vel;
				grid->setLastDensity(pos, grid->getDensity(p.x, p.y, p.z));
			}
		}
	}
	
	grid->swapDensities();

	// Advect the velocity field
	for (int z=0; z<m_gridZ; ++z) {
		int pos = z*m_slice;
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				int velIdx = x + y * (m_gridX+1) + z * m_velSlice;
				if (grid->isSolid(pos)) 
					continue;

				// Advect X velocities
				if (x < m_gridX-1 && !grid->isSolid(pos + 1)) {					
					fdl::Point3 p((x+1.0f)*m_dx, (y+.5f)*m_dx, (z+0.5f)*m_dx);
					fdl::Vector3 vel = grid->getVelocity(p + grid->getVelocity(p) * (-dt * 0.5f)) * -dt;
					p += vel;
					grid->getLastVelocity(0)[velIdx+1] = grid->getVelocity(p).x;
				}

				// Advect Y velocities
				if (y < m_gridY-1 && !grid->isSolid(pos + m_gridX)) {
					fdl::Point3 p((x+0.5f)*m_dx, (y+1.0f)*m_dx, (z+0.5f)*m_dx);
					fdl::Vector3 vel = grid->getVelocity(p + grid->getVelocity(p) * (-dt * 0.5f)) * -dt;
					p += vel;
					grid->getLastVelocity(1)[velIdx+m_gridX+1] = grid->getVelocity(p).y;
				}

				// Advect Z velocities
				if (z < m_gridZ-1 && !grid->isSolid(pos + m_slice)) {
					fdl::Point3 p((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+1.0f)*m_dx);
					fdl::Vector3 vel = grid->getVelocity(p + grid->getVelocity(p) * (-dt * 0.5f)) * -dt;
					p += vel;
					grid->getLastVelocity(2)[velIdx+m_velSlice] = grid->getVelocity(p).z;
				}
			}
		}
	}
	
	grid->swapVelocities();
}

// CONSTRUCT MATRIX USING CONSTANT DENSITY
/**
 * Populates vectors that contain the coefficients of the sparse linear 
 * system that we intend to solve. The approach is described in Bridson Chapter 4.
 *
 * @param dx cell width of the grid
 * @param dt delta time value to step forward
 * @param rho the global scaling factor accounting for density
 *
 */
void FluidSolver::constructMatrix(float dx, float dt, float rho)
{
	INFO() << "    Constructing pressure matrix for constant density";
	
	// clear the coefficients matrix
	m_ADiag *= 0.0f;
	m_APlusX *= 0.0f;
	m_APlusY *= 0.0f;
	m_APlusZ *= 0.0f;
	
	const float scale = dt / (rho * dx * dx);
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				bool fluid = !grid->isSolid(pos);
				bool fluidRight = (x != m_gridX-1) && !grid->isSolid(pos+1);
				bool fluidBelow = (y != m_gridY-1) && !grid->isSolid(pos+m_gridX);
				bool fluidBehind = (z != m_gridZ-1) && !grid->isSolid(pos+m_slice);

				if (fluid && fluidRight) {
					m_ADiag[pos] += scale;
					m_ADiag[pos+1] += scale;
					m_APlusX[pos] = -scale;
				}

				if (fluid && fluidBelow) {
					m_ADiag[pos] += scale;
					m_ADiag[pos+m_gridX] += scale;
					m_APlusY[pos] = -scale;
				}

				if (fluid && fluidBehind) {
					m_ADiag[pos] += scale;
					m_ADiag[pos+m_slice] += scale;
					m_APlusZ[pos] = -scale;
				}
			}
		}
	}
}

// VARIABLE DENSITY MATRIX
/*
void FluidSolver::constructMatrix(float dx, float dt, float rho, bool variable_density)
{
	INFO() << "    Constructing pressure matrix " << (variable_density? "for variable density": "");
	
	// clear the coefficients matrix
	m_ADiag *= 0.0f;
	m_APlusX *= 0.0f;
	m_APlusY *= 0.0f;
	m_APlusZ *= 0.0f;
	
	float scale = dt / (rho * dx * dx);
	float volume, t1, t2, t3, t4;
	float d1,d2,d3,d4,d5,d6;
	float v1,v2,v3,v4,v5,v6;
	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				bool fluid = !grid->isSolid(pos);
				bool fluidRight = (x != m_gridX-1) && !grid->isSolid(pos+1);
				bool fluidBelow = (y != m_gridY-1) && !grid->isSolid(pos+m_gridX);
				bool fluidBehind = (z != m_gridZ-1) && !grid->isSolid(pos+m_slice);
				
				if(variable_density && fluid){
					// for now we make an oversimplification for volume. Either zero or one.
					v1 = grid->isSolid(pos+1)? 1: 0;
					v2 = grid->isSolid(pos-1)? 1: 0;
					v3 = grid->isSolid(pos+m_gridX)? 1: 0;
					v4 = grid->isSolid(pos-m_gridX)? 1: 0;
					v5 = grid->isSolid(pos+m_slice)? 1: 0;
					v6 = grid->isSolid(pos-m_slice)? 1: 0;
					t1 = t2 = t3 = t4 = 0;
					
					d1 = std::max(grid->getDensity((float)x+0.5,y,z).density, 0.01f);
					d2 = std::max(grid->getDensity((float)x-0.5,y,z).density, 0.01f);
					d3 = std::max(grid->getDensity(x,(float)y+0.5,z).density, 0.01f);
					d4 = std::max(grid->getDensity(x,(float)y-0.5,z).density, 0.01f);
					d5 = std::max(grid->getDensity(x,y,(float)z+0.5).density, 0.01f);
					d6 = std::max(grid->getDensity(x,y,(float)z-0.5).density, 0.01f);
					t1 = v1/d1 + v2/d2 + v3/d3 + v4/d4 + v5/d5 + v6/d6;
					t1 *= m_pressure[pos];
					
					if(pos+1 < m_numPoints)
						t2 += v1/d1 * m_pressure[pos+1];
					if(pos-1 > 0)
						t2 += v2/d2 * m_pressure[pos-1];
					
					if(pos+m_gridX < m_numPoints)
						t3 += v3/d3 * m_pressure[pos+m_gridX];
					if(pos-m_gridX > 0)
						t3 += v4/d4 * m_pressure[pos-m_gridX];
					
					if(pos+m_slice < m_numPoints)
						t4 += v5/d5 * m_pressure[pos+m_slice];
					if(pos-m_slice > 0)
						t4 += v6/d6 * m_pressure[pos-m_slice];
					
					scale = (dt / (dx * dx)) * (t1 - t2 - t3 - t4);

					if(scale != scale){
						std::cerr << "scale at " << pos << " = nan" << std::endl;
						exit(1);
					}
				}

				if (fluid && fluidRight) {
					m_ADiag[pos] += scale;
					m_ADiag[pos+1] += scale;
					m_APlusX[pos] = -scale;
				}

				if (fluid && fluidBelow) {
					m_ADiag[pos] += scale;
					m_ADiag[pos+m_gridX] += scale;
					m_APlusY[pos] = -scale;
				}

				if (fluid && fluidBehind) {
					m_ADiag[pos] += scale;
					m_ADiag[pos+m_slice] += scale;
					m_APlusZ[pos] = -scale;
				}
			}
		}
	}
}
*/

/**
 * constructPreconditioner currently makes the modified incomplete cholesky preconditioner 
 * for a preconditioned conjugate gradient solve of the positive semi-definite pressure 
 * matrix.
 *
 * @param rho the global scaling factor accounting for density
 * @param tau precondiitoner "tuning parameter"
 *
 */
void FluidSolver::constructPreconditioner(float rho, float tau)
{	
	INFO() << "    Computing modified incomplete cholesky preconditioner";
	
	float termLeft = 0.0f;
	float termAbove = 0.0f;
	float termRight = 0.0f;
	float termLeft2 = 0.0f;
	float termAbove2 = 0.0f;
	float termRight2 = 0.0f;

	for (int z=0, pos=0; z<m_gridZ; ++z) {
		for (int y=0; y<m_gridY; ++y) {
			for (int x=0; x<m_gridX; ++x, ++pos) {
				if (grid->isSolid(pos))
					continue;
					
				termLeft = termAbove = termRight = termLeft2 = termAbove2 = termRight2 = 0.0f;
				
				if (x > 0 && !grid->isSolid(pos-1)) {
					termLeft  = m_APlusX[pos-1] * m_precond[pos-1];
					termLeft2 = m_APlusX[pos-1] * (m_APlusY[pos-1] + m_APlusZ[pos-1]) * m_precond[pos-1] * m_precond[pos-1];
				}

				if (y > 0 && !grid->isSolid(pos-m_gridX)) {
					termAbove  = m_APlusY[pos-m_gridX] * m_precond[pos-m_gridX];
					termAbove2 = m_APlusY[pos-m_gridX] * (m_APlusX[pos-m_gridX] + m_APlusZ[pos-m_gridX]) * m_precond[pos-m_gridX] * m_precond[pos-m_gridX];
				}

				if (z > 0 && !grid->isSolid(pos-m_slice)) {
					termRight  = m_APlusZ[pos-m_slice] * m_precond[pos-m_slice];
					termRight2 = m_APlusZ[pos-m_slice] * (m_APlusX[pos-m_slice] + m_APlusY[pos-m_slice]) * m_precond[pos-m_slice] * m_precond[pos-m_slice];
				}
				
				termRight = termRight*termRight;
				termAbove = termAbove*termAbove;
				termLeft = termLeft*termLeft;

				float e = m_ADiag[pos] - termLeft - termAbove - termRight - tau * (termLeft2 + termAbove2 + termRight2);
				if (e < rho * m_ADiag[pos]){
					e = m_ADiag[pos];
				}
				
				m_precond[pos] = 1.0f / std::sqrt(e);
			}
		}
	}
}


/**
 * cgSolve performs a conjugate gradient solve for the pressure matrix given two vectors 
 * x and b. 
 *
 * @param b the vector representing the negative divergence
 * @param x the resulting pressure gradient after the linear solve
 *
 * @return float value representing the residual after the iterative solve
 */
float FluidSolver::cgSolve(const Vector& b, Vector& x)
{	
	int k = 0;
	float beta = 0;
	float lastRho = 0;
	float rho = 0;
	float tolerance = tol_cg * tol_cg; // rho is a "norm squared" measurement.

	axpy_prod(x, m_tempR);
	m_tempR = b - m_tempR;
	rho = inner_prod(m_tempR, m_tempR);

	while (k < maxiter_cg && rho > tolerance) {
		if (k == 0) {
			noalias(m_tempP) = m_tempR;
		} else {
			beta = rho / lastRho;
			m_tempP = m_tempR + m_tempP * beta;
		}
		
		axpy_prod(m_tempP, m_tempW);
		float alpha = rho / inner_prod(m_tempP, m_tempW);
		noalias(x) += alpha * m_tempP;
		noalias(m_tempR) -= alpha * m_tempW;
		lastRho = rho;
		rho = inner_prod(m_tempR, m_tempR);
		k++;
	}

	INFO() << "  + CG: Residual after " << k << " iterations : " << std::sqrt(rho);
	return std::sqrt(rho);
}

/**
 * pcgSolve performs a preconditioned conjugate gradient solve for the pressure matrix 
 * given two vectors x and b.
 *
 * @param b the vector representing the negative divergence
 * @param x the resulting pressure gradient after the linear solve
 *
 * @return float value representing the residual after the iterative solve
 */
float FluidSolver::pcgSolve(const Vector& b, Vector& x)
{
float M = inner_prod(x, b);
if(M!=M){
	std::cerr << "M is nan!" << std::endl;
	exit(1);
}

	int k = 0;
	float beta = 0;
	float lastRho = 0;
	float rho = 0;
	float tolerance = tol_cg * tol_cg; // rho is a "norm squared" measurement.

	axpy_prod(x, m_tempR);
	m_tempR = b - m_tempR;
	solvePreconditioner(m_tempR, m_tempZ);
	rho = inner_prod(m_tempR, m_tempZ);
	if(rho!=rho){
		std::cerr << "rho is nan!" << std::endl;
		exit(1);
	}
	
	while (k < maxiter_cg && rho > tolerance) {
		if (k == 0) {
			noalias(m_tempP) = m_tempZ;
		} else {
			beta = rho / lastRho;
			m_tempP = m_tempZ + m_tempP * beta;
		}

		axpy_prod(m_tempP, m_tempW);
		float denom = inner_prod(m_tempP, m_tempW);
		float alpha = rho / denom;
		noalias(x) += alpha * m_tempP;
		noalias(m_tempR) -= alpha * m_tempW;

		solvePreconditioner(m_tempR, m_tempZ);
		lastRho = rho;
		rho = inner_prod(m_tempR, m_tempZ);
		k++;
	}
	
	INFO() << "  + PCG: Residual after " << k << " iterations : " << std::sqrt(rho);

	return std::sqrt(rho);
}


/**
 * axpy_prod computes y = Ax. This function provides an alternative to the BLAS function 
 * axpy_prod which uses the sparse matrix A represented by vectors only. <br/>
 * See <a href="http://www.boost.org/doc/libs/1_41_0/libs/numeric/ublas/doc/products.htm">the uBlas axpy_prod method</a> for more info.
 *
 * @param x vector to be multiplied by A
 * @param y result of Ax
 *
 */
void FluidSolver::axpy_prod(const Vector& x, Vector& y) const
{
	for (int i=0; i<m_slice; ++i) {
		float result = m_ADiag[i] * x[i]
					+ m_APlusX[i] * x[i+1]
					+ m_APlusY[i] * x[i+m_gridX]
					+ m_APlusZ[i] * x[i+m_slice];
		if (i-1 >= 0){
			result += m_APlusX[i-1] * x[i-1];
		}
		if (i-m_gridX >= 0){
			result += m_APlusY[i-m_gridX] * x[i-m_gridX];
		}
		if (i-m_slice >= 0){
			result += m_APlusZ[i-m_slice] * x[i-m_slice];
		}

		y[i] = result;
	}
	
	for (int i=m_numPoints-m_slice; i<m_numPoints; ++i) {
		float result = m_ADiag[i] * x[i]
					+ m_APlusX[i-1] * x[i-1]
					+ m_APlusY[i-m_gridX] * x[i-m_gridX]
					+ m_APlusZ[i-m_slice] * x[i-m_slice];

		if (i+1 < m_numPoints){
			result += m_APlusX[i] * x[i+1];
		}
		if (i+m_gridX < m_numPoints){
			result += m_APlusY[i] * x[i+m_gridX];
		}
		if (i+m_slice < m_numPoints){
			result += m_APlusZ[i] * x[i+m_slice];
		}

		y[i] = result;
	}
	
	for (int i=m_slice; i<m_numPoints-m_slice; ++i) {
		y[i] = m_ADiag[i] * x[i]
			+ m_APlusX[i] * x[i+1]
			+ m_APlusY[i] * x[i+m_gridX]
			+ m_APlusZ[i] * x[i+m_slice]
			+ m_APlusX[i-1] * x[i-1]
			+ m_APlusY[i-m_gridX] * x[i-m_gridX]
			+ m_APlusZ[i-m_slice] * x[i-m_slice];
	}
}

/**
 * solvePreconditioner 
 *
 * @param b
 * @param x
 *
 */
void FluidSolver::solvePreconditioner(const Vector& b, Vector& x)
{
	// Solve lower triangular system
	for (int i=0; i<m_numPoints; ++i) {
		if (grid->isSolid(i))
			continue;
			
		float temp = b[i];
		if (i > 0) {
			temp -= m_APlusX[i-1] * m_precond[i-1] * m_tempQ[i-1];
		}
		if (i > m_gridX) {
			temp -= m_APlusY[i-m_gridX] * m_precond[i-m_gridX] * m_tempQ[i-m_gridX];
		}
		if (i > m_slice) {
			temp -= m_APlusZ[i-m_slice] * m_precond[i-m_slice] * m_tempQ[i-m_slice];
		}
		
		m_tempQ[i] = temp * m_precond[i];
	}

	// Solve upper triangular system
	for (int i=m_numPoints-1; i>=0; --i) {
		if (grid->isSolid(i))
			continue;
			
		float temp = m_tempQ[i];
		if (i+1 < m_numPoints) {
			temp -= m_APlusX[i] * m_precond[i] * x[i+1];
		}
		if (i+m_gridX < m_numPoints) {
			temp -= m_APlusY[i] * m_precond[i] * x[i+m_gridX];
		}
		if (i+m_slice < m_numPoints) {
			temp -= m_APlusZ[i] * m_precond[i] * x[i+m_slice];
		}
		
		x[i] = temp * m_precond[i];
	}
}

}	// namespace fdl
