/**
 * @file grid.h
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

#ifndef __FDL_GRID_H
#define __FDL_GRID_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <math.h>

#include "core/common.h"
#include "core/vector.hpp"
#include "logger/logger.h"

namespace fdl {

struct Sample {
	float density;
	float smoke;
	float temperature;
	float volume;
	medium_T medium;

	inline Sample(float _r = 0.05f, float _s = 0.0f, float _t = 0.0f, float _v = 1.0f)
		: density(_r), smoke(_s), temperature(_t), volume(_v)  {
			medium = FLUID;
	}
	
	inline Sample operator+(const Sample &v) const
	{
		return Sample(density + v.density, smoke + v.smoke, temperature + v.temperature);
	}

	inline Sample operator-(const Sample &v) const
	{
		return Sample(density - v.density, smoke - v.smoke, temperature - v.temperature);
	}

	inline Sample& operator+=(const Sample &v)
	{
		density += v.density;
		smoke += v.smoke;
		temperature += v.temperature;
		return *this;
	}

	inline Sample& operator-=(const Sample &v)
	{
		density -= v.density;
		smoke -= v.smoke;
		temperature -= v.temperature;
		return *this;
	}

	inline Sample operator*(float f) const
	{
		return Sample(density*f, smoke*f, temperature*f);
	}

	inline Sample &operator*=(float f)
	{
		density *= f;
		smoke *= f;
		temperature *= f;
		return *this;
	}

	inline Sample operator-() const
	{
		return Sample(-density, -smoke, -temperature);
	}

	inline Sample operator/(float f) const
	{
		float r = 1.0f / f;
		return Sample(density * r, smoke * r, temperature * r);
	}

	inline Sample &operator/=(float f)
	{
		float r = 1.0f / f;
		density *= r;
		smoke *= r;
		temperature *= r;
		return *this;
	}

	inline bool operator==(const Sample &v) const
	{
		return (v.density == density && v.smoke == smoke && v.temperature == temperature && v.volume == volume);
	}

	inline bool operator!=(const Sample &v) const
	{
		return !operator==(v);
	}

	inline std::string toString() const
	{
		std::ostringstream oss;
		std::string mName = "";
		if(medium==SOLID) mName = "solid";
		if(medium==SMOKE) mName = "smoke";
		if(medium==FLUID) mName = "fluid";
		if(medium==AIR) mName = "air";
		oss << "[" << medium << ": " << density << ", " << smoke << ", " << temperature << ", " << volume << "]";
		return oss.str();
	}
};

template<class T>
class TGrid {
public:
	TGrid(int _x, int _y, int _z, float _dx) : m_gridX(_x), m_gridY(_y), m_gridZ(_z), m_dx(_dx)
	{
		int velocityPoints = (m_gridX+1) * (m_gridY+1) * (m_gridZ+1);
		m_numPoints = m_gridX * m_gridY * m_gridZ;
		m_slice = m_gridX * m_gridY;
		m_velSlice = (m_gridX+1) * (m_gridY+1);

		std::cout << "FluidSolver: Allocating " << (sizeof(float)*(m_numPoints*15 + velocityPoints*DIMENSIONS*3)
			+ m_numPoints*2*sizeof(double))/1024 << " KB for a " << m_gridX << "x" << m_gridY << "x" << m_gridZ << " MAC grid.." << std::endl;
		m_invDx = 1.0f / m_dx;

		/* Allocate dense vectors for all quantities, which are stored on the MAC grid. */
		m_d0.resize(m_numPoints);
		m_d1.resize(m_numPoints);
		for (int i=0; i<DIMENSIONS; ++i) {
			m_forces[i].resize(velocityPoints);
			m_u0[i].resize(velocityPoints);
			m_u1[i].resize(velocityPoints);
		}
		// m_solid = new bool[m_numPoints];
		// memset(m_solid, 0, sizeof(bool)*m_numPoints);

		// seed block of densities
		// int top_limit = (int) m_gridY/2;
		// for (int z=0; z<m_gridZ; ++z) {
		// 	for (int y=0; y<top_limit; ++y) {
		// 		for (int x=0; x<m_gridX; ++x) {
		// 			int pos = x + (y * m_gridX) + (z * m_slice);
		// 			m_d0[pos] = (double) z/m_gridZ;
		// 		}
		// 	}
		// }
	}

	~TGrid()
	{
	}


	/**
	 * Allocates a copy of the grid and returns it.
	 *
	 * @return a copy of the grid
	 */
	TGrid<T>* copy()
	{
		TGrid<T>* _g = new TGrid<T>(m_gridX, m_gridY, m_gridZ, m_dx);
		
		for(int i=0; i<DIMENSIONS; i++){
			_g->m_u0[i] = this->getVelocity(i);
			_g->m_u1[i] = this->getLastVelocity(i);
			_g->m_forces[i] = this->getForce(i);
		}
		_g->m_d0 = this->getDensity();
		_g->m_d1 = this->getLastDensity();
		return _g;
	}


	/**
	 * Swaps the internal density vectors
	 *
	 */
	void swapDensities()
	{
		m_d0.swap(m_d1);
	}

	/**
	 * Swaps the internal velocity vectors
	 *
	 */
	void swapVelocities()
	{
		m_u0[0].swap(m_u1[0]);
		m_u0[1].swap(m_u1[1]);
		m_u0[2].swap(m_u1[2]);
	}

	/**
	 * Sets the force field vectors to zero.
	 *
	 */
	void clearForces()
	{
		std::fill(m_forces[0].begin(), m_forces[0].end(), 0);
		std::fill(m_forces[1].begin(), m_forces[1].end(), 0);
		std::fill(m_forces[2].begin(), m_forces[2].end(), 0);
	}

	/**
	 * Sets the velocity field vectors to zero.
	 *
	 */
	void clearVelocities()
	{
		std::fill(m_u0[0].begin(), m_u0[0].end(), 0);
		std::fill(m_u0[1].begin(), m_u0[1].end(), 0);
		std::fill(m_u0[2].begin(), m_u0[2].end(), 0);

		std::fill(m_u1[0].begin(), m_u1[0].end(), 0);
		std::fill(m_u1[1].begin(), m_u1[1].end(), 0);
		std::fill(m_u1[2].begin(), m_u1[2].end(), 0);
	}

	/**
	 * Sets teh density field vectors to zero.
	 *
	 */
	void clearDensities()
	{
		std::fill(m_d0.begin(), m_d0.end(), 0);
		std::fill(m_d1.begin(), m_d1.end(), 0);
	}


	/**
	 * Updates a value of the force field at a position given a dimension.
	 *
	 * @param dimension the array index for the dimension
	 * @param index the array index within the field
	 * @param value the new force value
	 *
	 */
	void setForce(int dimension, int index, float value)
	{
		m_forces[dimension][index] = value;
	}

	/**
	 * Updates the centered grid-cell value for density at a position (using m_d0).
	 *
	 * @param index the array position to update
	 * @param value the new Sample value
	 *
	 */
	void setDensity(int index, T value)
	{
		m_d0[index] = value;
	}

	/**
	 * Updates the centered grid-cell value for density within the "old" density vector (m_d1).
	 *
	 * @param index the array position to update
	 * @param value the new Sample value
	 *
	 */
	void setLastDensity(int index, T value)
	{
		m_d1[index] = value;
	}

	/**
	 * Samples the velocity field at a point.
	 *
	 * @param p the point to sample velocity from
	 * 
	 * @return the velocity vector
	 *
	 */
	fdl::Vector3 getVelocity(fdl::Point3 p) const
	{
		/* Re-scale to [0, gridX] x [0,gridY] x [0,gridZ] */
		p.x /= m_dx;
		p.y /= m_dx;
		p.z /= m_dx;

		return fdl::Vector3(
			getVelocityComponent(p.x, p.y - 0.5f, p.z - 0.5f, 0),
			getVelocityComponent(p.x - 0.5f, p.y, p.z - 0.5f, 1),
			getVelocityComponent(p.x - 0.5f, p.y - 0.5f, p.z, 2)
		);
	}


	/**
	 * Samples the force field at a point.
	 *
	 * @param p the point to sample force from
	 * 
	 * @return the force vector
	 *
	 */
	fdl::Vector3 getForce(fdl::Point3 p) const
	{
		/* Re-scale to [0, gridX] x [0,gridY] x [0,gridZ] */
		p.x /= m_dx;
		p.y /= m_dx;
		p.z /= m_dx;

		return fdl::Vector3(
			getForceComponent(p.x, p.y - 0.5f, p.z - 0.5f, 0),
			getForceComponent(p.x - 0.5f, p.y, p.z - 0.5f, 1),
			getForceComponent(p.x - 0.5f, p.y - 0.5f, p.z, 2)
		);
	}


	/**
	 * Samples the velocity field at a point using linear interpolation.
	 *
	 * @param x the x coordinate to sample at
	 * @param y the y coordinate to sample at
	 * @param z the z coordinate to sample at
	 * @param c the array index for the dimension of the velocity vector
	 * 
	 * @return the computed velocity
	 *
	 */
	float getVelocityComponent(float x, float y, float z, int c) const
	{
		int i = (int) x;
		int j = (int) y;
		int k = (int) z;
		// TO-DO: bounds checking!!!
		int pos = i + j * (m_gridX+1) + k*m_velSlice;

		/* Return zero for positions outside of the grid */
		if (i < 0 || j < 0 || k < 0 || i >= m_gridX || j >= m_gridY || k >= m_gridZ){
			return 0;
		}

		float alpha = x-i;
		float beta = y-j;
		float gamma = z-k;
		float A1 = m_u0[c][pos];
		float B1 = m_u0[c][pos + 1];
		float C1 = m_u0[c][pos + m_gridX + 1];
		float D1 = m_u0[c][pos + m_gridX + 2];
		float A2 = m_u0[c][pos + m_velSlice];
		float B2 = m_u0[c][pos + m_velSlice + 1];
		float C2 = m_u0[c][pos + m_velSlice + m_gridX + 1];
		float D2 = m_u0[c][pos + m_velSlice + m_gridX + 2];

		return (1-gamma) * ((1-alpha) * (1-beta) * A1 + alpha * (1-beta) * B1 + (1-alpha) * beta * C1 + alpha*beta*D1)
			 	+ gamma * ((1-alpha) * (1-beta) * A2 + alpha * (1-beta) * B2 + (1-alpha) * beta * C2 + alpha*beta*D2);
	}


	/**
	 * Samples the force field at a point using linear interpolation.
	 *
	 * @param x the x coordinate to sample at
	 * @param y the y coordinate to sample at
	 * @param z the z coordinate to sample at
	 * @param c the array index for the dimension of the force vector
	 * 
	 * @return the computed force
	 *
	 */
	float getForceComponent(float x, float y, float z, int c) const
	{
		int i = (int) x;
		int j = (int) y;
		int k = (int) z;
		// TO-DO: bounds checking!!!
		int pos = i + j * (m_gridX+1) + k*m_velSlice;

		/* Return zero for positions outside of the grid */
		if (i < 0 || j < 0 || k < 0 || i >= m_gridX || j >= m_gridY || k >= m_gridZ){
			return 0;
		}

		float alpha = x-i;
		float beta = y-j;
		float gamma = z-k;
		float A1 = m_forces[c][pos];
		float B1 = m_forces[c][pos + 1];
		float C1 = m_forces[c][pos + m_gridX + 1];
		float D1 = m_forces[c][pos + m_gridX + 2];
		float A2 = m_forces[c][pos + m_velSlice];
		float B2 = m_forces[c][pos + m_velSlice + 1];
		float C2 = m_forces[c][pos + m_velSlice + m_gridX + 1];
		float D2 = m_forces[c][pos + m_velSlice + m_gridX + 2];

		return (1-gamma) * ((1-alpha) * (1-beta) * A1 + alpha * (1-beta) * B1 + (1-alpha) * beta * C1 + alpha*beta*D1)
			 	+ gamma * ((1-alpha) * (1-beta) * A2 + alpha * (1-beta) * B2 + (1-alpha) * beta * C2 + alpha*beta*D2);
	}


	/**
	 * Samples the density field (stored at cell centers) at a point using linear interpolation.
	 *
	 * @param x the x coordinate to sample at
	 * @param y the y coordinate to sample at
	 * @param z the z coordinate to sample at
	 * @param c the array index for the dimension of the density vector
	 * 
	 * @return the computed Sample value
	 *
	 */
	T getDensity(float x, float y, float z) const
	{
		// return getDensityTrilinear(x, y, z);
		return getDensityCatmullRom(x,y,z);
	}

	// bad behavior...
	float* getDensityArray() const
	{
		float* _density = (float*)malloc(sizeof(float) * m_numPoints);
		for (int i=0; i<m_numPoints; ++i) {
			_density[i] = (float) m_d0[i].density;
		}
		return _density;
	}

	/**
	 * Samples the density field (stored at cell centers) at a point using linear interpolation.
	 *
	 * @param x the x coordinate to sample at
	 * @param y the y coordinate to sample at
	 * @param z the z coordinate to sample at
	 * 
	 * @return the computed Sample value
	 *
	 */
	T getDensityTrilinear(float x, float y, float z) const
	{
		x = (x / m_dx) - 0.5f;
		y = (y / m_dx) - 0.5f;
		z = (z / m_dx) - 0.5f;

		int i = (int) x;
		int j = (int) y;
		int k = (int) z;
		int pos = i + (j * m_gridX) + (k * m_slice);

		if (i < 0 || j < 0 || k < 0 || i > m_gridX-1 || j > m_gridY-1 || k > m_gridZ-1){
			return 0;
		}

		float alpha = x-i;
		float beta = y-j;
		float gamma = z-k;

		T A1 = m_d0[pos];
		T B1 = (i+1<m_gridX) ? m_d0[pos+1] : 0;
		T C1 = (j+1<m_gridY) ? m_d0[pos+m_gridX] : 0;
		T D1 = (i+1<m_gridX && j+1<m_gridY) ? m_d0[pos+m_gridX+1] : 0;

		T A2, B2, C2, D2;
		if (k + 1 < m_gridZ) {
			A2 = m_d0[pos+m_slice];
			B2 = (i+1<m_gridX) ? m_d0[pos+1+m_slice] : 0;
			C2 = (j+1<m_gridY) ? m_d0[pos+m_gridX+m_slice] : 0;
			D2 = (i+1<m_gridX && j+1<m_gridY) ? m_d0[pos+m_gridX+m_slice+1] : 0;
		}

		return  (A1 * ((1-alpha) * (1-beta))
			  +  B1 * (   alpha  * (1-beta))
			  +  C1 * ((1-alpha) *    beta)
			  +  D1 *     alpha  *    beta) * (1-gamma)
		      + (A2 * ((1-alpha) * (1-beta))
			  +  B2 * (   alpha  * (1-beta))
			  +  C2 * ((1-alpha) *    beta)
			  +  D2 *     alpha  *    beta) * gamma;
	}


	/**
	 * Samples the density field (stored at cell centers) at a point using a Catmull-Rom 
	 * spline interpolation. We start with the z dimension and then use catmullRomX and 
	 * catmullRomY for x and y dimensions respectively.
	 *
	 * @param x the x coordinate to sample at
	 * @param y the y coordinate to sample at
	 * @param z the z coordinate to sample at
	 * 
	 * @return the computed Sample value
	 *
	 * @see TGrid::catmullRomY
	 * @see TGrid::catmullRomX
	 */
	T getDensityCatmullRom(float x, float y, float z) const
	{
		x = (x / m_dx) - 0.5f;
		y = (y / m_dx) - 0.5f;
		z = (z / m_dx) - 0.5f;

		int i = (int) x;
		int j = (int) y;
		int k = (int) z;

		if (i < 0 || j < 0 || k < 0 || i > m_gridX-1 || j > m_gridY-1 || k > m_gridZ-1){
			return 0;
		}

		float alpha = x-i;
		float beta = y-j;
		float gamma = z-k;

		// bounds checking is here...
		T A = (z-1>= 0)? catmullRomY(x, y, z-1, alpha, beta): 0;
		T B = catmullRomY(x, y, z, alpha, beta);
		T C = (z+1<m_gridZ)? catmullRomY(x, y, z+1, alpha, beta): 0;
		T D = (z+2<m_gridZ)? catmullRomY(x, y, z+2, alpha, beta): 0;

		float gamma2 = gamma*gamma;
		float gamma3 = gamma2*gamma;
		T d =	A * (-0.5f*gamma + gamma2 - 0.5f*gamma3) + B * (1.0f - gamma2*(5.0f/2.0f) + gamma3*(3.0f/2.0f)) +
					C * (0.5f*gamma + 2*gamma2 - gamma3*(3.0f/2.0f)) + D * (-0.5f*gamma2 + 0.5f*gamma3);

		/* Switch to trilinear interpolation in the case of an overshoot */
		if (d.density < std::min(B.density, C.density) || d.temperature < std::min(B.temperature, C.temperature) 
			|| d.density > std::max(B.density, C.density) || d.temperature > std::max(B.temperature, C.temperature)
			|| d.smoke < std::min(B.smoke, C.smoke) || d.smoke > std::max(B.smoke, C.smoke)) {
			return B*(1.0f - gamma) + C*gamma;
		}

		return d;
	}


	/**
	 * Finds the greatest velocity value in the field and returns it. This is used by the 
	 * FluidSolver to compute the maximum delta time value.
	 * 
	 * @return the float velocity value
	 *
	 * @see FluidSolver::computeMaxTimeStep
	 *
	 */
	float getMaximumVelocity() const
	{
		float maxVelocity = 0.0f;
		for (int z=0; z<m_gridZ; ++z) {
			for (int y=0; y<m_gridY; ++y) {
				for (int x=0; x<m_gridX; ++x) {
					fdl::Point3 p((x+0.5f)*m_dx, (y+0.5f)*m_dx, (z+0.5f)*m_dx);
					float length = getVelocity(p).length();
					if (length > maxVelocity)
						maxVelocity = length;
				}
			}
		}

		return maxVelocity;
	}
	
	int getGridSizeX() const { return m_gridX; }
	int getGridSizeY() const { return m_gridY; }
	int getGridSizeZ() const { return m_gridZ; }
	int getNumberOfGridCells() const { return m_numPoints; }
	int getDensityGridSlice() const { return m_slice; }
	int getVelocityGridSlice() const { return m_velSlice; }
	float getVoxelSize() const { return m_dx; }
	// const Vector& getVelocity(int i) const { return m_u0[i]; }
	// const Vector& getLastVelocity(int i) const { return m_u1[i]; }
	// const Vector& getForce(int i) const { return m_forces[i]; }
	// const std::vector<double>& getDensity() const { return m_d0; }
	
	Vector& getVelocity(int dimension) { return m_u0[dimension]; }
	Vector& getLastVelocity(int dimension) { return m_u1[dimension]; }
	Vector& getForce(int dimension) { return m_forces[dimension]; }
	const bool isSolid(int index) const { return m_d0[index].medium==SOLID; }
	const bool isFluid(int index) const { return m_d0[index].medium==FLUID; }
	const bool isSmoke(int index) const { return m_d0[index].medium==SMOKE; }
	const bool isAir(int index) const { return m_d0[index].medium==AIR; }
	const T& getDensity(int index) const { return m_d0[index]; }
	std::vector<T>& getDensity() { return m_d0; }
	std::vector<T>& getLastDensity() { return m_d1; }
	
	void setVelocityX(int index, float value) { m_u0[0][index] = value; }
	void setVelocityY(int index, float value) { m_u0[1][index] = value; }
	void setVelocityZ(int index, float value) { m_u0[2][index] = value; }
	void setVelocity(int dimension, int index, float value) { m_u0[dimension][index] = value; }
	void setLastVelocity(int dimension, int index, float value) { m_u1[dimension][index] = value; }
	
private:
	/* Cell centers - density, temperature, etc. */
	std::vector<T> m_d0;
	std::vector<T> m_d1;

	/* Velocities */
	Vector m_u0[DIMENSIONS];
	Vector m_u1[DIMENSIONS];
	
	/* Aggregated forces */
	Vector m_forces[DIMENSIONS];

	/* Grid resolution */
	int m_gridX;
	int m_gridY;
	int m_gridZ;
	int m_numPoints;
	int m_slice;
	int m_velSlice;
	float m_dx;
	float m_invDx;

	/**
	 * Returns a density sample using catmullRomX.
	 *
	 * @param x the x coordinate to sample at
	 * @param y the y coordinate to sample at
	 * @param z the z coordinate to sample at
	 * @param alpha the interpolation first weight
	 * @param beta the interpolation second weight
	 * 
	 * @return the computed Sample value
	 *
	 * @see Grid::catmullRomX
	 * @see Grid::getDensityCatmullRom
	 *
	 */
	inline T catmullRomY(int x, int y, int z, float alpha, float beta) const
	{		
		T A = (y-1>=0)? catmullRomX(x, y-1, z, alpha): 0;
		T B = catmullRomX(x, y, z, alpha);
		T C = (y+1<m_gridY)? catmullRomX(x, y+1, z, alpha): 0;
		T D = (y+2<m_gridY)? catmullRomX(x, y+2, z, alpha): 0;

		float beta2 = beta*beta;
		float beta3 = beta2*beta;
		T d =	A * (-0.5f*beta + beta2 - 0.5f*beta3) + B * (1.0f - beta2*(5.0f/2.0f) + beta3*(3.0f/2.0f)) + 
					C * (0.5f*beta + 2*beta2 - beta3*(3.0f/2.0f)) + D * (-0.5f*beta2 + 0.5f*beta3);

		/* Switch to trilinear interpolation in the case of an overshoot */
		if (d.density < std::min(B.density, C.density) || d.temperature < std::min(B.temperature, C.temperature)
			|| d.density > std::max(B.density, C.density) || d.temperature > std::max(B.temperature, C.temperature)
			|| d.smoke < std::min(B.smoke, C.smoke) || d.smoke > std::max(B.smoke, C.smoke)) {
			return B*(1.0f - beta) + C*beta;
		}

		return d;
	}

	/**
	 * Returns a density sample using m_d0 and the alpha interpolation weight.
	 *
	 * @param x the x coordinate to sample at
	 * @param y the y coordinate to sample at
	 * @param z the z coordinate to sample at
	 * @param alpha the interpolation weight
	 * 
	 * @return the computed Sample value
	 *
	 * @see Grid::catmullRomY
	 * @see Grid::getDensityCatmullRom
	 *
	 */
	inline T catmullRomX(int x, int y, int z, float alpha) const
	{
		int pos = x + y*m_gridX + z*m_slice;
		T A = (x-1 >= 0)? m_d0[pos-1]: 0;
		T B = m_d0[pos];
		T C = (x+1<m_gridX)? m_d0[pos+1]: 0;
	    T D = (x+2<m_gridX)? m_d0[pos+2]: 0;

		float alpha2 = alpha*alpha;
		float alpha3 = alpha2*alpha;
		T d =	A * (-0.5f*alpha + alpha2 - 0.5f*alpha3) + B * (1.0f - alpha2*(5.0f/2.0f) + alpha3*(3.0f/2.0f)) + 
					C * (0.5f*alpha + 2*alpha2 - alpha3*(3.0f/2.0f)) + D * (-0.5f*alpha2 + 0.5f*alpha3);

		/* Switch to trilinear interpolation in the case of an overshoot */
		if (d.density < std::min(B.density, C.density) || d.temperature < std::min(B.temperature, C.temperature)
			|| d.density > std::max(B.density, C.density) || d.temperature > std::max(B.temperature, C.temperature)
			|| d.smoke < std::min(B.smoke, C.smoke) || d.smoke > std::max(B.smoke, C.smoke)) {
			return B*(1.0f - alpha) + C*alpha;
		}

		return d;
	}	
	
};	// class Grid

typedef TGrid<Sample> Grid;
typedef TGrid<float> Gridf;
typedef TGrid<double> Gridd;

}	// namespace fdl

#endif // __FDL_GRID_H