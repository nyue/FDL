#ifndef __FDL_PARTICLE_SYSTEM_H
#define __FDL_PARTICLE_SYSTEM_H

#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <cfloat>
#include <cmath>

#include "core/vector.hpp"
#include "core/grid.hpp"
#include "core/particle.hpp"

namespace fdl {
	
/**
 * Maintains dynamic lists of Particle and Force objects, and provides
 * access to their state for numerical integration of dynamics.
 */
class ParticleSystem {
public:
	ParticleSystem();
	~ParticleSystem();
	
	// suggested doubling the resolution of the grid (2*2*2)
	void init(const Grid* grid, int resolution_factor=2);
	void reset();
	void step(double dt);
	unsigned long getParticleCount() { return m_numParticles; }
	
protected:
	/** Current simulation time. */
	double m_time;

	/** Array of Particle objects. */
	Particle** m_particles;
	unsigned long m_numParticles;
	const Grid* m_grid;
};

}	// namespace fdl

#endif	// __FDL_PARTICLE_SYSTEM_H