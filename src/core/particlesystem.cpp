#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <cfloat>
#include <cmath>

#include "core/grid.hpp"
#include "core/vector.hpp"
#include "core/particle.hpp"
#include "core/particlesystem.h"

namespace fdl {

ParticleSystem::ParticleSystem()
{
	m_time = 0;
}

ParticleSystem::~ParticleSystem()
{
	free(m_particles);
}

void ParticleSystem::init(const Grid* grid, int resolution_factor)
{
	m_grid = grid;
	
	if(resolution_factor < 1) resolution_factor = 1;
	
	int width = m_grid->getGridSizeX() * resolution_factor;
	int height = m_grid->getGridSizeY() * resolution_factor;
	int depth = m_grid->getGridSizeZ() * resolution_factor;
	float dx = m_grid->getVoxelSize();
	
	
	m_numParticles = width * height * depth;
	m_particles = (Particle**) malloc(m_numParticles * sizeof(Particle*));
	int index = 0;
	float pos[3];
	for(int i=0, index=0; i<width; ++i){
		for(int j=0; j<height; ++j){
			for(int k=0; k<depth; ++k, ++index){
				pos[0] = dx * i;
				pos[1] = dx * j;
				pos[2] = dx * k;
				m_particles[index] = new Particle(pos[0], pos[1], pos[2]);
			}
		}
	}
}


/** 
 *
 */
void ParticleSystem::reset()
{
	for(int i=0; i<m_numParticles; ++i){
		m_particles[i]->reset();
	}
}

/**
 * 
 */
void ParticleSystem::step(double dt)
{
	Vector3 vel;
	Vector3 pos;
	Vector3 accel;
	double half_dt = dt/2.0;
	// sample grid to get velocity and acceleration at the particle position, 
	// then update particle position using RK2 integrator.
	for(int i=0; i<m_numParticles; ++i){
		accel = m_grid->getForce((Point3) m_particles[i]->x) / m_particles[i]->mass;
		vel = m_grid->getVelocity((Point3) m_particles[i]->x) + accel * half_dt;
		pos = m_particles[i]->x + vel * half_dt;
		m_particles[i]->v = vel + accel * half_dt;
		m_particles[i]->x = pos + m_particles[i]->v * half_dt;
	}
	
	m_time += dt;
}

}	// namespace fdl