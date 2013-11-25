#ifndef __FDL_PARTICLE_H
#define __FDL_PARTICLE_H

#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <cfloat>
#include <cmath>

#include "core/vector.hpp"

#define PARTICLE_MASS 1.0;

namespace fdl {

template<typename scalar_t>
class TParticle {
public:
	/** Default mass. */
	scalar_t mass;

	/** Deformed Position. */
	TVec3<scalar_t> x;

	/** Velocity. We probably don't need this actually. */
	TVec3<scalar_t> v;

	/** Force accumulator. We probably don't need this actually.  */
	TVec3<scalar_t> f;

	/** particle color. */
	TVec3<scalar_t> c;

	/** Density. Necessary? */
	scalar_t rho;
	
	/** OpenGL display list name identifier */
	// unsigned int name;
	
	/** 
	 * Constructs particle with default data.
	 * 
	 */
	TParticle() 
	{
		timeToLive = 1.0;
		x *= 0;
		v *= 0;
		f *= 0;
		// name = -1;
		mass = PARTICLE_MASS;
	}
	
	/** 
	 * Constructs particle with the specified material/undeformed
	 * coordinate, p0.
	 */
	TParticle(const TVec3<scalar_t>& x0) 
	{
		timeToLive = 1.0;
		x = x0;
		// name = -1;
		mass = PARTICLE_MASS;
	}
	
	/** 
	 * Constructs particle with the specified x, y, z initial position.
	 */
	TParticle(const scalar_t _x, const scalar_t _y, const scalar_t _z) 
	{
		timeToLive = 1.0;
		x.x = _x;
		x.y = _y;
		x.z = _z;
		// name = -1;
		mass = PARTICLE_MASS;
	}
	
	~TParticle() {}
	
	// Particle(const Particle& x0) 
	// {
	//		x = p.x;
	//		v = p.v;
	//		f = p.f;
	//		c = p.c;
	//		rho = p.rho;
	//		m = p.mass;
	// }
	
	void reset()
	{
		timeToLive = 1.0;
		x *= 0;
		v *= 0;
		f *= 0;
		mass = PARTICLE_MASS;
	}
	
	TParticle& operator=(const TParticle<scalar_t>& p)
	{
		x = p.x;
		v = p.v;
		f = p.f;
		c = p.c;
		rho = p.rho;
		mass = p.mass;
		return *this;
	}

	bool isAlive(double deltaT)
	{
		timeToLive -= deltaT;
		if(timeToLive > 0){
			return true;
		}
		else{
			return false;
		}
	}

	// void draw()
	// {
	// 	if(PARTICLE_DISPLAY_LIST < 0) {// MAKE DISPLAY LIST:
	// 		int displayListIndex = glGenLists(1);
	// 		glNewList(displayListIndex, GL_COMPILE);
	// 		drawParticle(0.0, 0.0, 0.0);	// for now we translate with a matrix
	// 		glEndList();
	// 		PARTICLE_DISPLAY_LIST = displayListIndex;
	// 	}
	// 
	// 	glColor3f(c.x, c.y, c.z);
	// 
	// 	/// DRAW ORIGIN-CIRCLE TRANSLATED TO "p":
	// 	glLoadName(name);
	// 	glPushMatrix();
	// 	glTranslated(x.x, x.y, x.z);
	// 	glCallList(PARTICLE_DISPLAY_LIST);
	// 	glPopMatrix();
	// }

private:
	/** age of the particle */
	scalar_t timeToLive;

	/** 
	 * Draws a canonical circular particle.
	 */
	// static void drawParticle(float z, float y, float z)
	// {
	// 	glPointSize(6);
	// 	glBegin(GL_POINTS);
	// 	glVertex3d(x, y, z);
	// 	glEnd();
	// }
};

typedef TParticle<double> Particle;
typedef TParticle<float> Particlef;

}	// namespace fdl

#endif	// __FDL_PARTICLE_H
