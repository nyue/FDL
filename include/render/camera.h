#ifndef __CAMERA_H
#define __CAMERA_H

#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <assert.h>
#include <math.h>

//#include "cg/vecmath/vec2.hpp"
//#include "cg/vecmath/vec3.hpp"
//#include "cg/vecmath/vec4.hpp"
//#include "cg/vecmath/mat4.hpp"
//#include "cg/vecmath/mat3.hpp"
#include <OpenEXR/ImathMatrix.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

namespace fdl {

class Camera {
public:
	/**
	 * The default camera position
	 */
	Camera();

	/**
	 * Create a camera with a specific position
	 * 
	 * @param eye
	 * @param target
	 * @param up
	 * @param camName
	 */
	Camera(	const Imath::V3f& eye, const Imath::V3f& target,
			const Imath::V3f& up, float yFov, float near, float far);
	
	~Camera();
	
	/**
	 * Update the projection and view matrices
	 */
	void updateMatrices();

	/**
	 * Return the orientation matrix
	 */
	Imath::M44f getOrientationMatrix() const {
		return m_orientation;
	}

	/**
	 * Return the inverse of the orientation matrix
	 */
	Imath::M44f getInverseOrientationMatrix() const {
		return m_inverseOrientation;
	}
	
	/**
	 * Return the view matrix
	 */
	Imath::M44f getViewMatrix() const {
		return m_view;
	}

	/**
	 * Return the inverse of the view matrix
	 */
	Imath::M44f getInverseViewMatrix() const {
		return m_inverseView;
	}

	/**
	 * Return the projection matrix
	 */
	Imath::M44f getProjectionMatrix() const {
		return m_projection;
	}

	/**
	 * Return the inverse of the projection matrix
	 */
	Imath::M44f getInverseProjectionMatrix() const {
		return m_inverseProjection;
	}

	/**
	 * Returns the eye of the camera.
	 * 
	 * @return
	 */
	Imath::V3f getEye() const;

	/**
	 * Sets the eye of the camera.
	 */
	void setEye(const Imath::V3f& eyeVec);

	/**
	 * Sets the target of the camera.
	 */
	void setTarget(const Imath::V3f& targetVec);
	

	/**
	 * Moves the camera forwards and backwards relative to the target.
	 * 
	 * @param last
	 * @param cur
	 */
	void dolly(float amount);

	/**
	 * Orbits the camera around the target.
	 * 
	 * @param lastMousePoint
	 * @param currMousePoint
	 */
	void orbit(const Imath::V2f& lastMousePoint, const Imath::V2f& currMousePoint);

	void setAspect(float d);
	
	void setDirection(const Imath::V3f& direction);
	
	void setUp(const Imath::V3f& up);

protected:
	/** The YFOV in degrees */
	float m_fov;

	/** The near plane distance */
	float m_near;

	/** The far plane distance */
	float m_far;

	/** The aspect ratio */
	float m_aspect;

	// The position and orientation of the camera
	Imath::V3f m_eye;
	Imath::V3f m_target;
	Imath::V3f m_up;

	// Some temporary space for computing camera motions and transformations
	Imath::V3f m_u;
	Imath::V3f m_v;
	Imath::V3f m_w;
	Imath::V3f m_q;
	Imath::M33f m_basis;
	
private:
	/** Projection matrix */
	Imath::M44f m_projection;
	Imath::M44f m_inverseProjection;

	/** Orientation matrix (view w/o translate) */
	Imath::M44f m_orientation;
	Imath::M44f m_inverseOrientation;

	/** View matrix */
	Imath::M44f m_view;
	Imath::M44f m_inverseView;
	Imath::M44f m_translate;
	
	/**
	 * Returns 0 if a is smallest, 1 if b is smallest, 2 if c is smallest.
	 */
	static int argmin(double a, double b, double c) {
		return a < b ? (a < c ? 0 : 2) : (b < c ? 1 : 2);
	}

	/**
	 * Returns a vector that is not nearly parallel to v.
	 */
	static Imath::V3f nonParallelVector(const Imath::V3f& v);
	
};

}

#endif	// __CAMERA_H
