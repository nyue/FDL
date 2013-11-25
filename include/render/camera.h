#ifndef __CAMERA_H
#define __CAMERA_H

#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <assert.h>
#include <math.h>

#include "cg/vecmath/vec2.hpp"
#include "cg/vecmath/vec3.hpp"
#include "cg/vecmath/vec4.hpp"
#include "cg/vecmath/mat4.hpp"
#include "cg/vecmath/mat3.hpp"

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
	Camera(	const cg::vecmath::Vector3f& eye, const cg::vecmath::Vector3f& target, 
			const cg::vecmath::Vector3f& up, float yFov, float near, float far);
	
	~Camera();
	
	/**
	 * Update the projection and view matrices
	 */
	void updateMatrices();

	/**
	 * Return the orientation matrix
	 */
	cg::vecmath::Matrix4f getOrientationMatrix() const {
		return m_orientation;
	}

	/**
	 * Return the inverse of the orientation matrix
	 */
	cg::vecmath::Matrix4f getInverseOrientationMatrix() const {
		return m_inverseOrientation;
	}
	
	/**
	 * Return the view matrix
	 */
	cg::vecmath::Matrix4f getViewMatrix() const {
		return m_view;
	}

	/**
	 * Return the inverse of the view matrix
	 */
	cg::vecmath::Matrix4f getInverseViewMatrix() const {
		return m_inverseView;
	}

	/**
	 * Return the projection matrix
	 */
	cg::vecmath::Matrix4f getProjectionMatrix() const {
		return m_projection;
	}

	/**
	 * Return the inverse of the projection matrix
	 */
	cg::vecmath::Matrix4f getInverseProjectionMatrix() const {
		return m_inverseProjection;
	}

	/**
	 * Returns the eye of the camera.
	 * 
	 * @return
	 */
	cg::vecmath::Vector3f getEye() const;

	/**
	 * Sets the eye of the camera.
	 */
	void setEye(const cg::vecmath::Vector3f& eyeVec);

	/**
	 * Sets the target of the camera.
	 */
	void setTarget(const cg::vecmath::Vector3f& targetVec);
	

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
	void orbit(const cg::vecmath::Point2f& lastMousePoint, const cg::vecmath::Point2f& currMousePoint);

	void setAspect(float d);
	
	void setDirection(const cg::vecmath::Vector3f& direction);
	
	void setUp(const cg::vecmath::Vector3f& up);

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
	cg::vecmath::Vector3f m_eye;
	cg::vecmath::Vector3f m_target;
	cg::vecmath::Vector3f m_up;

	// Some temporary space for computing camera motions and transformations
	cg::vecmath::Vector3f m_u;
	cg::vecmath::Vector3f m_v;
	cg::vecmath::Vector3f m_w;
	cg::vecmath::Vector3f m_q;
	cg::vecmath::Matrix3f m_basis;
	
private:
	/** Projection matrix */
	cg::vecmath::Matrix4f m_projection;
	cg::vecmath::Matrix4f m_inverseProjection;

	/** Orientation matrix (view w/o translate) */
	cg::vecmath::Matrix4f m_orientation;
	cg::vecmath::Matrix4f m_inverseOrientation;

	/** View matrix */
	cg::vecmath::Matrix4f m_view;
	cg::vecmath::Matrix4f m_inverseView;
	cg::vecmath::Matrix4f m_translate;
	
	/**
	 * Returns 0 if a is smallest, 1 if b is smallest, 2 if c is smallest.
	 */
	static int argmin(double a, double b, double c) {
		return a < b ? (a < c ? 0 : 2) : (b < c ? 1 : 2);
	}

	/**
	 * Returns a vector that is not nearly parallel to v.
	 */
	static cg::vecmath::Vector3f nonParallelVector(const cg::vecmath::Vector3f& v);
	
};

}

#endif	// __CAMERA_H