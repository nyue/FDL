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

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include "render/camera.h"


namespace fdl {
	
static const double PI = 3.14159265;

    // using namespace cg::vecmath;

/**
 * The upper limit of the theta angle allowed when orbiting. Avoids gimbal lock.
 */
static const double THETA_LIMIT = 89.0 * PI / 180.0;

Camera::Camera()
{
	m_fov = 25.0;
	m_near = 0.1;
	m_far = 1000.0;
	m_aspect = 1;
	
	m_eye.setValue(5.0, 3.0, -3.0);
	m_target.setValue(0, 0, 0);
	m_up.setValue(0.0, 1.0, 0.0);
}

/**
 * Create a camera with a specific position
 * 
 * @param eye
 * @param target
 * @param up
 * @param camName
 */
Camera::Camera(	const Imath::V3f& eye, 
				const Imath::V3f& target, 
				const Imath::V3f& up, 
				float yFov, float near, float far)
{
	m_eye = eye;
	m_target = target;
	m_up = up;
	m_fov = yFov;
	m_near = near;
	m_far = far;
}

/**
 * Update the projection and view matrices
 */
void Camera::updateMatrices()
{
	/**
	 * Generate a perspective projection matrix for OpenGL
	 * (See: http://www.opengl.org/documentation/specs/man_pages/hardcopy/GL/html/glu/perspective.html)
	 * Avoid setting near too close to 0 or much Z-buffer precision will be lost!
	 */
	float thetaY = (float) PI / 180.0 * (m_fov / 2);
	float tanThetaY = (float) tan(thetaY);
	float h = 1.0 / tanThetaY;
	float ww = h / m_aspect;
	float qq = -(m_far + m_near) / (m_far - m_near);
	float qn = -2.0 * (m_far * m_near) / (m_far - m_near);

	m_projection = 0;
	// m_projection.m00 = ww;
	// m_projection.m11 = h;
	// m_projection.m22 = qq;
	// m_projection.m23 = qn;
	// m_projection.m32 = -1.0f;
	m_projection[0][0] = ww;
	m_projection[1][1] = h;
	m_projection[2][2] = qq;
	m_projection[2][3] = qn;
	m_projection[3][2] = -1.0;
	m_inverseProjection = m_projection.inverse();

	/**
	 * Generate a view transformation matrix for OpenGL
	 * (See: http://www.opengl.org/documentation/specs/man_pages/hardcopy/GL/html/glu/lookat.html)
	 */
	m_u = m_target - m_eye;
	m_u.normalize();
	m_v = m_up;
	m_v.normalize();
	m_w = m_u.cross(m_v);
	m_q = m_w.cross(m_u);
	m_w.normalize();
	m_q.normalize();

	m_orientation.makeIdentity();
	m_translate.makeIdentity();
	m_orientation[0][0] = m_w.x;	m_orientation[0][1] = m_w.y;	m_orientation[0][2] = m_w.z;
	m_orientation[1][0] = m_q.x;	m_orientation[1][1] = m_q.y;	m_orientation[1][2] = m_q.z;
	m_orientation[2][0] = -m_u.x;	m_orientation[2][1] = -m_u.y;	m_orientation[2][2] = -m_u.z;
	m_translate[0][3] = -m_eye.x;
	m_translate[1][3] = -m_eye.y;
	m_translate[2][3] = -m_eye.z;
	m_inverseOrientation = m_orientation.transpose();
	m_view = m_orientation * m_translate;
	m_inverseView = m_view.inverse();
}

/**
 * Returns the eye of the camera.
 * 
 * @return
 */
Imath::V3f Camera::getEye() const {
	return m_eye;
}

/**
 * Sets the eye of the camera.
 */
void Camera::setEye(const Imath::V3f& eyeVec) {
	m_eye = eyeVec;
}

/**
 * Sets the target of the camera.
 */
void Camera::setTarget(const Imath::V3f& targetVec) {
	m_target = targetVec;
}


/**
 * Moves the camera forwards and backwards relative to the target.
 * 
 * @param last
 * @param cur
 */
void Camera::dolly(float amount)
{
	Imath::V3f diff;
	diff = m_eye - m_target;
	m_eye = (diff * amount) + m_eye;
}

/**
 * Orbits the camera around the target.
 * 
 * @param lastMousePoint
 * @param currMousePoint
 */
void Camera::orbit(const Imath::V2f& lastMousePoint, const Imath::V2f& currMousePoint)
{
	Imath::V2f mouseDelta(currMousePoint);
	mouseDelta -= lastMousePoint;

	// Construct an arbitrary frame at the target with the z-axis the up
	// vector
	m_w = m_up;
	m_w.normalize();
	m_u = nonParallelVector(m_w);
	m_v = m_w.cross(m_u);
	m_v.normalize();
	m_u = m_v.cross(m_w);
	/*
	m_basis.setColumn(0, m_u);
	m_basis.setColumn(1, m_v);
	m_basis.setColumn(2, m_w);
	*/
	m_basis = Imath::M33f(m_u.x,m_v.x,m_w.x,
	        m_u.y,m_v.y,m_w.y,
	        m_u.z,m_v.z,m_w.z);
	Imath::M44f frame(m_basis, m_target);
	Imath::M44f frameInv;
	frameInv = frame.inverse();

	// write eye in that frame
	Imath::V4f e(m_eye.x, m_eye.y, m_eye.z, 0.0);
	// e = frameInv * e;
    e = e * frameInv;

	// write e in spherical coordinates
	float r = e.length();
	float phi = atan2(e.y, e.x);
	float theta = asin(e.z / r);

	// increment phi and theta by mouse motion
	phi += -PI / 2.0 * mouseDelta.x;
	theta += -PI / 2.0 * mouseDelta.y;
	if (theta > THETA_LIMIT) {
		theta = THETA_LIMIT;
	}
	if (theta < -THETA_LIMIT) {
		theta = -THETA_LIMIT;
	}

	// write e back in cartesian world coords
	e = Imath::V4f(	(float)(r * cos(phi) * cos(theta)),
			(float)(r * sin(phi) * cos(theta)),
			(float)(r * sin(theta)), 1.0);

	m_eye.normalize();
	// e = frame * e;
    e = e * frame;
	m_eye.setValue(e.x,e.y,e.z);
}

void Camera::setAspect(float d)
{
	m_aspect = d;
}

void Camera::setDirection(const Imath::V3f& direction)
{
	m_target = m_eye;
	m_target += direction;
}

void Camera::setUp(const Imath::V3f& up)
{
	m_up = up;
}

Imath::V3f Camera::nonParallelVector(const Imath::V3f& v)
{
	int i = argmin(abs(v.x), abs(v.y), abs(v.z));
	Imath::V3f u;
	if (i == 0) {
		u.x = 1.0;
	} else if (i == 1) {
		u.y = 1.0;
	} else if (i == 2) {
		u.z = 1.0;
	}
	return u;
}

}
