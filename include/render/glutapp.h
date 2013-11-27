#ifndef __GLUTAPP_H
#define __GLUTAPP_H
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <assert.h>
#include <math.h>

#include <boost/bind.hpp>
#include <boost/mem_fn.hpp>
#include <boost/function.hpp>
#include <boost/thread.hpp>
#include <boost/signal.hpp>

#include "cg/vecmath/vec2.hpp"
#include "cg/vecmath/vec3.hpp"

#include "render/camera.h"

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#ifndef _WIN32
#include <sys/time.h>
#else
#include <ctime>
#include <windows.h>

int gettimeofday (struct timeval *tv, void* tz)
{
	union
	{
		long long ns100;
		FILETIME ft;
	} now;

	GetSystemTimeAsFileTime (&now.ft);
	tv->tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
	tv->tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL);
	return (0);
}
#endif

namespace fdl {

static int window;
static int die;

class GlutApp {
public:
	GlutApp(std::string title="Glut Window", int width=800, int height=600);
	virtual ~GlutApp();
	int init();
	int getWidth() { return m_width; }
	int getHeight() { return m_height; }
	void windowToViewport(cg::vecmath::Point2f& p) const;
	void renderCamera(const cg::vecmath::Vector3f& eye);
	
protected:
	void start();
	int m_width;
	int m_height;
	bool m_button1Down;
	bool m_button2Down;
	int m_lastX;
	int m_lastY;
	std::string m_windowTitle;
	boost::mutex m_mutex;
	Camera* m_camera;
	cg::vecmath::Point2f* lastMousePt;
	cg::vecmath::Point2f* currMousePt;
	
	virtual int render();
	virtual int resize(int width, int height);
	virtual int keyPressed(unsigned char key);
	virtual int keyReleased(unsigned char key);
	virtual int mousePressed(int button, int x, int y);
	virtual int mouseReleased(int button, int x, int y);
	virtual int mouseMoved(int x, int y);

private:	
	static void render_cb();
	static void resize_cb(int width, int height);
	static void keyPressed_cb(unsigned char key, int x, int y);
	static void mouse_cb(int button, int state, int x, int y);
	static void motion_cb(int x, int y);
	bool m_running;

};	// GlutApp

}	// namespace fdl

#endif	// __GLUTAPP_H
