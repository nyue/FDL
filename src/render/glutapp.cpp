#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <assert.h>
#include <math.h>

#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/signal.hpp>

#include "render/glutapp.h"

namespace fdl {
	
static GlutApp* self;

/**
 * Constructor allocates and initializes data
 * 
 */
GlutApp::GlutApp(std::string title, int width, int height) : m_windowTitle(title), m_width(width), m_height(height)
{
	die = 0;
	m_button1Down = false;
	m_button2Down = false;
	m_lastX = 0;
	m_lastY = 0;
	m_running = false;
	
	m_camera = new Camera();
	lastMousePt = new Imath::V2f();
	currMousePt = new Imath::V2f();
	
	self = this;
}

/**
 * Destructor de-allocates data
 * 
 */
GlutApp::~GlutApp()
{
	m_running = false;
}

/**
 * The init function begins the GLUT loopback
 * 
 */
int GlutApp::init()
{
	if(!m_running){
		//t1 = new boost::thread(boost::bind(&GlutApp::start,this));
		//t1->join();
		m_running = true;
		this->start();
	}
	
	return 0;
}

int GlutApp::render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//     glLoadIdentity();  
	// gluLookAt (0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);

	glColor3f(1.0, 0.0, 1.0);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	Imath::V3f* eye = new Imath::V3f();
	renderCamera(*eye);
	
	glutSolidTeapot(1.0);
	// glBegin(GL_QUADS);
	// 
	// glColor3f(0.5, 0.0, 1.0);
	// glVertex3f(-2, -2, 0);
	// glVertex3f(2, -2, 0);
	// glColor3f(1.0, 0.5, 0.0);
	// glVertex3f(2, 2, 0);
	// glVertex3f(-2, 2, 0);
	// 
	// glEnd();

	// swap drawing buffers
	glutSwapBuffers();
}

void GlutApp::windowToViewport(Imath::V2f& p) const
{
	p.setValue((float)((2.0f * p.x - m_width) / m_width),
	        (float)((2.0f * (m_height - p.y - 1.0) - m_height) / m_height));
}

void GlutApp::renderCamera(const Imath::V3f& eye)
{
	m_camera->updateMatrices();
	//eye.set(mainCamera.getEye());
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	float proj[16];
	float view[16];
	/*
	m_camera->getProjectionMatrix().toArray(proj);
	m_camera->getViewMatrix().toArray(view);
	glMultMatrixf(proj);
	glMultMatrixf(view);
	*/
    glMultMatrixf(m_camera->getProjectionMatrix().getValue());
    glMultMatrixf(m_camera->getViewMatrix().getValue());
//	rotationGizmo.setCamera(mainCamera);
}

int GlutApp::resize(int width, int height)
{
	m_width = width;
	m_height = height;

	glLoadIdentity();
	glViewport(0,0,width,height);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//gluPerspective (60, width / height, 1.0, 100.0);
	//glOrtho(0, self->m_width, self->m_height, 0, -1.0f, 1.0f);
	// glMatrixMode(GL_MODELVIEW);
	m_camera->setAspect(width / height);
}

int GlutApp::keyPressed(unsigned char key)
{
	// for escape key ...
	if (key == 27) {
		die = 1;
		glutDestroyWindow(window);
	}else{
		std::cout << "key pressed: " << key << std::endl;
	}
}

int GlutApp::keyReleased(unsigned char key)
{
}

int GlutApp::mousePressed(int button, int x, int y)
{
	// std::cout << "mouse pressed: " << x << ", " << y << std::endl;
	lastMousePt->setValue(x,y);
	windowToViewport(*lastMousePt);
	m_button1Down = true;
}

int GlutApp::mouseReleased(int button, int x, int y)
{
	// std::cout << "mouse released: " << x << ", " << y << std::endl;
	m_button1Down = false;
}

int GlutApp::mouseMoved(int x, int y)
{
	// std::cout << "mouse moved: " << x << ", " << y << std::endl;
	currMousePt->setValue(x,y);
	windowToViewport(*currMousePt);
	m_camera->orbit(*lastMousePt, *currMousePt);
	*lastMousePt = *currMousePt;
}

/**
 * The init function finds the device, and starts a thread that begins the input processing
 * 
 */
void GlutApp::render_cb()
{
	self->render();
}


/**
 * GLUT callback for window resizing
 *
 */
void GlutApp::resize_cb(int width, int height)
{
	self->resize(width,height);
}

/**
 * Process keyboard input callback
 *
 */
void GlutApp::keyPressed_cb(unsigned char key, int x, int y)
{	
	self->keyPressed(key);
}


/**
 * Mouse click callback function
 * 
 */
void GlutApp::mouse_cb(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN && button == GLUT_LEFT_BUTTON) {
		self->m_button1Down = true;
		self->m_lastX = x;
		self->m_lastY = y;
		self->mousePressed(button, x, y);
	}
	else if (state == GLUT_UP && button == GLUT_LEFT_BUTTON) {
		self->m_button1Down = false;
		self->mouseReleased(button, x, y);
	}
	else if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON) {
		self->m_button2Down = true;
		self->m_lastX = x;
		self->m_lastY = y;
		self->mousePressed(button, x, y);
	}
	else if (state == GLUT_UP && button == GLUT_RIGHT_BUTTON) {
		self->m_button2Down = false;
		self->mouseReleased(button, x, y);
	}
}


/**
 * Mouse movement callback function
 * 
 */
void GlutApp::motion_cb(int x, int y)
{
	if (self->m_button1Down) {
		double dx = (x - self->m_lastX)/2.0f;
		double dy = (y - self->m_lastY)/2.0f;
		self->m_lastX = x;
		self->m_lastY = y;
		glutPostRedisplay();
	}
	else if (self->m_button2Down) {
		// double dy = (y - lastY)/10.0f;
		// m_distZ -= dy;
		// m_lastY = y;
		// glutPostRedisplay();
	}
	
	self->mouseMoved(x, y);
}


/**
 * 
 * 
 */
void GlutApp::start()
{
	m_mutex.lock();
	
	int argc = 1;
	char* argv = (char*) m_windowTitle.c_str();
	
	glutInit(&argc, &argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
	glutInitWindowSize(self->m_width, self->m_height);
	glutInitWindowPosition(0, 0);

	// create X11 window
	window = glutCreateWindow(m_windowTitle.c_str());

	// pass callback handler function pointers to GLUT
	glutDisplayFunc(&render_cb);
	glutIdleFunc(&render_cb);
	glutReshapeFunc(&resize_cb);
	glutKeyboardFunc(&keyPressed_cb);
	glutMouseFunc(&mouse_cb);
	glutMotionFunc(&motion_cb);

	// init GL
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClearDepth(1.0);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glShadeModel(GL_SMOOTH);

	// begin loop
	glutMainLoop();
	
	m_running = false;

	m_mutex.unlock();
}

}
