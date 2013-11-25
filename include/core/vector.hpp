#ifndef __FDL_VECTOR_H
#define __FDL_VECTOR_H

#include <iostream>
#include <sstream>
#include <math.h>

namespace fdl {

template<class T>
class TVec3 {
public:
	T x, y, z;

	inline TVec3(T _x = 0.0f, T _y = 0.0f, T _z = 0.0f)
		: x(_x), y(_y), z(_z) {
	}
	
	inline TVec3 operator+(const TVec3 &v) const {
		return TVec3<T>(x + v.x, y + v.y, z + v.z);
	}

	inline TVec3 operator-(const TVec3 &v) const {
		return TVec3<T>(x - v.x, y - v.y, z - v.z);
	}

	inline TVec3& operator+=(const TVec3 &v) {
		x += v.x; y += v.y; z += v.z;
		return *this;
	}

	inline TVec3& operator-=(const TVec3 &v) {
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}

	inline TVec3 operator*(T f) const {
		return TVec3<T>(x*f, y*f, z*f);
	}

	inline TVec3 &operator*=(T f) {
		x *= f;
		y *= f;
		z *= f;
		return *this;
	}

	inline TVec3 operator-() const {
		return TVec3<T>(-x, -y, -z);
	}

	inline TVec3 operator/(T f) const {
		T r = 1.0f / f;
		return TVec3<T>(x * r, y * r, z * r);
	}

	inline TVec3 &operator/=(T f) {
		T r = 1.0f / f;
		x *= r;
		y *= r;
		z *= r;
		return *this;
	}

	inline T operator[](int i) const {
		return (&x)[i];
	}

	inline T &operator[](int i) {
		return (&x)[i];
	}

	inline T lengthSquared() const {
		return x*x + y*y + z*z;
	}

	inline T length() const {
		return sqrt(lengthSquared());
	}

	inline TVec3 getSquareRoot() const {
		T sqrt_x = sqrt(x);
		T sqrt_y = sqrt(y);
		T sqrt_z = sqrt(z);
		return TVec3<T>(sqrt_x, sqrt_y, sqrt_z);
	}

	inline bool operator==(const TVec3 &v) const {
		return (v.x == x && v.y == y && v.z == z);
	}
	
	inline bool operator!=(const TVec3 &v) const {
		return !operator==(v);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << ", " << z << "]";
		return oss.str();
	}
};

template<class T>
class TVec2 {
public:
	T x, y;

	inline TVec2(T _x = 0.0f, T _y = 0.0f)
		: x(_x), y(_y) {
	}
	
	inline TVec2 operator+(const TVec2 &v) const {
		return TVec2<T>(x + v.x, y + v.y);
	}

	inline TVec2 operator-(const TVec2 &v) const {
		return TVec2<T>(x - v.x, y - v.y);
	}

	inline TVec2& operator+=(const TVec2 &v) {
		x += v.x; y += v.y;
		return *this;
	}

	inline TVec2& operator-=(const TVec2 &v) {
		x -= v.x; y -= v.y;
		return *this;
	}

	inline TVec2 operator*(T f) const {
		return TVec2<T>(x*f, y*f);
	}

	inline TVec2 &operator*=(T f) {
		x *= f;
		y *= f;
		return *this;
	}

	inline TVec2 operator-() const {
		return TVec2<T>(-x, -y);
	}

	inline TVec2 operator/(T f) const {
		T r = 1.0f / f;
		return TVec2<T>(x * r, y * r);
	}

	inline TVec2 &operator/=(T f) {
		T r = 1.0f / f;
		x *= r; y *= r;
		return *this;
	}

	inline T operator[](int i) const {
		return (&x)[i];
	}

	inline T &operator[](int i) {
		return (&x)[i];
	}

	inline T lengthSquared() const {
		return x*x + y*y;
	}

	inline T length() const {
		return sqrt(lengthSquared());
	}

	inline bool operator==(const TVec2 &v) const {
		return (v.x == x && v.y == y);
	}
	
	inline bool operator!=(const TVec2 &v) const {
		return !operator==(v);
	}

	inline std::string toString() const {
		std::ostringstream oss;
		oss << "[" << x << ", " << y << "]";
		return oss.str();
	}
};

template<class T>
inline T dot(const TVec3<T> &v1, const TVec3<T> &v2) {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

template<class T>
inline T dot(const TVec2<T> &v1, const TVec2<T> &v2) {
	return v1.x*v2.x + v1.y*v2.y;
}

template<class T>
inline TVec3<T> cross(const TVec3<T> &v1, const TVec3<T> &v2) {
	return TVec3<T>(
		(v1.y * v2.z) - (v1.z * v2.y), 
		(v1.z * v2.x) - (v1.x * v2.z),
		(v1.x * v2.y) - (v1.y * v2.x)
	);
}

typedef TVec3<double> Point3;
typedef TVec3<float> Point3f;
typedef TVec3<double> Vector3;
typedef TVec3<float> Vector3f;

typedef TVec2<double> Point2;
typedef TVec2<float> Point2f;
typedef TVec2<double> Vector2;
typedef TVec2<float> Vector2f;


}	// namespace fdl

#endif // __FDL_VECTOR_H

