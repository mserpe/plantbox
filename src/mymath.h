// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MYMATH_H_
#define MYMATH_H_

/**
 * My own minimalistic non-generic fixed dimension
 * not operator overloading vector matrix classes
 */

#include <cmath>
#include <sstream>
#include <assert.h>
#include <vector>
#include <functional>
#include <limits>


namespace CPlantBox {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/**
 * Vector2i stores two int values
 */
class Vector2i
{
public:

	Vector2i(): x(0), y(0) { } ///< Default constructor
	Vector2i(int x_, int y_): x(x_),y(y_) { } ///< Constructor passing two ints
	Vector2i(const Vector2i& v): x(v.x), y(v.y) { } ///< Copy constructor
	Vector2i(const std::vector<int>& xx): x(xx.at(0)), y(xx.at(1)) { } ///< for Python coupling

	std::string toString() const {
		std::ostringstream strs;
		strs << "( "<<x<<", "<<y<<" )";
		return strs.str();
	} ///< creates a string representing the two doubles

	int x; ///< number 1
	int y; ///< number 2

};



/**
 * Vector2d stores two double values
 */
class Vector2d
{
public:

	Vector2d() = default;
	Vector2d(double x_, double y_): x(x_),y(y_) { } ///< Constructor passing two doubles
	Vector2d(const Vector2d& v): x(v.x), y(v.y) { } ///< Copy Constructor
	Vector2d(const std::vector<double>& xx): x(xx.at(0)), y(xx.at(1)) { } ///< for Python coupling

	std::string toString() const {
		std::ostringstream strs;
		strs << "( "<<x<<", "<<y<<" )";
		return strs.str();
	} ///< creates a string representing the two doubles

	double x{0.0}; ///< number 1
	double y{0.0}; ///< number 2

};


/**
 * Vector3d stores three double values
 */
class Vector3d
{
public:

	Vector3d(): x(0),y(0),z(0) { } ///< Default constructor
	Vector3d(double x_, double y_, double z_): x(x_), y(y_), z(z_) { } ///< Constructor passing three doubles
	Vector3d(const Vector3d& v): x(v.x), y(v.y), z(v.z) { } ///< Copy Constructor
	Vector3d(const std::vector<double>& xx): x(xx.at(0)), y(xx.at(1)), z(xx.at(2)) { } ///< for Python coupling

	static Vector3d rotAB(double a, double b) { ///< first column of Rx(b)*Rz(a)
		double sa = sin(a);
		return Vector3d(cos(a), sa*cos(b), sa*sin(b) );
	};

	void normalize() { double l=length(); x/=l; y/=l; z/=l; } ///< normalizes the vector
  Vector3d normalized() const { double l=length(); return Vector3d(x/l,y/l,z/l); } ///< returns the normalized vector
  Vector3d scaledTo(double s) const { double l=length(); return Vector3d(x*s/l,y*s/l,z*s/l); } ///< returns the vector scaled to a given length

	double times(const Vector3d& v) const { return v.x*x+v.y*y+v.z*z; } ///< inner product (dot product, officially)
	double length() const { return sqrt(x*x+y*y+z*z); } ///< returns the Euclidian length

	Vector3d times(const double s) const { return Vector3d(s*x,s*y,s*z); } ///< returns the vector multiplied by a scalar value
	Vector3d plus(const Vector3d& v) const { return Vector3d(x+v.x,y+v.y,z+v.z); } ///< adds a vector and returns the result
	Vector3d minus(const Vector3d& v) const { return Vector3d(x-v.x,y-v.y,z-v.z); } ///< subtracts a vector and returns the result
	Vector3d cross(const Vector3d& v) const { return Vector3d(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); } ///< takes the cross product

  Vector3d operator+(const Vector3d& v) const { return plus(v); } ///< adds a vector and returns the result
  Vector3d operator-(const Vector3d& v) const { return minus(v); } ///< subtracts a vector and returns the result
  Vector3d operator*(const double s) const { return times(s); } ///< returns the vector multiplied by a scalar value
  Vector3d operator/(const double s) const { return times(1.0/s); } ///< returns the vector divided by a scalar value
  Vector3d operator-() const { return times(-1.0); } ///< returns the vector multiplied by -1
  friend Vector3d operator*(const double s, const Vector3d& v) { return v.times(s); } ///< commutative multiplication
  friend Vector3d operator/(const double s, const Vector3d& v) { return v.times(1.0/s); } ///< commutative division

	std::string toString() const {
		std::ostringstream strs;
		strs << "( "<<x<<", "<<y<<", "<<z<<" )";
		return strs.str();
	} ///< creates a string representing the three doubles

  

	double x; ///< double number 1
	double y; ///< double number 2
	double z; ///< double number 3

};

inline bool operator==(const Vector3d& lhs, const Vector3d& rhs){ return ((lhs.x==rhs.x) && (lhs.y==rhs.y) && (lhs.z==rhs.z)); } // needed for boost python indexing suite
inline bool operator!=(const Vector3d& lhs, const Vector3d& rhs){ return !(lhs == rhs); }

/**
 * 3x3 Matrix class, compatible with Vector3d for basic linear algebra
 * (i.e. exactly the operations needed for CPlantBox)
 */
class Matrix3d
{
public:

	Matrix3d(): r0(1,0,0), r1(0,1,0), r2(0,0,1) { } ///< Default constructor (identity matrix)
	Matrix3d(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32, double m33):
		r0(m11,m12,m13), r1(m21,m22,m23), r2(m31,m32,m33) { } ///< Constructor passing nine doubles as 3x3 matrix entries
	Matrix3d(const Vector3d& c1, const Vector3d& c2, const Vector3d& c3) :
		r0(c1.x,c2.x,c3.x), r1(c1.y,c2.y,c3.y), r2(c1.z,c2.z,c3.z) { } ///< Constructs the matrix from three column vectors
	Matrix3d(const Matrix3d& m): r0(m.r0), r1(m.r1), r2(m.r2) { } ///< Copy constructor

	static Matrix3d rotX(double a) {
		double ca = cos(a);
		double sa = sin(a);
		return Matrix3d(1,0,0,0,ca,-sa,0,sa,ca);
	} ///< Creates a rotation matrix around the X-axis
	static Matrix3d rotY(double a) {
		double ca = cos(a);
		double sa = sin(a);
		return Matrix3d(ca,0,sa,0,1,0,-sa,0,ca);
	} ///< Creates a rotation matrix around the Y-axis
	static Matrix3d rotZ(double a) {
		double ca = cos(a);
		double sa = sin(a);
		return Matrix3d(ca,-sa,0,sa,ca,0,0,0,1);
	} ///< Creates a rotation matrix around the Z-axis
	static Matrix3d rotAB(double a, double b) { ///< Rx(b)*Rz(a)
		auto rxb = rotX(b);
		rxb.times(rotZ(a));
		return rxb;
	};

	/**
	 * Creates an orthonormal system (ONS) around the vector v
	 *
	 * Remark: This is not unique, in Rootbox the ONS is rotated around a random angle to make it truly random.
	 * This is likely to be not optimal, but I guess it is fast enough anyway.
	 *
	 * @param v           vector, that will be normalized and then be the first column of the resulting ONS
	 *
	 * \return            three orthonormal column vectors
	 */
	static Matrix3d ons(Vector3d& v) {
		Vector3d v2;
		Vector3d v3;
		if ((std::abs(v.x)>=std::abs(v.y)) && (std::abs(v.x)>=std::abs(v.z))) { // choose x and z
			v2 = Vector3d(-v.z, 0, v.x);
			v3 = v.cross(v2);
		} else if ((std::abs(v.y)>=std::abs(v.x)) && (std::abs(v.y)>=std::abs(v.z))) { // choose y and z
			v2 = Vector3d(0,-v.z, v.y);
			v3 = v.cross(v2);
		} else if ((std::abs(v.z)>=std::abs(v.x)) && (std::abs(v.z)>=std::abs(v.y))) { // choose x and z
			v2 = Vector3d(-v.z, 0, v.x);
			v3 = v.cross(v2);
		}
		v.normalize();
		v2.normalize();
		v3.normalize();
		return Matrix3d(v,v2,v3);
	} ///< Creates an orthonormal system (ONS) around the vector v

	double det() const {
		return  r0.x*(r1.y*r2.z-r2.y*r1.z)-r0.y*(r1.x*r2.z-r1.z*r2.x)+r0.z*(r1.x*r2.y-r1.y*r2.x);
	} ///< determinant of the matrix

	Matrix3d inverse() const {
		double d = det();
		double idet = 1. / d;
		Matrix3d A;
		A.r0.x = (r1.y*r2.z-r2.y*r1.z)*idet;
		A.r0.y = (r0.z*r2.y-r0.y*r2.z)*idet;
		A.r0.z = (r0.y*r1.z-r0.z*r1.y)*idet;
		A.r1.x = (r1.z*r2.x-r1.x*r2.z)*idet;
		A.r1.y = (r0.x*r2.z-r0.z*r2.x)*idet;
		A.r1.z = (r1.x*r0.z-r0.x*r1.z)*idet;
		A.r2.x = (r1.x*r2.y-r2.x*r1.y)*idet;
		A.r2.y = (r2.x*r0.y-r0.x*r2.y)*idet;
		A.r2.z = (r0.x*r1.y-r1.x*r0.y)*idet;
		return A;
	} ///< calculates the inverse of the matrix

	Vector3d column(int i) const {
		assert((i>=0) && (i<3));
		switch(i) {
		case 0: return Vector3d(r0.x,r1.x,r2.x);
		case 1: return Vector3d(r0.y,r1.y,r2.y);
		case 2: return Vector3d(r0.z,r1.z,r2.z);
		}
		throw 0; // just to not produce a warning
	} ///< returns the i-th column of the matrix (i=0..2)

	Vector3d row(int i) const {
		assert((i>=0) && (i<3));
		switch(i) {
		case 0: return r0;
		case 1: return r1;
		case 2: return r2;
		}
		throw 0; // just to not produce a warning
	} ///< returns the i-th row of the matrix (i=0..2)

	void times(const Matrix3d& m) {
		r0 = Vector3d(r0.times(m.column(0)), r0.times(m.column(1)), r0.times(m.column(2)) );
		r1 = Vector3d(r1.times(m.column(0)), r1.times(m.column(1)), r1.times(m.column(2)) );
		r2 = Vector3d(r2.times(m.column(0)), r2.times(m.column(1)), r2.times(m.column(2)) );
	} ///< Multiplies matrix m from right

	Vector3d times(const Vector3d& v) const {
		return Vector3d(r0.times(v), r1.times(v), r2.times(v));
	} ///<  Multiplies vector v from right

	std::string toString() const {
		std::ostringstream strs;
		strs << r0.toString() << "\n" << r1.toString() << "\n" << r2.toString();
		return strs.str();
	} ///< creates a string representing the 3x3 matrix

	Vector3d r0; ///< row 1
	Vector3d r1; ///< row 2
	Vector3d r2; ///< row 3

};

/**
 * Quaternion class with common functions for CG
 * this is for general Quaternions, but we use it for H0 and H1
*/
class Quaternion {
  public:

    double w{};
    Vector3d v{};

    Quaternion() = default;
    ~Quaternion() = default;
    Quaternion(double w, Vector3d v) : w(w), v(v) {}
    Quaternion(const Quaternion& q) : w(q.w), v(q.v) {}
    Quaternion(Quaternion&& q) : w(q.w), v(q.v) {}
    Quaternion(std::initializer_list<double> l) {
      assert(l.size() == 4);
      auto it = l.begin();
      w = *it;
      ++it;
      v = Vector3d(*it, *(it+1), *(it+2));
    }

    Quaternion& operator=(const Quaternion& q) { w = q.w; v = q.v; return *this; }

    Quaternion& operator=(Quaternion&& q) { w = q.w; v = q.v; return *this; }

    Quaternion operator+(const Quaternion& q) const { return Quaternion(w + q.w, v + q.v); }

    Quaternion operator-(const Quaternion& q) const { return Quaternion(w - q.w, v - q.v); }

    Quaternion operator*(const Quaternion& q) const {
      return Quaternion(w * q.w - v.times(q.v), q.v.times(w) + v.times(q.w) + v.cross(q.v));
    }

    std::string toString() const {
      std::ostringstream strs;
      strs << w << " + " << v.toString();
      return strs.str();
    }

    friend Quaternion operator*(double s, const Quaternion& q) { return Quaternion(q.w * s, q.v * s); }

    Quaternion operator*(double s) const { return Quaternion(w * s, v * s); }

    Quaternion operator/(double s) const { return Quaternion(w / s, v / s); }

    Quaternion operator-() const { return Quaternion(-w, -v); }

    double& operator[](int i) {
      assert((i >= 0) && (i < 4));
      switch (i) {
        case 0: return w;
        case 1: return v.x;
        case 2: return v.y;
        case 3: return v.z;
      }
      throw 0; // just to not produce a warning
    }

    Quaternion& operator+=(const Quaternion& q) {
      w += q.w;
      v = v + q.v;
      return *this;
    }

    Quaternion& operator-=(const Quaternion& q) {
      w -= q.w;
      v = v - q.v;
      return *this;
    }

    Quaternion& operator*=(const Quaternion& q) {
      w = w * q.w - v.times(q.v);
      v = w * q.v + v * q.w + v.cross(q.v);
      return *this;
    }

    Quaternion& operator*=(double s) {
      w *= s;
      v = v * s;
      return *this;
    }

    Quaternion& operator/=(double s) {
      w /= s;
      v = v / s;
      return *this;
    }

    Quaternion conjugate() const {
      return Quaternion(w, -v);
    }

    Quaternion inverse() const {
      return conjugate() / norm2();
    }

    Quaternion inverse2() const {
      return conjugate() / (w * w + v.times(v));
    }

    void invert() {
      *this = inverse();
    }
    void conjugateInPlace() {
      v = -v;
    }

    inline double norm() const {
      return std::sqrt(w * w + v.times(v));
    }

    inline double norm2() const {
      return w * w + v.times(v);
    }

    inline Quaternion normalized() const {
      return *this / norm();
    }

    void normalize() {
      *this /= norm();
    }

    double dot(const Quaternion& q) const {
      return w * q.w + v.times(q.v);
    }

    // a method that calculates a quaternion from a forward vector.
    static Quaternion FromForward(const Vector3d& forward) {
      Vector3d up(0, 0, 1);
      Vector3d right = up.cross(forward);
      up = forward.cross(right);
      right.normalize();
      up.normalize();
      return Quaternion::FromMatrix3d(Matrix3d(right, up, forward));
    }

    /**
     * @brief Rotates a vector by the quaternion
     * @param v vector to rotate
     * @return rotated vector
     * @note the formula is v' = q * v * q^-1
    */
    Vector3d Rotate(const Vector3d& v) const {
      Quaternion standin = *this;
      standin.w = 1.0;
      return (standin * Quaternion(0, v) * standin.conjugate()).v;
    }

    /**
     * @brief A method that rotates a vector, ensuring uniform scale
     * @param v vector to rotate
     * @return rotated vector
     */
    Vector3d RotateUniform(const Vector3d& v) const {
      Quaternion standin = *this;
      standin.w = 1.0;
      return (standin * Quaternion(0, v) * standin.conjugate()).v * norm();
    }

    /**
     * @brief Converts the quaternion to a 3x3 rotation matrix
     * @return 3x3 rotation matrix
     * @note this is essentially the same as rotating the axis vectors by the quaternion
    */
    Matrix3d ToMatrix3d() const {
      return {{1 - 2 * v.y * v.y - 2 * v.z * v.z, 2 * v.x * v.y - 2 * w * v.z, 2 * v.x * v.z + 2 * w * v.y},
      {2 * v.x * v.y + 2 * w * v.z, 1 - 2 * v.x * v.x - 2 * v.z * v.z, 2 * v.y * v.z - 2 * w * v.x},
      {2 * v.x * v.z - 2 * w * v.y, 2 * v.y * v.z + 2 * w * v.x, 1 - 2 * v.x * v.x - 2 * v.y * v.y}};
    }

    static Quaternion FromMatrix3d(const Matrix3d& m)
    {
      double w = std::sqrt(1 + m.r0.x + m.r1.y + m.r2.z) / 2.0;
      double w4 = (4 * w);
      double x = (m.r0.y - m.r1.z) / w4;
      double y = (m.r0.z - m.r2.x) / w4;
      double z = (m.r1.x - m.r0.y) / w4;
      return Quaternion(w, Vector3d(x, y, z));
    }

    // a method that returns a look-at-direction representing the local x axis
    // TODO: PLEASE check with daniel about coordinate axis conventions!!
    Vector3d Forward() const 
    {
      return Rotate(Vector3d(1,0,0));
    }

    // a method that returns a look-at-direction representing the local y axis
    // TODO: SOON check with daniel about coordinate axis conventions!!
    Vector3d Right() const 
    {
      return Rotate(Vector3d(0,1,0));
    }

    // a method that returns a look-at-direction representing the local z axis
    // TODO: NOW check with daniel about coordinate axis conventions!!
    Vector3d Up() const 
    {
      return Rotate(Vector3d(0,0,1));
    }

    // a method that computes the normed linear sum of axis vectors
    // as I always use it with a set length, I am not normalizing the result
    Vector3d Axis(int x = 0, int y = 0, int z = 0) const
    {
      return (Forward() * x + Right() * y + Up() * z);
    }

    // a method that selects the shortest arc between two quaternions
    static Quaternion ShortestArc(const Quaternion& a, const Quaternion& b)
    {
      if (a.dot(b) < 0)
        return a * -1.0;
      else
        return a;
    }

    /**
     * @brief Static method to compute look-at-direction between two points
     * @param from the point to look from
     * @param to the point to look at
     * @param up the up direction
    */
    static Quaternion LookAt(const Vector3d& from, const Vector3d& to, const Vector3d& up)
    {
      // compute the forward direction
      auto forward = (to - from).normalized();
      // compute the right direction
      auto right = forward.cross(up).normalized();
      // compute the up direction
      auto up2 = right.cross(forward).normalized();
      // compute the rotation matrix
      Matrix3d m = {right, up2, -forward};
      // return the quaternion
      return FromMatrix3d(m);
    }

    // This Method calculates the Quaternion that rotates a to b
    // this is basically the same as the implementation in the Unreal Engine
    inline static Quaternion geodesicRotation(Vector3d a, Vector3d b)
    {
      double w = 1.0;
      // normalize the vectors
      a.normalize();
      b.normalize();
      // if the vectors are parallel, the rotation axis is undefined
      if(w < std::numeric_limits<double>::epsilon() && std::abs(a.x) <= std::abs(a.z))
        return Quaternion(w, Vector3d(0, -a.z, a.y));
      // case 2 of the above
      else if(w < std::numeric_limits<double>::epsilon() && std::abs(a.x) > std::abs(a.z))
        return Quaternion(w, Vector3d(-a.y, a.x, 0));
      // otherwise, the rotation axis is the cross product of a and b
      else
        return Quaternion(w, a.cross(b)).normalized();
    }


    // this method provides a spherical interpolation between two quaternions
    inline static Quaternion SphericalInterpolation(Quaternion a, Quaternion b, double t = 0.5)
    {
      double angle = a.dot(b);
      if (angle >= 1.0) ///< if angle is 1.0, the quaternions are the same
      {
        return a;
      }
      double half = acos(angle);
      double sinHalf = sin(half);
      if(std::abs(sinHalf) < std::numeric_limits<double>::epsilon()) //< if angle is 0.0, the quaternions are opposite
      {
        return {0.5*a.w + 0.5*b.w, 0.5*a.v + 0.5*b.v};
      }
      //< otherwise, interpolate
      double ratioA = sin((1.0 - t) * half) / sinHalf;
      double ratioB = sin(t * half) / sinHalf;
      return {ratioA * a.w + ratioB * b.w, ratioA * a.v + ratioB * b.v};
    }
};

/**
 * usefull
 */
class Function {
public:

	/**
	 * trapezoidal rule  (needed in dumux-rosi schroeder)
	 */
	static double quad(std::function<double(double)> f, double a, double b, int n) {
		double h = (b-a)/n;
		double s = h*(0.5*f(a) + 0.5*f(b));
		for (int i=1; i<n; i++) {
			s+=h*f(a+i*h);
		}
		return s;
	};

	/**
	 *  mimics numpy.linspace, from stack overflow
	 *  not in the binding, use numpy.linspace
	 */
	static std::vector<double> linspace(double start, double end, int num) {
		std::vector<double> linspaced;
		if (num == 0) {
			return linspaced;
		}
		if (num == 1) {
			linspaced.push_back(start);
			return linspaced;
		}
		double delta = (end - start) / (num - 1);
		for(int i=0; i < num-1; ++i) {
			linspaced.push_back(start + delta * i);
		}
		linspaced.push_back(end); // I want to ensure that start and end are exactly the same as the input
		return linspaced;
	}

	/**
	 * linearly interpolates a single value @param x
	 * not in the binding, use scipy instead
	 */
	static double interp1(double x, std::vector<double> x_,std::vector<double> y_) {
		if (x > x_.back()) { // check bounds
			return y_.back();
		}
		if (x < x_[0]) {
			return y_[0];
		}

		double lookUpIndex = std::distance(x_.begin(), std::lower_bound(x_.begin(), x_.end(), x));     // if we are within bounds find the index of the lower bound
		if (lookUpIndex == 0) {
			return y_[0];
		}
		double ipLinear = (x - x_[lookUpIndex-1])/(x_[lookUpIndex] - x_[lookUpIndex-1]);
		return (1.-ipLinear)*y_[lookUpIndex-1] + (ipLinear)*y_[lookUpIndex];
	}
  
  /**
   * Checks if a point is inside a cone frustum
   * (which is a cylinder with two different radii)
   * @param p point
   * @param x0 cone base center
   * @param r0 cone base radius
   * @param x1 cone top center
   * @param r1 cone top radius
  */
 bool isInsideFrustum(Vector3d p, Vector3d x0, double r0, Vector3d x1, double r1)
 {
   // compiler should be able to optimize this
    Vector3d v = x1-x0;
    Vector3d w = p-x0;
    double c1 = w.times(v);
    double c2 = v.times(v);
    double b = c1/c2;
    double r = r0 + b*(r1-r0);
    return (p-x0-b*v).length() < r;
 }
};

/**
 * Catmull-Rom spline interpolation
 * It is used to store the Catmull-Rom splines of depth 4
*/
class CatmullRomSpline
{
  public:
  CatmullRomSpline() = default;
  CatmullRomSpline(Vector3d y0, Vector3d y1, Vector3d y2, Vector3d y3, double t0, double t1) : y0(y0), y1(y1), y2(y2), y3(y3), t0(t0), t1(t1) {}
  CatmullRomSpline(std::vector<Vector3d> y, double t0, double t1) : y0(y[0]), y1(y[1]), y2(y[2]), y3(y[3]), t0(t0), t1(t1) {}
  Vector3d operator() (double t) const {
    double t_ = (t-t0)/(t1-t0);
    return 0.5 * ((2.0*y1) + (-y0 + y2) * t_ + (2.0*y0 - 5.0*y1 + 4.0*y2 - y3) * t_ * t_ + (-y0 + 3.0*y1 - 3.0*y2 + y3) * t_ * t_ * t_);
  }

  double getT0() const { return t0; }
  double getT1() const { return t1; }

  Vector3d derivative(double t) const {
    double t_ = (t-t0)/(t1-t0);
    return 0.5 * ((-y0 + y2) + (2.0*y0 - 5.0*y1 + 4.0*y2 - y3) * 2 * t_ + (-y0 + 3.0*y1 - 3.0*y2 + y3) * 3.0 * t_ * t_);
  }

  Vector3d operator[](int i) {
    switch(i) {
      case 0: return y0;
      case 1: return y1;
      case 2: return y2;
      case 3: return y3;
      default: throw std::runtime_error("CatmullRomSpline: index out of bounds");
    }
  }

  Quaternion computeOrientation(double t) const {
    Vector3d v = derivative(t);
    return Quaternion::FromForward(v);
  }

  private:
  // the control points of the spline
  Vector3d y0, y1, y2, y3;
  // start and end time of the spline
  double t0, t1;
};

/**
 * A manager class for the Catmull-Rom spline interpolation
 * It is used to store the Catmull-Rom splines of depth 4
 * We use a weighted average of the splines to get a smooth transition
*/
class CatmullRomSplineManager
{
  public:
  CatmullRomSplineManager() = default;
  CatmullRomSplineManager(std::vector<Vector3d> y) : y(y) {
    computeT();
  }

  CatmullRomSpline spline(int i) const {
    return splines[i];
  }
  
  // a method that selects the most suitable spline depending on most equal distance to start and finish
  CatmullRomSpline selectSpline(double t) const {
    double min = std::numeric_limits<double>::max();
    int index = 0;
    for(int i = 0; i < splines.size(); i++)
    {
      double d = std::abs(splines[i].getT0() - t) + std::abs(splines[i].getT1() - t);
      if(d < min)
      {
        min = d;
        index = i;
      }
    }
    return splines[index];
  }

  Vector3d operator() (double t) const {
    Vector3d p(0,0,0);
    int sum = 0;
    for(int i = 0; i < splines.size(); i++)
    {
      // whether the spline has t in its interval
      bool in = t >= splines[i].getT0() && t <= splines[i].getT1();
      // we add the spline if it is in the interval
      if(in)
      {
        p = p + splines[i](t);
        sum++;
      }
    }
    return p / static_cast<double>(sum);
  }
  void setY(std::vector<Vector3d> y) {
    this->y = y;
    computeT();
  }


  private:

  void computeT()
  {
    yt.clear();
    yt.push_back(0);
    for(int i = 1; i < y.size(); i++)
    {
      yt.push_back(yt[i-1] + (y[i]-y[i-1]).length());
    }
    for(int i = 0; i < yt.size(); i++)
    {
      yt[i] /= yt.back();
    }
    splines.clear();
    for(int i = 0; i < y.size()-3; i++)
    {
      splines.push_back(CatmullRomSpline({y[i], y[i+1], y[i+2], y[i+3]}, yt[i], yt[i+3]));
    }
  }

  std::vector<Vector3d> y; // control points
  std::vector<double> yt; // t values
  std::vector<CatmullRomSpline> splines;
};

} // end namespace CPlantBox


#endif
