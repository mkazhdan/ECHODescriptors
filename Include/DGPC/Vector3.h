/************************************************************
 * This file is part of the DGPC library. The library computes
 * Discrete Geodesic Polar Coordinates on a polygonal mesh.
 *
 * More info:
 *   http://folk.uio.no/eivindlm/dgpc/
 *
 * Authors: Eivind Lyche Melvær and Martin Reimers
 * Centre of Mathematics and Department of Informatics
 * University of Oslo, Norway, 2012
 ************************************************************/
#ifndef VECTOR3_H
#define VECTOR3_H

namespace DGPC {
  template<class real_t>
  class Vector3 {
    real_t p_[3];
    void set(real_t x, real_t y, real_t z) {
      p_[0] = x;
      p_[1] = y;
      p_[2] = z;
    }
  public:
    typedef real_t value_type;
    typedef Vector3<real_t> vector_type;

    static const size_t size_ = 3;

    Vector3() {
      p_[0] = p_[1] = p_[2] = std::numeric_limits<real_t>::max();
    }

    Vector3(real_t x, real_t y, real_t z) {
      set(x, y, z);
    }

    Vector3(const real_t v[]) {
      set(v[0], v[1], v[2]);
    }

    Vector3(const Vector3<real_t>& p) {
      set(p.x(), p.y(), p.z());
    }

    const real_t& x() const { return p_[0]; };
          real_t& x()       { return p_[0]; };
    const real_t& y() const { return p_[1]; };
          real_t& y()       { return p_[1]; };
    const real_t& z() const { return p_[2]; };
          real_t& z()       { return p_[2]; };

    const real_t& operator [] (int i) const { return p_[i]; };
    real_t& operator [] (int i)             { return p_[i]; };

    Vector3<real_t>& operator= (const Vector3<real_t>& p) {
      set(p.x(), p.y(), p.z());
      return *this;
    }

    Vector3<real_t> operator- (const Vector3<real_t>& v) const {
      return Vector3<real_t>(x() - v.x(), y() - v.y(), z() - v.z());
    }

    Vector3<real_t> operator+ (const Vector3<real_t>& v) const {
      return Vector3<real_t>(x() + v.x(), y() + v.y(), z() + v.z());
    }

    real_t operator* (const Vector3<real_t>& v) const {
      return x()*v.x() + y()*v.y() + z()*v.z();
    }

    Vector3<real_t> operator* (real_t d) const {
      return Vector3<real_t>(x()*d, y()*d, z()*d);
    }

    real_t dist(const Vector3<real_t>& v) const {
      return sqrt(dist2(v));
    }

    real_t dist2(const Vector3<real_t>& v) const {
      real_t dx = x() - v.x();
      real_t dy = y() - v.y();
      real_t dz = z() - v.z();
      return dx*dx + dy*dy + dz*dz;
    }

    real_t length() const {
      return sqrt(length2());
    }

    real_t length2() const {
      return x()*x() + y()*y() + z()*z();
    }

    Vector3<real_t> crossProd(const Vector3<real_t>& v) const {
      return Vector3(y()*v.z() - z()*v.y(), 
                    z()*v.x() - x()*v.z(), 
                    x()*v.y() - y()*v.x());
    }

    Vector3<real_t>& normalize() {
      const real_t len2 = length2();
      if (len2) {
        const real_t len = sqrt(len2);
        p_[0] /= len;
        p_[1] /= len;
        p_[2] /= len;
      }
      return *this;
    }
  };

  template<class real_t>
    Vector3<real_t> operator* (real_t d, const Vector3<real_t>& v)
    {
      return v*d;
    }

} //End namespace DGPC

#endif //VECTOR3_H
