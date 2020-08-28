#ifndef VECTOR3D_H
#define VECTOR3D_H

/// Simple 3-vector struct with mathematical operators
struct vector3d {
    /// x value
    double x;
    /// y value
    double y;
    /// z value
    double z;

    /// Constructor for a vector3d object
    /// @param x0 the initial x value
    /// @param y0 the initial y value
    /// @param z0 the initial z value
    vector3d(double x0 = 0, double y0 = 0, double z0 = 0) noexcept : x(x0), y(y0), z(z0) {}
    /// Copy constructor
    vector3d(const vector3d &obj) noexcept : x(obj.x), y(obj.y), z(obj.z) {}
    
// Addition (and in-place)
    vector3d operator + (const vector3d &obj) const { return {x + obj.x, y + obj.y, z + obj.z}; }
    void operator += (const vector3d &obj) { x += obj.x; y += obj.y; z += obj.z; }

// Subtraction (and in-place)
    vector3d operator - (const vector3d &obj) const { return {x - obj.x, y - obj.y, z - obj.z}; }
    void operator -= (const vector3d &obj) { x -= obj.x; y -= obj.y; z -= obj.z; }

// Componentwise Multiplication (and in-place)
    vector3d operator * (const double &val) const { return {x * val, y * val, z * val}; }
    void operator *= (const double &val) { x *= val; y *= val; z *= val; }

// Componentwise Division (and in-place)
    vector3d operator / (const double &val) const { return {x / val, y / val, z / val}; }
    void operator /= (const double &val) { x /= val; y /= val; z /= val; }

// Equality
    bool operator == (const vector3d &obj) const { return x == obj.x && y == obj.y && z == obj.z; }
    bool operator != (const vector3d &obj) const { return not (*this == obj); }
    
// inner product
    double operator * (const vector3d &obj) const { return x*obj.x + y*obj.y + z*obj.z; }
};

#endif
