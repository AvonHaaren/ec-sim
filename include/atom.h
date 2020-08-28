#ifndef ATOM_H
#define ATOM_H

#include "vector3d.h"       // own simple 3-vector struct with operators         

/// Simple struct that contains the position, velocity and current movement timestep for an atom
struct atom {
    /// position vector
    vector3d r;
    /// velocity vector
    vector3d v;
    /// movement timestep
    double dtMove;
    /// mass
    double mass;
    /// unique identifier
    uint64_t ID;
    /// identifier of the atom last collided with
    uint64_t lastCollided;
    /// Cell index
    unsigned int cellIndex;
    
    /// Constructor for an atom object
    /// @param r0 a 3-vector with the initial position
    /// @param v0 a 3-vector with the initial velocity
    /// @param dt0 a double with the initial timestep
    /// @param m a double with the mass
    /// @param N an integer identifier (unique)
    /// @param lastCollided_ID the unique identifier of the atom with which the last collision took place
    atom(vector3d r0, vector3d v0, double dt0, double m, uint64_t N, uint64_t lastCollided_ID = 0) noexcept :
    r(r0),
    v(v0),
    dtMove(dt0),
    mass(m),
    ID(N),
    lastCollided(lastCollided_ID) {};
    /** Copy Constructor */
    atom(const atom &obj) noexcept :
    r(obj.r),
    v(obj.v),
    dtMove(obj.dtMove),
    mass(obj.mass),
    ID(obj.ID),
    lastCollided(obj.lastCollided),
    cellIndex(obj.cellIndex)
    {};

    // Equality and inequality operators
    bool operator == (const atom& obj) const { return ID == obj.ID; }
    bool operator != (const atom& obj) const { return ID != obj.ID; }
};


#endif // atom_h
