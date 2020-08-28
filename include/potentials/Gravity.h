//
//  Gravity.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 20/05/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef Gravity_h
#define Gravity_h

#include "potential.h"

class Gravity : public Potential
{
private:
    double m;                               // mass
    static constexpr double g = 9.81;       // acceleration (m/s^2)
    
public:
    // Constructor gives m the friendly name "mass"
    Gravity(bool active, const dict<double>& pars) :
    Potential(active, {{"mass", m}}, pars) {}
    // Virtual constructor
    static std::unique_ptr<Potential> Create(bool active, const dict<double>& pars) {
        return std::make_unique<Gravity>(active, pars);
    }
    
    // Potential energy at position r
    double potential(const vector3d& r) const { return m*g*r.z; }
    // Force on a particle at position r
    vector3d force(const vector3d& r) const { return vector3d(0,0,-1)*m*g; }
    // Gravity does not have a boundary
    bool withinBoundary(const vector3d& r) const { return false; }
};


#endif /* Gravity_h */
