//
//  Harmonic.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 20/05/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef Harmonic_h
#define Harmonic_h

#include "potential.h"

class Harmonic : public Potential
{
private:
    // trap frequencies, mass, potential depth
    double omega_x, omega_y, omega_z, m, U0;
    
    double square(const double& x) const { return x*x; }
    
public:
    // Constructor assigns friendly names to the parameters
    Harmonic(bool active, const dict<double>& pars) :
    Potential(
              active,
              {
        {"omega_x", omega_x},
        {"omega_y", omega_y},
        {"omega_z", omega_z},
        {"mass", m},
        {"depth", U0}
    },
              pars
              ) {}
    // Virtual constructor
    static std::unique_ptr<Potential> Create(bool active, const dict<double>& pars) {
        return std::make_unique<Harmonic>(active, pars);
    }
    
    // potential depth at position r
    double potential(const vector3d& r) const {
        return 0.5*m*(square(r.x * omega_x) +
                      square(r.y * omega_y) +
                      square(r.z * omega_z));
    }
    // force on a particle at position r
    vector3d force(const vector3d& r) const {
        return vector3d(r.x * square(omega_x),
                        r.y * square(omega_y),
                        r.z * square(omega_z)) * (-1*m);
    }
    // the boundary is set by the cutoff depth
    bool withinBoundary(const vector3d& r) const {
        return potential(r) < U0;
    }
};


#endif /* Harmonic_h */
