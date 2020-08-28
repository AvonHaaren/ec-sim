//
//  PowerLaw.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 08/08/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef PowerLaw_h
#define PowerLaw_h

#include "potential.h"

class PowerLaw : public Potential
{
private:
    double x0,y0,z0;
    double Ux,Uy,Uz;
    double ax,ay,az;
    
public:
    PowerLaw(bool active, const dict<double>& pars) :
    Potential(active, {
        {"x0", x0},
        {"y0", y0},
        {"z0", z0},
        {"Ux", Ux},
        {"Uy", Uy},
        {"Uz", Uz},
        {"exp_x", ax},
        {"exp_y", ay},
        {"exp_z", az}
    }, pars) {}
    // Virtual constructor
    static std::unique_ptr<Potential> Create(bool active, const dict<double>& pars) {
        return std::make_unique<PowerLaw>(active, pars);
    }
    
    // Potential energy at position r
    double potential(const vector3d& r) const {
        return Ux*std::pow(std::abs(r.x/x0), ax) + Uy*std::pow(std::abs(r.y/y0), ay) + Uz*std::pow(std::abs(r.z/z0), az);
    }
    // Force on a particle at position r
    vector3d force(const vector3d& r) const {
        if (r == vector3d{0,0,0}) return vector3d{0,0,0};
        return vector3d(
                        -Ux*ax/x0/x0*r.x*std::pow(std::abs(r.x/x0), ax - 2),
                        -Uy*ay/y0/y0*r.y*std::pow(std::abs(r.y/y0), ay - 2),
                        -Uz*az/z0/z0*r.z*std::pow(std::abs(r.z/z0), az - 2)
                        );
    }
    // never eject from this type of potential
    bool withinBoundary(const vector3d& r) const { return true; }
};

#endif /* PowerLaw_h */
