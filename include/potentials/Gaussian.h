//
//  Gaussian.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 20/05/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef Gaussian_h
#define Gaussian_h


#include "potential.h"

class Gaussian : public Potential
{
private:
    // beam waist, rayleigh range and potential depth (at the origin)
    double w0, zR, U0;
    
    double square(const double& x) const { return x*x; }
    
public:
    // Constructor assigns friendly names to the parameters
    Gaussian(bool active, const dict<double>& pars) :
    Potential(
              active,
              {
        {"waist", w0},
        {"rayleigh_range", zR},
        {"depth", U0}
    },
              pars
              ) {}
    
    // Virtual constructor
    static std::unique_ptr<Potential> Create(bool active, const dict<double>& pars) {
        return std::make_unique<Gaussian>(active, pars);
    }
    
    // potential depth at position r
    double potential(const vector3d& r) const {
        return U0/(1+square(r.z/zR)) * exp(-2*(square(r.x)+square(r.y)) / (square(w0) * (1+square(r.z/zR))));
    }
    // force on a particle at position r
    vector3d force(const vector3d& r) const {
        double VarZ = (1 + square(r.z/zR));
        double Exp = exp(-2*(square(r.x)+square(r.y))/(square(w0)*VarZ));
        
        return vector3d(
                        4*r.x*Exp/(square(w0)*VarZ),
                        4*r.y*Exp/(square(w0)*VarZ),
                        2*r.z*Exp/(square(zR)*pow(VarZ, 3)) * (1 - 2 * (square(r.x)+square(r.y))/square(w0) + square(r.z))
                        ) * U0;
    }
    // the bounding box is a cylinder with a
    // radius of twice the beam waist and a
    // height of twice the rayleigh length
    bool withinBoundary(const vector3d& r) const {
        return (r.z < 2*zR && square(r.x) + square(r.y) < 4*square(w0));
    }
};
#endif /* Gaussian_h */
