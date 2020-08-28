//
//  CrossedBeam.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 19/08/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef CrossedBeam_h
#define CrossedBeam_h


#include "potential.h"
#include <tuple>


class CrossedBeam : public Potential
{
private:
    // ring radius, potential depth, exponent, angle between beams
    double R, U0, O, alpha;
    
    double square(const double& x) const { return x*x; }
    std::tuple<double, double> transform(const double& x, const double& y) const {
        double s = std::sin(alpha/2), c = std::cos(alpha/2);
        return {(-s*x + c*y)/R, (s*x + c*y)/R};
    }
    std::tuple<double, double> radiusXY(const vector3d& r) const {
        auto[xN,yN] = transform(r.x, r.y);
        return {sqrt(square(xN) + square(r.z/R)), sqrt(square(yN) + square(r.z/R))};
    }
    
public:
    // Constructor assigns friendly names to all the parameters
    CrossedBeam(bool active, const dict<double>& pars) :
    Potential(
              active,
              {
        {"radius", R},
        {"depth", U0},
        {"order", O},
        {"angle", alpha}
    },
              pars
              ) {}
    // Virtual constructor
    static std::unique_ptr<Potential> Create(bool active, const dict<double>& pars) {
        return std::make_unique<CrossedBeam>(active, pars);
    }
    
    // potential height at position r
    double potential(const vector3d& r) const {
        auto[rx, ry] = radiusXY(r);
        return U0 * (pow(rx, O) + pow(ry, O));
    }
    // force on a particle at position r
    vector3d force(const vector3d& r) const {
        auto [xN,yN] = transform(r.x,r.y);
        auto [rx,ry] = radiusXY(r);
        double D1 = pow(rx, O - 2), D2 = pow(ry, O - 2);
    
        return vector3d(std::sin(alpha/2) * (-xN*D1 + yN*D2),
                        std::cos(alpha/2) * (xN*D1 + yN*D2),
                        r.z/R * (D1 + D2)) * (-U0*O/R);
    }
    // the intersection of two cylinders creates this bounding box
    bool withinBoundary(const vector3d& r) const {
        auto[rx, ry] = radiusXY(r);
        return (rx < 1 && ry < 1);
    }
};



#endif /* CrossedBeam_h */
