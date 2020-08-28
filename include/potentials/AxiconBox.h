//
//  AxiconBox.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 20/05/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef AxiconBox_h
#define AxiconBox_h

#include "potential.h"
#include <tuple>


class AxiconBox : public Potential
{
private:
    // ring radius, wall thickness (1/e^2), potential depth, plateau order
    double R, W, U0, O;

    double square(const double& x) const { return x*x; }
    std::tuple<double, double> radiusXY(const vector3d& r) const {
        return {sqrt(square(r.y) + square(r.z)), sqrt(square(r.x) + square(r.z))};
    }
    int sgn(double x) const { return (x > 0) - (x < 0); }
                     
public:
    // Constructor assigns friendly names to all the parameters
    AxiconBox(bool active, const dict<double>& pars) :
    Potential(
              active,
    {
        {"radius", R},
        {"wall_thickness", W},
        {"depth", U0},
        {"order", O}
    },
              pars
    ) {}
    // Virtual constructor
    static std::unique_ptr<Potential> Create(bool active, const dict<double>& pars) {
        return std::make_unique<AxiconBox>(active, pars);
    }
    
    // potential height at position r
    double potential(const vector3d& r) const {
        auto[rx, ry] = radiusXY(r);
        double powX = -2 * pow(abs(rx - R) / W, O), powY = -2 * pow(abs(ry - R) / W, O);
        return U0 * ((powX < -60 ? 0 : exp(powX)) + (powY < -60 ? 0 : exp(powY)));
    }
    // force on a particle at position r
    vector3d force(const vector3d& r) const {
        auto[rx, ry] = radiusXY(r);
        double beamx, beamy;
        double helperX = abs(rx-R)/W, helperY = abs(ry-R)/W;
        double powX = -2. * pow(helperX, O - 1), powY = -2. * pow(helperY, O - 1);
        beamx = O/W * powX * sgn(rx-R) * (powX * helperX < -60 ? 0 : exp(powX*helperX)) / rx;
        beamy = O/W * powY * sgn(ry-R) * (powY * helperY < -60 ? 0 : exp(powY*helperY)) / ry;

        return vector3d(r.x*beamy, r.y*beamx, r.z*(beamx + beamy)) * (-1.*U0);
    }
    // the intersection of two cylinders creates this bounding box
    bool withinBoundary(const vector3d& r) const {
        auto[rx, ry] = radiusXY(r);
        return (rx < R && ry < R);
    }
};


#endif /* AxiconBox_h */
