//
//  rungeKutta.cpp
//  EvaporationSim
//
//  Created by Andreas von Haaren on 27/02/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#include <algorithm>            // needed for Visual Studio
#include <cmath>                // needed for xCode

#include "rungeKutta.h"

void rungeKuttaStep(atom& a, const std::vector<std::unique_ptr<Potential>>& potentials, double timeStep)
{
    static constexpr double tolerance = 1.E-7;              // maximum relative change in energy permitted
    double timePassed = 0.;                                 // time after start of this step (in seconds)
    int timeStep_in_ns = (int)(round(timeStep * 1e9));      // timeStep in nanoseconds

    // Optimal initial timestep
    double stepSize = std::min(a.dtMove, timeStep);
    // empty if statement (can be used for debugging)
    // will be optimised away by compiler in release mode
    if (stepSize < 1.E-9) {
        
    }
    
    vector3d K[4] = {vector3d(), vector3d(), vector3d(), vector3d()};   // intermediate results 1
    vector3d L[4] = {vector3d(), vector3d(), vector3d(), vector3d()};   // intermediate results 2
    
    vector3d newR, newV;    // new position, velocity
    vector3d oldR = a.r, oldV = a.v;
    double E0, E, dE;       // old energy, new energy, relative change
    
    // keep going until timeStep is completed
    while ((int)(round(timePassed*1E9)) < timeStep_in_ns) {
        K[0] = {0,0,0};
        for (const auto& p : potentials)
            if (p->isActive())
                K[0] += p->force(a.r);
        K[0] *= stepSize / a.mass;
        L[0] = a.v*stepSize;
        
        K[1] = {0,0,0};
        for (const auto& p : potentials)
            if (p->isActive())
                K[1] += p->force(a.r + L[0]*0.5);
        K[1] *= stepSize / a.mass;
        L[1] = (a.v + K[0]*0.5)*stepSize;
        
        K[2] = {0,0,0};
        for (const auto& p : potentials)
            if (p->isActive())
                K[2] += p->force(a.r + L[1]*0.5);
        K[2] *= stepSize / a.mass;
        L[2] = (a.v + K[1]*0.5)*stepSize;
        
        K[3] = {0,0,0};
        for (const auto& p : potentials)
            if (p->isActive())
                K[3] += p->force(a.r + L[2]);
        K[3] *= stepSize / a.mass;
        L[3] = (a.v + K[2])*stepSize;
        
        newR = a.r + (L[0] + L[1]*2 + L[2]*2 + L[3])/6.;
        newV = a.v + (K[0] + K[1]*2 + K[2]*2 + K[3])/6.;
        
        if (newR.x > 1) {
            
        }
        
        // add kinetic and potential energy
        E0 = 0.5*(a.v * a.v);
        E = 0.5*(newV*newV);
        for (const auto& p : potentials) {
            if (p->isActive()) {
                E0 += p->potential(a.r) / a.mass;
                E += p->potential(newR) / a.mass;
            }
        }
        
        // E0 should never be 0, but check just in case
        // If E0 = 0, make sure that this current step is not repeated
        dE = E0 != 0 ? abs((E-E0)/E0) : tolerance*0.999;
        
        if (dE < tolerance) {
            timePassed += stepSize;
            a.r = newR;
            a.v = newV;
        }
        
        // Calculate new stepsize
        stepSize = dE != 0 ? 0.9*stepSize*pow(tolerance/dE, 0.25) : 10*stepSize;
        if (stepSize == 0 or !isfinite(stepSize)) stepSize = 1e-8;
        // If the current stepsize would surpass timeStep,
        // adjust it to keep time synchronous
        if ((uint64_t)(round((timePassed + stepSize)*1E9)) > timeStep_in_ns) {
            a.dtMove = stepSize;
            stepSize = timeStep - timePassed;
        }
    }
}
