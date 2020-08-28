//
//  rungeKutta.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 27/02/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef rungeKutta_h
#define rungeKutta_h

#include <memory>           // smart pointers
#include <vector>           // vector container

#include "atom.h"           // struct for one particle
#include "potential.h"      // potential factory and base class

/// Calculates the new position and velocity of a particle after a certain time
/// @discussion Propagates the particle passed to the function using one or multiple RK4 steps until timeStep is reached
/// Energy conservation is guaranteed to a certain accuracy by having an internal timestep that changes if the change in energy surpasses a certain maximum
/// @param a a reference to the particle that shall be moved
/// @param potentials a const reference to a vector container of unique pointers to the potentials that are relevant for the propagation
/// @param timeStep a double representing the timestep in s
void rungeKuttaStep(atom& a, const std::vector<std::unique_ptr<Potential>>& potentials, double timeStep);

#endif /* rungeKutta_h */
