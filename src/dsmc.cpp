//
//  DSMC.cpp
//  EvaporationSim
//
//  Created by Andreas von Haaren on 25/02/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#include <iostream>
#include <thread>
#include <fstream>
#include <string>

#include "dsmc.h"           // Definitions for this class
#include "Timer.h"          // Timing of function evaluation
#include "typename.hpp"     // Method to get a string representation of a type at runtime

// Initialize static members
bool DSMC::randomize = false;
bool DSMC::animate = false;
double DSMC::frametime = 0;
bool DSMC::overwrite = false;
bool DSMC::continueForHighPSD = false;
bool DSMC::use_subCells = false;
fs::path DSMC::outputFolder;
fs::path DSMC::configFile;

/* ---------------------- CONSTRUCTOR ---------------------- */

// Constructor default initializes all members
DSMC::DSMC(const fs::path& config) :
cloud(),
_cellStructure({4,4,4}),
cells(),
_nPart(0),
_initMultiplicity(1),
_cloudShape(DSMC::DensDistr::Uniform),
_velocityDist(DSMC::VelDistr::MaxwellBoltzmann),
_cloudLengthDimension(0.),
_Duration(0),
_mass(0),
_atomcounter(1), // 0 means empty (unsigned),
_meanProbability(0),
_maxProbability(0),
_temperature(0),
_domainVolume(0),
_cellVolume(0),
_TIME(0),
_timeStep(0),
_totalCollisions(0),
//m_totalExpectedCollisions(0),
_nPairs(0),
_scatteringLength(0),
_scatteringLengthTrajectory(0),
_lossCoefficients({0.,0.,0.}),
log()
{
    // Read, check and apply all settings
    {
        Timer t = Timer(std::string("Importing settings from ") + config.filename().string());
        readConfig(config);
    }
    // Fill the cloud and compute initial time step
    // and other variables
    {
        Timer t = Timer("Initialisation");
        init();
    }
}

/* -------------------- END CONSTRUCTOR -------------------- */





/* -------------------- INITIALISATION --------------------- */

void DSMC::readConfig(const fs::path& config)
{
    // Assume success at first
    // If reading any value fails, this will be set to false
    bool success = true;
    
    // import the configuration file as a json object
    std::ifstream f(config);
    json j;
    f >> j;
    
    // Friendly names for the initial density distribution
    std::map<std::string, std::pair<DSMC::DensDistr, std::string>> shapes = {
        {"uniform", std::pair<DSMC::DensDistr, std::string>(DSMC::DensDistr::Uniform, "side_length")},
        {"harmonic", std::pair<DSMC::DensDistr, std::string>(DSMC::DensDistr::Harmonic, "radius")}
    };
    
    std::string context = "";
    
    // Lambda to read a value from a json
    // Handles exceptions and logs errors
    auto readValue = [this, &context]<typename T>(const json& j, const std::string& key, T& value) {
        try {
            value = j.at(key).get<T>();
        } catch (const json::type_error& err) {
            log.Error((context + "Value of '" + key + "' should be of type ").append(type_name<T>()) + ".\n");
            return false;
        } catch (const json::out_of_range& err) {
            log.Error(context + "No entry given for '" + key + "'.\n");
            return false;
        }
        return true;
    };
    
    
    
    if (not readValue(j, "runtime", _Duration)) success = false;
    if (not readValue(j, "mass", _mass)) success = false;
    if (not readValue(j, "cell_structure", _cellStructure)) success = false;
    
    json start;
    if (not readValue(j, "starting_conditions", start)) {
        success = false;
    } else {
        context = "Initial conditions: ";
        if (not readValue(start, "atom_number", _nPart)) success = false;
        if (not readValue(start, "temperature", _temperature)) success = false;
        if (not readValue(start, "scattering_length", _scatteringLength)) success = false;
        if (success) _scatteringLengthTrajectory = Trajectory(_scatteringLength);
        json cloud;
        if (not readValue(start, "cloud", cloud)) {
            success = false;
        } else {
            context = "Initial conditions - cloud: ";
            std::string shape;
            if (not readValue(cloud, "shape", shape)) {
                success = false;
            } else {
                if (shapes.find(shape) == shapes.end()) {
                    log.Error(context + "'" + shape + "' is not a valid cloud shape.\n");
                    success = false;
                } else {
                    context = "Initial conditions - cloud (" + shape + "): ";
                    _cloudShape = shapes[shape].first;
                    if (not readValue(cloud, shapes[shape].second, _cloudLengthDimension)) success = false;
                }
            }
            std::string velocityDistribution;
            if (cloud.find("velocity_distribution") != cloud.end()) {
                if (readValue(cloud, "velocity_distribution", velocityDistribution)) {
                    if (velocityDistribution != "delta" && velocityDistribution != "maxwell-boltzmann") {
                        log.Error(context + "'" + velocityDistribution + "' is not a valid velocity distribution.\nThe only valid options are 'maxwell-boltzmann' or 'delta'.\n");
                        success = false;
                    } else {
                        _velocityDist = velocityDistribution == "delta" ? VelDistr::Delta : VelDistr::MaxwellBoltzmann;
                    }
                } else {
                    success = false;
                }
            } else {
                _velocityDist = VelDistr::MaxwellBoltzmann;
            }
        }
        
        context = "Initial conditions: ";
        std::array<std::string, 3> losses_names = {"K1","K2","K3"};
        if (start.find("losses") != start.end()) {
            // losses are not a necessary entry, default is 0
            for (int i = 0; i < 3; ++i) {
                try {
                    _lossCoefficients[i] = start["losses"].value(losses_names[i], 0.);
                    _lossCoefficientTrajectories.emplace_back(_lossCoefficients[i]);
                } catch (const json::type_error& e) {
                    log.Error("Initial conditions - losses: Entries should be numbers.");
                    success = false;
                }
            }
        } else {
            for (int i = 0; i < 3; ++i) {
                _lossCoefficients[i] = 0.;
                _lossCoefficientTrajectories.emplace_back(0.);
            }
        }
    }
    
    std::vector<json> pot;
    std::vector<ActiveSwitch> temporarySwitches;
    if (not readValue(j, "potentials", pot)) {
        success = false;
    } else {
        for (const auto& p : pot) {
            bool thisSuccess = true;    // temporary success boolean
            std::string potential_type;
            std::string potential_id;
            if (not readValue(p, "type", potential_type)) {
                success = false;
                thisSuccess = false;
            } else {
                context = "Potential (type '" + potential_type + "'): ";
                bool active;
                try {
                    active = p.value("active", true);
                    if (not readValue(p, "id", potential_id)) {
                        success = false;
                        thisSuccess = false;
                    } else {
                        if (_potentialIDs.find(potential_id) != _potentialIDs.end()) {
                            log.Error(context + "Every potential needs a unique identifier. '" + potential_id + "' is already used.\n");
                            success = false;
                            thisSuccess = false;
                        }
                    }
                    dict<double> pars;
                    if (potential_type == "gravity") {
                        pars["mass"] = _mass;
                    } else if (not readValue(p, "parameters", pars)) {
                        success = false;
                        thisSuccess = false;
                    }
                    
                    if (thisSuccess) {
                        _potentials.push_back(Potentials::Use(potential_type, active, pars));
                        _potentialIDs[potential_id] = (int) _potentials.size() - 1;
                        dict<Trajectory> temp;
                        for (const auto&[name,value] : pars) temp.emplace(name, value);
                        _potentialTrajectories.push_back(std::move(temp));
                        _potentialSwitches.emplace_back(active);
                        temporarySwitches.emplace_back(active);
                    }
                } catch (const std::out_of_range& e) {
                    log.Error(context + e.what());
                    log.Error(context + "Couldn't be created.\n");
                    success = false;
                } catch (const json::type_error& e) {
                    log.Error(context + "'active' must be either true or false.\n");
                    log.Error(context + "Couldn't be created.\n");
                    success = false;
                }
            }
        }
    }
    
    // Lambda to read the variables neccessary to
    // add a model to a trajectory,
    // create the model and add it
    auto readAndAddModel = [&readValue, this](const json& j, double tEnd, const std::string& keyName, Trajectory& t) {
        bool returnValue = true;
        double start = std::nan(0), end;
        std::string modelType;
        dict<double> pars = {};
        if (not readValue(j, "end", end)) returnValue = false;                  // final y value
        if (j.find("start") != j.end())
            start = j.at("start").get<double>();                                // start y value (optional)
        if (not readValue(j, "model", modelType)) returnValue = false;          // type of the model (e.g. lin, exp)
        pars = j.value("parameters", dict<double>());                           // dictionary with additional parameters
        try {
            t.addModel(modelType, tEnd, end, pars, start);
        } catch (const std::invalid_argument& inv) {
            log.Error(std::string("Timepoint t = ") + std::to_string(tEnd) + "; " + keyName + ": " + inv.what());
            returnValue = false;
        } catch (const std::out_of_range& oor) {
            log.Error(std::string("Timepoint t = ") + std::to_string(tEnd) + "; " + keyName + ": " + oor.what());
            return false;
        }
        
        
        return returnValue;
    };
    
    context = "";
    std::vector<json> trajectories;
    std::vector<double> times({0.});
    try {
        trajectories = j.value("timepoints", std::vector<json>());
    } catch (const json::type_error& e) {
        log.Error("The 'timepoints' parameter should contain a list of json objects.\n");
        success = false;
    }
    int trajectoriesIndex = -1;
    
    for (const auto& timepoint : trajectories) {
        ++trajectoriesIndex;
        double tEnd;
        context = "";
        if (not readValue(timepoint, "time", tEnd)) {
            success = false;
            log.Error(std::string("No time provided at position ") + std::to_string(trajectoriesIndex) + " in timepoints list.\n");
            break;
        } else {
            context = "Timepoint t = " + std::to_string(tEnd);
            if (tEnd < times[times.size() - 1]) {
                success = false;
                log.Error(std::string("Time provided at position ") + std::to_string(trajectoriesIndex) + " in timepoints list is lower than the last. Make sure that all points are in the correct order.\n");
                break;
            }
            
            // Scattering length trajectory
            if (timepoint.find("scattering_length") != timepoint.end()) {
                if (not readAndAddModel(timepoint.at("scattering_length"), tEnd, "Scattering Length", _scatteringLengthTrajectory)) success = false;
            }
            
            
            // Losses trajectories
            if (timepoint.find("losses") != timepoint.end()) {
                json l = timepoint.at("losses");
                std::array<std::string, 3> k = {"K1", "K2", "K3"};
                for (int i = 0; i < 3; ++i) {
                    if (l.find(k[i]) != l.end()) {
                        if (not readAndAddModel(l.at(k[i]), tEnd, std::string("Losses (") + k[i] + ")", _lossCoefficientTrajectories.at(i))) success = false;
                    }
                }
            }
            
            
            // Potential trajectories
            if (timepoint.find("potentials") != timepoint.end()) {
                dict<json> ps;
                if (not readValue(timepoint, "potentials", ps)) success = false;
                else {
                    for (auto&[id,p] : ps) {
                        if (_potentialIDs.find(id) == _potentialIDs.end()) {
                            log.Error(context + ": There is no potential in the simulation that has the ID '" + id + "'.\n");
                            success = false;
                            continue;
                        } else {
                            int index = _potentialIDs.at(id);
                            if (p.contains("active")) {
                                try {
                                    bool active = p.value("active", true);
                                    if (temporarySwitches[index].check(tEnd) != active) {
                                        temporarySwitches[index].switchAt(tEnd);
                                        _potentialSwitches[index].switchAt(tEnd);
                                    }
                                } catch (const json::type_error& e) {
                                    success = false;
                                    log.Error(context + "; ID='" + id + "': 'active' has to be a boolean value.\n");
                                }
                                
                                p.erase("active"); // All other parameters should be numbers after this
                            }
                            try {
                                dict<json> parameters = p.get<dict<json>>();
                                for (const auto& parameter : parameters) {
                                    try {
                                        if (not readAndAddModel(p.at(parameter.first), tEnd, std::string("ID='") + id + "', parameter='" + parameter.first + "'", (*std::next(_potentialTrajectories.begin(), index)).at(parameter.first)))
                                            success = false;
                                    } catch (const std::out_of_range& oor) {
                                        log.Error(context + ": '" + parameter.first + "' is not a valid parameter for potential '" + id + "'.\n");
                                    }
                                }
                            } catch (const json::type_error& e) {
                                log.Error(context + " Parameters for potential '" + id + "' are in an incorrect format.\n");
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (not success) {
        std::cout << std::endl;
        log.Error("There were errors while reading the config. The simulation is unable to run...\n");
        exit(1);
    }
}

void DSMC::set_initialProbabilities()
{
    _maxProbability = 0;
    _meanProbability = 0;
    uint64_t N = 0;
    double prob = 0;
    double vRel = 0;
    vector3d vRelVector = {0,0,0};
    
    // compute cross section * relative velocity for every possible pair of atoms
    for (auto it = cloud.begin(); it != cloud.end(); ++it) {
        for (auto jt = it+1; jt != cloud.end(); ++jt, ++N) {
            vRelVector = (*it).v - (*jt).v;
            vRel = sqrt(vRelVector * vRelVector);
            prob = crossSection(vRel) * vRel;
            _meanProbability += prob;
            if (prob > _maxProbability) _maxProbability = prob;
        }
    }
    
    _meanProbability /= N;
}


void DSMC::init()
{
    // Instantiate a random engine and initialise it
    MyRandomGenerator generator(randomize * std::random_device{}());
    
    // set the initial multiplicity
    while (_nPart > 4*_nDoublingThreshold) {
        _nPart /= 2;
        _initMultiplicity *= 2;
    }
    _Multiplicity = _initMultiplicity;
    
    
    // Fill the cloud according to the initial distribution
    cloud.reserve(_nPart);
    vector3d r,v;                                                                   // position, velocity
    std::normal_distribution<double> vel(0., sqrt(_temperature*kB/_mass));        // maxwell-boltzmann distribution
    double vMostProbable = sqrt(2*_temperature*kB/_mass);
    switch(_cloudShape){
        case DSMC::DensDistr::Harmonic: {
            std::normal_distribution<double> pos(0., _cloudLengthDimension);
            for (int i=0; i < _nPart; ++i) {
                r = {pos(generator), pos(generator), pos(generator)};
                v = _velocityDist == VelDistr::MaxwellBoltzmann ? vector3d(vel(generator), vel(generator), vel(generator)) : randomDirection(generator)*vMostProbable;
                cloud.emplace_back(r,v,0,_mass,_atomcounter);
                ++_atomcounter;
            }
            break;
        }
        case DSMC::DensDistr::Uniform: {
            std::uniform_real_distribution<double> pos(-_cloudLengthDimension/2., _cloudLengthDimension/2.);
            for (int i=0; i < _nPart; ++i) {
                r = {pos(generator), pos(generator), pos(generator)};
                v = _velocityDist == VelDistr::MaxwellBoltzmann ? vector3d(vel(generator), vel(generator), vel(generator)) : randomDirection(generator)*vMostProbable;
                cloud.emplace_back(r,v,0,_mass,_atomcounter);
                ++_atomcounter;
            }
            break;
        }
    }
    
    // apply boundary conditions once
    boundaries();
    // log initial atom number
    log.Info("Starting with ");
    *log.logger << NAtoms() << " atoms\n";
    // log the capture rate of the potential
    log.Info("Capture Rate: ");
    *log.logger << 100.*NAtoms()/((double)_nPart*_initMultiplicity) << " %\n\n";
    // If the boundary conditions reduced the number of particles,
    // apply the doubling procedure until it is above the threshold
    while (cloud.size() < _nDoublingThreshold) {
        if (not doubling()) {
            log.Warn("Almost no atoms left: ");
            *log.logger << cloud.size()*_Multiplicity << "\n";
            break;
        }
    }
    
    
    _averageParticlesPerCell.resize(_cellStructure[0] * _cellStructure[1] * _cellStructure[2]);
    cells.resize(_cellStructure[0] * _cellStructure[1] * _cellStructure[2]);
    _nPairsCells = std::vector<double>(_cellStructure[0] * _cellStructure[1] * _cellStructure[2], 0.);
    _nCollCells = std::vector<int>(_cellStructure[0] * _cellStructure[1] * _cellStructure[2], 0);
    // fill cells once for losses in first iteration
    fillCells();
    
    // Calculate the initial domain volume
    update_domainVolume();
    
    // Calculate the value of the maximum of
    // cross section * relative velocity over all particle pairs
    log.Info("Calculating initial Maximum Collision Probability ...\n");
    set_initialProbabilities();
    log.Info("Done. Maximum Collision Probability is ");
    *log.logger << _maxProbability << " m^3/s" << std::endl;
    log.Info("Expected mean collision time is ");
    int tCollMean = static_cast<int>(round(1E6*getMeanCollisionTime()));
    *log.logger << tCollMean << " us" << std::endl;
    
    // Calculate the initial time step
    _timeStep = (int)round(_ratio_timeStep_meanCollisionTime * getMeanCollisionTime() / timeUnit);
    log.Info("The initial time step is ");
    *log.logger << _timeStep << " us\n" << std::endl;
    
    // initialise the movement timestep for all particles
    for (auto& a : cloud) {
        a.dtMove = _timeStep*timeUnit;
    }
    
    // Initialise random number generators for the collision procedure
    for (int i = 0; i < cells.size(); ++i) {
        _collisionRNGs.emplace_back(i + randomize * std::random_device{}());
    }
    
    
    // Create appropriate output folders
    // At this point, the user has either agreed to overwriting the contents of the output folder
    // Or the overwrite option is set anyway
    fs::create_directory(outputFolder);
    fs::copy(configFile, outputFolder / "config.json", fs::copy_options::overwrite_existing);
    if (animate) {
        animationFolder = outputFolder / "animation";
        _frameCounter = 0;
        // remove all items from the animation folder if it already exists
        if (!fs::create_directory(animationFolder)) {
            for (auto& path : fs::directory_iterator(animationFolder)) {
                fs::remove_all(path);
            }
        }
    } else {
        if (fs::exists(outputFolder / "animation")) {
            fs::remove_all(outputFolder / "animation");
        }
    }
    
    
    // Generate 3d integer coordinates, sorted by euclidean distance
    // there are at most 4*doublingThreshold particles present.
    // we want ... particles per subcell
    // if all particles are in one cell (not gonna happen but whatever, lets make sure)
    // there will be this many subcells -> cube root in each direction
    // e.g for doubling at 5000 and 2 atoms per subcell, N = 22. Then the filling and sorting
    // will take < 5ms.
    if (use_subCells) {
        int N = std::pow(4*_nDoublingThreshold/_nPerSubCells, 1./3) + 1;
        for (int a = -N; a < N; ++a)
            for (int b = -N; b < N; ++b)
                for (int c = -N; c < N; ++c)
                    _3d_coords_sorted.push_back({a,b,c});
        
        std::sort(_3d_coords_sorted.begin(), _3d_coords_sorted.end(),
                  []( const std::array<int, 3>& lhs, const std::array<int, 3>& rhs )
        { return lhs[0]*lhs[0] + lhs[1]*lhs[1] + lhs[2]*lhs[2] < rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]; }
                  );
    }
}

/* ------------------ END INITIALISATION ------------------- */





/* ------------------------ HELPERS ------------------------ */

inline unsigned int DSMC::cellIndex(const vector3d &r) const
{
    return ((int)((r.z - _cellOrigin.z)/_cellBasis.z))*_cellStructure[0]*_cellStructure[1] +
           ((int)((r.y - _cellOrigin.y)/_cellBasis.y))*_cellStructure[0] +
            (int)((r.x - _cellOrigin.x)/_cellBasis.x);
}

bool DSMC::collisionEvent(atom* a1, atom* a2, const double& local_maxProbability, double &prob, MyRandomGenerator& generator) const
{
    std::uniform_real_distribution<> randomDouble(0.,1.);
    vector3d vREL = a1->v - a2->v;          // relative velocity of a1 and a2
    double v = sqrt(vREL * vREL);           // absolute value of vREL
    double RAND = randomDouble(generator);  // random number between 0 and 1
    
    prob = crossSection(v)*v;
    
    if (prob/local_maxProbability < RAND)   // check if a collision happens
        return false;
    
    vector3d vCOM = (a1->v + a2->v)/2.;     // centre of mass velocity
    vREL = randomDirection(generator)*v;
    
    a1->v = vCOM + vREL/2.;                 // new velocity of a1
    a2->v = vCOM - vREL/2.;                 // new velocity of a2
    
    
    // update last collision partner
    a1->lastCollided = a2->ID;
    a2->lastCollided = a1->ID;
    
    return true;
}

inline double DSMC::crossSection(const double &vRel) const
{
    double k = _mass/hbar/2.*vRel;     // wave vector
    double kSquared = k*k;
    double aSquared = _scatteringLength * _scatteringLength;
    return 8*PI*aSquared/(1.+ aSquared * kSquared);
}

inline double DSMC::getMeanCollisionTime() const
{
    return _domainVolume/(NAtoms()*_meanProbability);
}


inline double DSMC::get_Mean_Velocity() const
{
    double mean = 0;
    for (const auto& a : cloud) {
        mean += sqrt(a.v * a.v);
    }
    return mean/NParticles();
}

inline double DSMC::get_temperatureEst() const
{
    double v_squared = 0;
    for (const auto& a : cloud) {
        v_squared += a.v * a.v;
    }
    return v_squared / (3. * NParticles()) * (_mass / kB);
}

double DSMC::get_totalEnergy() const
{
    double Ekin = 0, Epot = 0;
    for (const auto& a : cloud) {
        Ekin += 0.5 * _mass * (a.v * a.v);
        for (const auto& p : _potentials)
            if (p->isActive())
                Epot += p->potential(a.r);
    }
    return Ekin + Epot;
}

inline vector3d DSMC::randomDirection(MyRandomGenerator& generator) const
{
    std::uniform_real_distribution<> randomDouble(0.,1.);
    
    double RAND1 = randomDouble(generator), RAND2 = randomDouble(generator);
    double phi = 2*PI*RAND1;
    double cosTheta = 2*RAND2 -1;
    double sinTheta = sqrt(1.-cosTheta*cosTheta);
    return vector3d(cosTheta, sinTheta*cos(phi), sinTheta*sin(phi));
}

inline void DSMC::update_domainVolume()
{
    uint64_t temp = 0;
    for (const auto& cell : cells) {
        temp += cell.size() * cell.size();
    }
    
    _domainVolume = (_cellVolume * NParticles() * NParticles()) / temp;
}

void DSMC::update_parameters(const double& currentTime)
{
    // Update the scattering length
    _scatteringLength = _scatteringLengthTrajectory.eval(currentTime);
    
    // Update the loss coefficients
    for (int i = 0; i < 3; ++i)
        _lossCoefficients[i] = _lossCoefficientTrajectories[i].eval(currentTime);
    
    // Update the parameters for every potential
    auto it = _potentialTrajectories.begin();
    for (int i = 0; i < _potentialIDs.size(); ++i, ++it) {
        dict<double> newParameters;
        for (const auto&[name, trajectory] : *it) {
            newParameters[name] = trajectory.eval(currentTime);
        }
        
        _potentials[i]->setParameters(newParameters);
        _potentials[i]->setActive(_potentialSwitches[i].check(currentTime));
    }
    
    fillCells();
    update_domainVolume();
}

/* ---------------------- END HELPERS ---------------------- */





/* -------------------- MAIN DSMC STEPS -------------------- */

bool DSMC::doubling()
{
    if (_Multiplicity == 1) {
        // If every particle already represents a single atom,
        // the cloud cannot be doubled any further
        return false;
    } else {
        // Loop over the 'old' cloud and mirror every atom
        vector3d rNew, vNew;
        int N0 = (int)cloud.size();
        auto it = cloud.begin();
        
        for (int i = 0; i < N0; ++i) {
            // Only mirror x and y in case gravity is on
            rNew = (*it).r; vNew = (*it).v;
            rNew.x *= -1; rNew.y *= -1;
            vNew.x *= -1; vNew.y *= -1;
            cloud.emplace_back(rNew, vNew, (*it).dtMove, (*it).mass, _atomcounter, 0);
            ++_atomcounter;
            ++it;
        }
        
        _Multiplicity /= 2;

        return true;
    }
}

void DSMC::moveStep()
{
    constexpr int nThreadsMove = nThreads;                                // number of threads to use
    
    static ThreadPool threadpool(nThreadsMove);                     // threads to use for the movement
    static std::array<std::future<void>, nThreadsMove> futures;     // Void futures for each thread to synchronise the program
    
    
    // function to run in each thread
    auto moveThreadFn = [this](int index, int max_index, const double& timeStep) {
        // Do one Runge-Kutta step for each atom in this part of the cloud
        for (int i = index; i < max_index; ++i) {
            rungeKuttaStep(cloud[i], _potentials, timeStep);
        }
    };
    
    unsigned int len = NParticles() / nThreadsMove;         // number of particles in every thread
    double dt = _timeStep * timeUnit;                       // timeStep in seconds
    
    // execute the thread function on separate parts of the cloud
    for (int i = 0; i < nThreadsMove; ++i)
        futures[i] = threadpool.push(moveThreadFn, i * len, (i + 1) * len, dt);
    
    moveThreadFn(nThreadsMove * len, NParticles(), dt); // move the rest of the cloud
    for (auto& f : futures)
        f.get(); // Synchronize
}

void DSMC::boundaries()
{
    // lambda that checks whether an atom should stay or be removed
    auto removalCondition = [this](atom& a) {
        bool withinBoundary = false;
        // an atom is considered at captured if ANY of the potential
        // withinBoundary functions returns true for its position
        for (const auto& pot : _potentials) {
            // only check a potential if it is not inactive
            if (pot->isActive()) {
                withinBoundary = withinBoundary or pot->withinBoundary(a.r);
            }
        }
        return (not withinBoundary);
    };
    
    // remove all atoms from the cloud that are outside the potential boundaries
    std::erase_if(cloud, removalCondition);
}

void DSMC::losses()
{
    fillCells();
    // Loss probabilities for each cell
    static std::vector<double> m_lossProbabilities_Cells;
    
    double density;
    m_lossProbabilities_Cells.resize(cells.size());
    // calculate the loss probability in each cell
    for (int i = 0; i < cells.size(); ++i) {
        density = cells[i].size() / _cellVolume * _Multiplicity;
        m_lossProbabilities_Cells[i] = _timeStep*timeUnit*(
                                                           _lossCoefficients[0] +
                                                           _lossCoefficients[1] * density +
                                                           _lossCoefficients[2] * density * density
                                                           );
    }
    
    // Initialise a RNG for the losses
    static MyRandomGenerator gen(randomize*std::random_device{}());
    static std::uniform_real_distribution<> rand(0.,1.);

    // Remove a particle if the loss probability is greater
    // than a random number between 0 and 1
    std::erase_if(cloud, [this](const atom& a) { return m_lossProbabilities_Cells[a.cellIndex] > rand(gen); });
}

inline void DSMC::fillCells()
{
    // Remove all atoms from the cells first
    for (auto& cell : cells) {
        cell.clear();
        cell.reserve(2*NParticles()/cells.size());
    }
    static int counter = 0;
    
    // Get the outermost coordinates
    vector3d max(0,0,0);
    _cellOrigin = max;
    for (const auto& a : cloud) {
        if (a.r.x > max.x) max.x = a.r.x;
        if (a.r.x < _cellOrigin.x) _cellOrigin.x = a.r.x;
        if (a.r.y > max.y) max.y = a.r.y;
        if (a.r.y < _cellOrigin.y) _cellOrigin.y = a.r.y;
        if (a.r.z > max.z) max.z = a.r.z;
        if (a.r.z < _cellOrigin.z) _cellOrigin.z = a.r.z;
    }
    // expand by a little to make sure every atom is inside and not at the edge
    _cellOrigin *= 1.001;
    max *= 1.001;
    
    // get the basis vector
    _cellBasis = max - _cellOrigin;
    _cellBasis.x /= _cellStructure[0];
    _cellBasis.y /= _cellStructure[1];
    _cellBasis.z /= _cellStructure[2];
    
    _cellVolume = _cellBasis.x * _cellBasis.y * _cellBasis.z;
    
    // Check the index for each atom and
    // add it to the corresponding cell
    for (auto& a : cloud) {
        auto index = cellIndex(a.r);
        a.cellIndex = index;
        cells[index].push_back(&a);
    }
    for (int i = 0; i < cells.size(); ++i) {
        _averageParticlesPerCell[i] = (_averageParticlesPerCell[i]*counter + cells[i].size()) / (counter + 1);
    }
    
    ++counter;
}

void DSMC::collisionStep()
{
    // threadpool that is kept in memory
    static ThreadPool threadpool(nThreads);
    
    // counters for the number of collision pairs chosen
    // and the number of actually happened collisions
    int nPairs = 0, collisionCount = 0;
    double new_meanProbability = 0;
    
    static std::vector<std::future<singleCell_return>> futures(cells.size());
    std::vector<std::promise<singleCell_return>> promises(cells.size());
    
    for (auto& x : _nCollCells) x = 0;
    
    for (int i = 0; i < cells.size(); ++i) {
        // execute collideSingleCell for each cell
        futures[i] = promises[i].get_future();
        threadpool.push(&DSMC::collideSingleCell, this, std::ref(cells[i]), std::ref(promises[i]), std::ref(_collisionRNGs[i]));
    }
    
    std::chrono::time_point<std::chrono::steady_clock> tStart = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < cells.size(); ++i) {
        // Capture the return values from the single cells
        auto[nPairsResult, nCollResult, maxProbResult, meanProbResult] = futures[i].get();
        nPairs += nPairsResult;
        collisionCount += nCollResult;

        if (maxProbResult > _maxProbability)
            _maxProbability = maxProbResult;
        
        new_meanProbability += nPairsResult * meanProbResult;
        
    }
    new_meanProbability /= nPairs;
    _meanProbability = new_meanProbability;
    
    std::chrono::duration<float> runtime = std::chrono::high_resolution_clock::now() - tStart;
    

    
    // update the timeStep
    _timeStep = (int)round(_ratio_timeStep_meanCollisionTime * getMeanCollisionTime() / timeUnit);
    // time steps larger than 100 milliseconds should not exist
    if (_timeStep > 0.1 / timeUnit or _timeStep < 0) _timeStep = (int)(0.1 / timeUnit);

    // Logging of collision number
    _nPairs += nPairs;
    _totalCollisions += _Multiplicity*(unsigned long long)collisionCount;
    

    if (runtime.count() > 1) {
        log.Error("The collision procedure took too long. ");
        *log.logger << "(" << runtime.count() << " s)\n";
        log.Info("N = " + std::to_string(NAtoms()) + "\n");
        log.Info("max. Prob. = ");
        *log.logger << _maxProbability << "\n";
        log.Info("Cell counts: ");
        for (const auto& cell : cells)
            *log.logger << cell.size() << " ";
        *log.logger << "\n";
        log.Info("Pairs selected in this run: ");
        *log.logger << nPairs;
        std::cout << std::endl;
        //exit(1);
    }
}

void DSMC::collideSingleCell(const Cell &cell, std::promise<singleCell_return>& promise, MyRandomGenerator& generator)
{
    // If the cell is empty, has only one atom or
    // has two atoms but they are their last collision partner,
    // return from the function without doing anything
    if (cell.size() <= 1 or (cell.size() == 2 and cell[0]->lastCollided == cell[1]->ID and cell[1]->lastCollided == cell[0]->ID)) {
        promise.set_value({ 0,0,0.,0. });
        return;
    }
    
    /* ------------- subdivide the cell into multiple subcells ------------- */
    int n_subCells = 1; // number of subcells in each direction
    
    if (use_subCells) {
        n_subCells = (int)std::pow(cell.size()/_nPerSubCells, 1./3) + 1;
    }
    
    std::vector<Cell> subCells(n_subCells * n_subCells * n_subCells);
    vector3d boundLow;
    {
        int thisIndex = cell[0]->cellIndex;
        int x,y,z;
        z = thisIndex / _cellStructure[0] / _cellStructure[1];
        y = (thisIndex % (_cellStructure[0] * _cellStructure[1])) / _cellStructure[0];
        x = (thisIndex % (_cellStructure[0] * _cellStructure[1])) % _cellStructure[0];
        boundLow = _cellOrigin + vector3d(x * _cellBasis.x, y * _cellBasis.y, z * _cellBasis.z);
    }
    auto subCellIndex = [this,&boundLow,&n_subCells](const vector3d& r) {
        vector3d diff = r - boundLow;
        return (int)(diff.z / (_cellBasis.z / n_subCells)) * n_subCells * n_subCells +
               (int)(diff.y / (_cellBasis.y / n_subCells)) * n_subCells +
               (int)(diff.x / (_cellBasis.x / n_subCells));
    };
    for (auto& a : cell) {
        subCells[subCellIndex(a->r)].push_back(a);
    }
    /* -------------------------- end subdivision -------------------------- */
    
    
    // Random index inside the cell
    std::uniform_int_distribution<> cellRand(0, static_cast<int>(cell.size())-1);
    int nColl = 0; // number of collisions that took place
    
    
    double local_maxProbability = _maxProbability;
    
    // Number of particle pairs to check for a collision
    // Explanation in any textbook that covers DSMC
#define UseAverageParticlesPerCell 1
#if UseAverageParticlesPerCell
#define SecondMultiplier _averageParticlesPerCell[cell[0]->cellIndex]
#else
#define SecondMultiplier (cell.size() - 1)
#endif

    double nPairs_fp = _Multiplicity * cell.size() * SecondMultiplier * _maxProbability * _timeStep*timeUnit / (_cellVolume * 2.);
    
    _nPairsCells[cell[0]->cellIndex] += nPairs_fp;

    int nPairs = static_cast<int>(round(nPairs_fp));
    
    // we want at least ... pairs tested for statistical significance
    if (nPairs > _minPairsPerCell) {
        local_maxProbability = _maxProbability;
    } else {
        local_maxProbability = _minPairsPerCell / nPairs_fp * _maxProbability;
        nPairs = _minPairsPerCell;
    }
    
    double new_maxProbability = 0;     // maximum of cross section * relative velocity
    double mean_Probability = 0;        // mean value of cross section * relative velocity
    atom *a1, *a2;          // two atom pointers
    
    auto findCollisionPartner = [this,&subCells,&a1,&generator,&n_subCells](int a1_subCellIndex) {
        int coords_offset;
        int new_subCellIndex;
        std::vector<bool> searched(subCells.size(), false);
        for (int counter = 1, coords_index = 1; counter < subCells.size(); ++coords_index) {
            // coords_index = 0 is a 0,0,0 offset -> skip
            const auto& coords = _3d_coords_sorted[coords_index];
            
            // convert the 3 coordinates into a 1D index for the subCell vector
            coords_offset = coords[2] * n_subCells * n_subCells + coords[1] * n_subCells + coords[0];
            
            // the new subCellIndex must actually be within the subCell vector
            new_subCellIndex = a1_subCellIndex + coords_offset;
            if (new_subCellIndex >= 0 and new_subCellIndex < subCells.size()) {
                if (!searched[new_subCellIndex]) {
                    ++counter; // new subcell is being searched
                    searched[new_subCellIndex] = true;
                } else {
                    continue;
                }
                
                const auto& new_subCell = subCells[new_subCellIndex];
                if (new_subCell.size() > 0) {
                    if (new_subCell.size() == 1) {
                        if (new_subCell[0]->lastCollided != a1->ID or a1->lastCollided != new_subCell[0]->ID) {
                            return new_subCell[0];
                        }
                    } else {
                        // size is >1, there has to be a valid partner here!
                        std::uniform_int_distribution<> subCellRand(0, (int)new_subCell.size() - 1);
                        atom* partner = new_subCell[subCellRand(generator)];
                        while (partner->lastCollided == a1->ID and a1->lastCollided == partner->ID) {
                            partner = new_subCell[subCellRand(generator)];
                        }
                        return partner;
                    }
                } // subCell is empty -> next
            } // new Index is out of bounds for subCells
               
        }
        
        // this should never be reached
        return (atom*)nullptr; // auto lambda needs same return type in all cases
    };
    
    // Select nPairs pairs of atoms
    int i = nPairs;
    for (; i > 0; --i) {
        a1 = cell[cellRand(generator)]; // random atom from the entire cell
        // now we need to select another atom
        int thisSubIndex = subCellIndex(a1->r);
        const auto& currentSubCell = subCells[thisSubIndex];
        
        // if the cell has only two atoms, the check at the beginning
        // may succeed, but if they collide then, it will be stuck
        if (cell.size() == 2 and cell[0]->lastCollided == cell[1]->ID and cell[1]->lastCollided == cell[0]->ID) {
            break;
        }
        
        if (currentSubCell.size() > 1) {
            // try in the same subcell first.
            if (currentSubCell.size() == 2) {
                // if there are only two atoms, try those, but only if
                // they are not their last collision partner
                a2 = currentSubCell[0];
                if (a2->ID == a1->ID) a2 = currentSubCell[1];
                if (a2->lastCollided == a1-> ID and a1->lastCollided == a2->ID) {
                    // if a2 collided last with a1 and vice versa, it is not physical
                    // to let them collide again -> search a new partner
                    a2 = findCollisionPartner(thisSubIndex);
                }
            } else {
                // if there are more than 2 atoms present, one of them must be valid
                std::uniform_int_distribution<> subCellRand(0, (int)currentSubCell.size() - 1);
                a2 = currentSubCell[subCellRand(generator)];
                while (a2->ID == a1->ID or (a2->lastCollided == a1->ID and a1->lastCollided == a2->ID)) {
                    a2 = currentSubCell[subCellRand(generator)];
                }
            }
        } else {
            // if the current subcell has only one atom, we need to search around for a partner
            a2 = findCollisionPartner(thisSubIndex);
        }
        
        if (a2 == nullptr) {
            log.Error("The collision partner search returned a nullptr. Exiting...");
            exit(404);
        }
        
        double thisProb;    // value of cross section * relative velocity for this collision
        // check if a collision between a1 and a2 actually happens
        if (collisionEvent(a1, a2, local_maxProbability, thisProb, generator)) {
            ++nColl;
            _nCollCells[a1->cellIndex]++;
            if (thisProb > new_maxProbability)
                new_maxProbability = thisProb;
        }
        mean_Probability += thisProb;
    }
    mean_Probability /= (nPairs - i);
    
    promise.set_value({nPairs, nColl, new_maxProbability, mean_Probability});
}

/* ------------------ END MAIN DSMC STEPS ------------------ */






void DSMC::run()
{
    int runs = 0;  // counter for the number of steps
    
    namespace time = std::chrono;
    using clock = time::high_resolution_clock;
    using my_duration = time::duration<float>;
    my_duration tMOV(0), tCOL(0), tFILL(0), tLOSS(0), tDOUBLE(0), tBOUND(0), tREST(0);
    std::vector<std::pair<std::string, time::duration<float>&>> durations = {
        {"move atoms", tMOV},
        {"collide atoms", tCOL},
        {"fill cells", tFILL},
        {"calculate losses", tLOSS},
        {"perform the doubling", tDOUBLE},
        {"apply boundary conditions", tBOUND},
        {"do everything else", tREST}
    };
    
    
    std::ofstream output(outputFolder / "data.csv");
    output << "t,dt,N,V,Vc,T,PSD,NColl,multi,maxprob,mfp,nPair0,nPairC,nColl0,nCollC,nCollStep,NLoss,a,K1,K2,K3";
    
    
    if (animate) {
        std::ofstream ani(animationFolder / "0");
        print_Cloud(ani, true);
        ++_frameCounter;
    }
    
    auto it = _potentialTrajectories.begin();
    for (const auto&[id,_] : _potentialIDs) {
        for (const auto&[name,_] : *it) {
            output << "," << id + ":" + name;
        }
        output << "," << id + ":active";
        ++it;
    }
    output << "\n";
    
    fillCells();
    update_domainVolume();
    
    unsigned long lostFromTrajectory = 0;      // count the number of atoms that appear outside the boundary after the potential parameters are updated
    double lostFromTrajectory_Rel = 0.;        // the relative atom loss due to moving potential walls
    unsigned long maxParticles_lostFromTrajectory = 0;
    double maxRate_lostFromTrajectory = 0.;
    
    unsigned long lostFromLosses = 0;               // count the number of atoms lost from background and 2-/3-body collisions
    double lostFromLosses_Rel = 0.;                 // the relative atom loss from background and 2-/3-body collisions
    unsigned long n0;                               // temporary variables
    auto temp = get_temperatureEst();
    double lambda = sqrt(2*PI*hbar*hbar/(_mass*kB*temp));                      // thermal wavelength
    double PSD = NAtoms()/_domainVolume*std::pow(lambda,3); // phase space density
    double mfp = _domainVolume * get_Mean_Velocity()/(NAtoms() * _meanProbability); // mean free path
    
    double Kn = mfp / std::pow(_domainVolume, 1./3); // Knudsen number
    double min_Kn = Kn;
    
    double rCell_mfp = std::pow(_cellVolume, 1./3) / mfp; // ratio of cell size to mfp
    double max_rCell_mfp = rCell_mfp;
    
    log.Info("Initial temperature: " + std::to_string(temp*1e6) + " uK\n");
    log.Info("Initial PSD: " + std::to_string(PSD) + "\n");
    log.Info("Initial Knudsen number: " + std::to_string(Kn) + "\n");
    log.Info("Initial ratio of cell size to mean free path: " + std::to_string(rCell_mfp) + "\n");
    log.Info("Initial energy per particle: " + std::to_string(get_totalEnergy()/NParticles()/kB*1e6) + " kB * uK\n\n");
    double PSD0 = PSD;
    int N0 = NAtoms();
    double T0 = temp;
    double V0 = _domainVolume;
    
    std::vector<double> nPairsCells_Diff = std::vector(cells.size(), 0.);
    
//    _Duration = 1000*getMeanCollisionTime();
    
    double old_timeStep = _timeStep;
    double maxDev = 0;
    log.Info("Running Simulation... \n");
    { // begin statusBar
        statusBar statusBarInstance(0, _Duration);
    
        /* -------------------- MAIN SIMULATION LOOP -------------------- */
        time::time_point<time::steady_clock> tStart, tEnd;
        // run until m_Duration is reached, all atoms are gone or a phase space density of 1 is reached
        for (int i = 0; _TIME*timeUnit < _Duration and NParticles() >= _nDoublingThreshold and (PSD < PSD_Limit or continueForHighPSD) and _timeStep > 0; ++i, ++runs, statusBarInstance.update(_TIME * timeUnit)) {
            auto it = _potentialTrajectories.begin();
            maxDev = 0;
            for (int i = 0; i < _potentialIDs.size(); ++i, ++it) {
                if (!(_potentials[i]->isActive())) continue;
                double temp, comp;
                dict<double> oldParameters = _potentials[i]->getParameters();
                for (const auto& [name, trajectory] : *it) {
                    temp = trajectory.eval((_TIME + _timeStep) * timeUnit);
                    if (void(comp = std::abs((oldParameters[name] - temp) / oldParameters[name])), comp > maxDev) maxDev = comp;
                }
            }

            if (maxDev > _maxParameterChange) _timeStep = (int)std::round((_maxParameterChange / maxDev) * _timeStep);
            if (_timeStep == 0) _timeStep = 1;
            old_timeStep = _timeStep;


            tStart = clock::now();
            double currentTime = _TIME * timeUnit;
            auto temp = get_temperatureEst();
            
            lambda = sqrt(2*PI*hbar*hbar/(_mass*kB*temp));
            PSD = NAtoms()/_domainVolume*std::pow(lambda,3);
            
            mfp = _domainVolume * get_Mean_Velocity()/(NAtoms() * _meanProbability);
            if (void(Kn = mfp / std::pow(_domainVolume, 1./3)), Kn < min_Kn) min_Kn = Kn;
            if (void(rCell_mfp = std::pow(_cellVolume, 1./3) / mfp), rCell_mfp > max_rCell_mfp) max_rCell_mfp = rCell_mfp;
            
            
            
            int nCollSum = 0;
            for (const auto& x : _nCollCells) nCollSum += x;
            
            
            // Log parameters to file
            output << _TIME*timeUnit;
            output << "," << _timeStep*timeUnit;
            output << "," << NAtoms();
            output << "," << _domainVolume;
            output << "," << _cellVolume;
            output << "," << temp;
            output << "," << PSD;
            output << "," << _totalCollisions;
            output << "," << _Multiplicity;
            output << "," << _maxProbability;
            output << "," << mfp;
            output << "," << nPairsCells_Diff[0];
            output << "," << nPairsCells_Diff[cells.size()/2];
            output << "," << _nCollCells[0];
            output << "," << _nCollCells[1];
            output << "," << nCollSum;
            output << "," << lostFromLosses;
            output << "," << _scatteringLength;
            for (const auto& K : _lossCoefficients) output << "," << K;
            for (int i = 0; i < _potentialIDs.size(); ++i) {
                auto pars = _potentials[i]->getParameters();
                for (auto[_,value] : pars) output << "," << value;
                output << "," << (int)(_potentials[i]->isActive());
            }
            output << "\n";
            
            if (animate) {
                int nFramesToPrint = ((_TIME + _timeStep) * timeUnit - (_frameCounter - 1) * frametime) / frametime;
                for (int i = 0; i < nFramesToPrint; ++i, ++_frameCounter) {
                    std::ofstream ani(animationFolder / std::to_string(_frameCounter));
                    print_Cloud(ani, true);
                }
            }
            
            tEnd = clock::now();
            tREST += tEnd - tStart;

            
            // MOVE
            tStart = tEnd;
            moveStep();
            tEnd = clock::now();
            tMOV += tEnd - tStart;

            
            
            // BOUNDARIES
            tStart = tEnd;
            boundaries();
            tEnd = clock::now();
            tBOUND += tEnd - tStart;
            
            
            // UPDATE PARAMETERS
            tStart = tEnd;
            n0 = NAtoms();
            update_parameters(currentTime);
            tEnd = clock::now();
            tREST += tEnd - tStart;
            
            tStart = tEnd;
            boundaries();
            lostFromTrajectory += n0 - NAtoms();
            lostFromTrajectory_Rel += (n0 - NAtoms())/((double)(n0));
            maxRate_lostFromTrajectory = std::max(maxRate_lostFromTrajectory, (n0 - NAtoms())/(n0*_timeStep*timeUnit));
            maxParticles_lostFromTrajectory = std::max(maxParticles_lostFromTrajectory, n0/_Multiplicity - NParticles());
            tEnd = clock::now();
            tBOUND += tEnd - tStart;
            
            
            
            // LOSSES
            tStart = tEnd;
            n0 = NAtoms();
            losses();
            lostFromLosses += n0 - NAtoms();
            lostFromLosses_Rel += (n0 - NAtoms())/((double)(n0));
            tEnd = clock::now();
            tLOSS += tEnd - tStart;
            
            
            // DOUBLING
            tStart = tEnd;
            if (NParticles() < _nDoublingThreshold) {
                if (not doubling()) {
                    log.Warn("Minimum Atom Number reached");
                }
            }
            tEnd = clock::now();
            tDOUBLE += tEnd - tStart;
            
            
            
            // FILL CELLS
            tStart = tEnd;
            fillCells();
            tEnd = clock::now();
            tFILL += tEnd - tStart;

            
            // COLLISIONS
            auto old = _nPairsCells;
            tStart = tEnd;
            collisionStep();
            tEnd = clock::now();
            tCOL += tEnd - tStart;
            
            
            for (int i = 0; i < cells.size(); ++i) {
                nPairsCells_Diff[i] = _nPairsCells[i] - old[i];
            }
            
            
            _TIME += old_timeStep;
        }
        /* ------------------ END MAIN SIMULATION LOOP ------------------ */

    statusBarInstance.update(_Duration); // set status bar to 100%
    } // end statusBar
    
    
    temp = get_temperatureEst();
    mfp = _domainVolume * get_Mean_Velocity()/(NAtoms() * _meanProbability);
    lambda = sqrt(2*PI*hbar*hbar/(_mass*kB*temp));
    PSD = NAtoms()/_domainVolume*std::pow(lambda,3);
    
    int nCollSum = 0;
    for (const auto& x : _nCollCells) nCollSum += x;
    
    output << _TIME*timeUnit;
    output << "," << old_timeStep*timeUnit;
    output << "," << NAtoms();
    output << "," << _domainVolume;
    output << "," << _cellVolume;
    output << "," << temp;
    output << "," << PSD;
    output << "," << _totalCollisions;
    //                output << "," << m_totalExpectedCollisions;
    output << "," << _Multiplicity;
    output << "," << _maxProbability;
    output << "," << mfp;
    output << "," << nPairsCells_Diff[0];
    output << "," << nPairsCells_Diff[cells.size()/2];
    output << "," << _nCollCells[0];
    output << "," << _nCollCells[1];
    output << "," << nCollSum;
    output << "," << lostFromLosses;
    output << "," << _scatteringLength;
    for (const auto& K : _lossCoefficients) output << "," << K;
    for (int i = 0; i < _potentialIDs.size(); ++i) {
        auto pars = _potentials[i]->getParameters();
        for (auto[_,value] : pars) output << "," << value;
        output << "," << (int)(_potentials[i]->isActive());
    }
    output << "\n";
    
    output.close();
    
    
    double efficiency = -1 * std::log(PSD/PSD0) / std::log(NAtoms() / (double)(N0));


    log.Info("Simulated Time: " + std::to_string(_TIME*timeUnit*1000) + " ms\n");
    log.Info("Runs calculated: " + std::to_string(runs) + "\n");
    log.Info("The final time step is " + std::to_string(old_timeStep) + " us\n\n");

    log.Info("Inelastic atom loss: " + std::to_string(lostFromLosses) + " (" + std::to_string((100.*lostFromLosses)/(N0 - NAtoms())) + " % of all losses)\n");
    log.Info("Atom loss due to moving potential walls: " + std::to_string(lostFromTrajectory) + " (" + std::to_string((100.*lostFromTrajectory)/(N0 - NAtoms())) + " % of all losses)\n");
    log.Info("Maximum no. of particles lost in one step from moving potential walls: " + std::to_string(maxParticles_lostFromTrajectory) + "\n\n");
    
    log.Info("Minimum Knudsen number: " + std::to_string(min_Kn) + "\n");
    log.Info("Maximum ratio of cell size to mean free path: " + std::to_string(max_rCell_mfp) + "\n\n");
    
    log.Info("Number of pairs selected for collisions: " + std::to_string(_nPairs) + "\n");
    log.Info("Number of recorded collisions: " + std::to_string(_totalCollisions) + "\n\n");
/*
    log.Info("Number of recorded/expected collisions: " + std::to_string(m_totalCollisions) + "/" + std::to_string(m_totalExpectedCollisions) + "\n");
    log.Info("Relative deviation of collision number: " + std::to_string((100.*m_totalCollisions)/m_totalExpectedCollisions - 100) + " %\n\n");
*/
    
    int maxLen = 0;
    for (const auto&[name,_] : durations) {
        maxLen = std::max(maxLen, (int)name.size());
    }
    for (const auto&[name,duration] : durations) {
        std::string fill(maxLen - (int)name.size(), ' ');
        log.Info("Time to " + name + fill + ": " + std::to_string(duration.count()) + " s\n");
    }
    
    std::cout << std::endl;
    
    if (animate) {
        log.Info(std::to_string(_frameCounter) + " frames were output to " + animationFolder.string() + "\n\n");
    }
    
    
    log.Info("Particles remaining: " + std::to_string(NAtoms()) + "\n");
    log.Info("Final temperature: " + std::to_string(temp*1e6) + " uK\n");
    log.Info("Final PSD: " + std::to_string(PSD) + "\n");
    log.Info("Average evaporation efficiency: " + std::to_string(efficiency) + "\n");
    auto E = get_totalEnergy()/NParticles();
    log.Info("Final energy per particle: " + std::to_string(E/kB*1e6) + " kB * uK\n\n");
    
    double normedEfficiency = NAtoms() / (_TIME * timeUnit) * std::log(PSD / PSD0);
    
    json summary;
    summary["atom_number"] = {N0,NAtoms()};
    summary["temperature"] = {T0,temp};
    summary["energy_per_particle"] = E;
    summary["volume"] = {V0,_domainVolume};
    summary["phase_space_density"] = {PSD0,PSD};
    summary["runtime"] = _TIME * timeUnit;
    summary["inelastic_loss"] = lostFromLosses;
    summary["collisions"] = _totalCollisions;
    summary["efficiency"] = efficiency;
    summary["normed_efficiency"] = normedEfficiency;
    std::ofstream sum(outputFolder / "summary.json");
    sum << summary;
    sum.close();
}




/*
// UNUSED

inline vector3d DSMC::get_Mean_Position() const
{
    vector3d mean = 0;
    for (const auto& a : cloud) {
        mean += a.r;
    }
    return (mean/NParticles());
}

inline vector3d DSMC::get_RMS_Position() const
{
    vector3d rms = 0;
    for (const auto& a : cloud) {
        rms += {a.r.x * a.r.x, a.r.y * a.r.y, a.r.z * a.r.z};
    }
    rms /= NParticles();
    rms.x = sqrt(rms.x);
    rms.y = sqrt(rms.y);
    rms.z = sqrt(rms.z);
    return rms;
}

inline std::tuple<vector3d, vector3d> DSMC::get_Mean_RMS_Position() const
{
    vector3d rMean = 0, rRMS = 0;
    for (const auto &a : cloud) {
        rMean += a.r;
        rRMS += {a.r.x * a.r.x, a.r.y * a.r.y, a.r.z * a.r.z};
    }
    rMean /= NParticles();
    rRMS /= NParticles();
    rRMS.x = sqrt(rRMS.x);
    rRMS.y = sqrt(rRMS.y);
    rRMS.z = sqrt(rRMS.z);
    
    return {rMean, rRMS};
}

inline double DSMC::get_RMS_Velocity() const
{
    double rms = 0;
    for (const auto& a : cloud) {
        rms += a.v * a.v;
    }
    return sqrt(rms/NParticles());
}

inline std::tuple<double, double> DSMC::get_Mean_RMS_Velocity() const
{
    double vMean = 0, vRMS = 0;
    double vSq;
    for (const auto &a : cloud) {
        vSq = a.v * a.v;
        vMean += sqrt(vSq);
        vRMS += vSq;
    }
    vMean /= NParticles();
    vRMS = sqrt(vRMS/NParticles());
    
    return {vMean, vRMS};
}
 
 */
