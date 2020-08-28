//
//  dsmc.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 25/02/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef dsmc_h
#define dsmc_h

#include <array>                    // array container
#include <cmath>                    // standard math library
#include <filesystem>               // filesystem
namespace fs = std::filesystem;
#include <future>                   // promise + future for asynchronous code
#include <list>                     // list container
#include <memory>                   // for smart pointers
#include <random>                   // random number generators and distributions
#include <tuple>                    // tuple container for multiple return values
#include <vector>                   // vector container

#include "active-inactive.hpp"      // time dependent de-/activation of potentials
#include "atom.h"                   // struct for one particle
#include "json.hpp"                 // open source json library for configuration
using json = nlohmann::json;
#include "log.h"                    // definition of logging class and status bar
#include "potential.h"              // potential factory and base class
#include "rungeKutta.h"             // RK4 method with energy conservation
#include "threadpool.hpp"           // threadpool class
#include "trajectory.hpp"           // time evolution functionality for parameters
#include "vector3d.h"               // own simple 3-vector struct with operators


// Constants
constexpr auto u = 1.66054E-27;         // kg
constexpr auto kB = 1.38065E-23;        // J/K
constexpr auto hbar = 1.05457E-34;      // Js
constexpr auto A85 = 84.91179;
constexpr auto A87 = 86.90918;
constexpr auto a_0 = 5.29177E-11;       // m
constexpr auto NSEC = 1.E-9;            // s
constexpr auto USEC = 1.E-6;            // s
constexpr auto PI = 3.14159265;

/// A type for a random number generator.
typedef std::mt19937 MyRandomGenerator;

/// A type to store the atomic cloud.
/// @discussion This is the main storage type for the simulation
typedef std::vector<atom> Cloud;

/// A type to store smaller parts of the cloud.
/// @discussion This is used to access the atomic cloud for the collision procedure
typedef std::vector<atom*> Cell;


/// number of threads to run
constexpr unsigned int nThreads = 64;

/// Choosing a unit for time propagation (micro-/nanosecond)
constexpr double timeUnit = USEC;


/// Main simulation class.
/// @discussion Everything happens within this class
class DSMC
{
public: // member variables
    
    
    
// ---------------------------------------------------
//          COMMAND LINE OPTIONS/ARGUMENTS
//
    /// Static variable to select randomness.
    /// @discussion If set to true, all random number generators will be properly seeded.
    /// If set to false, results will repeat with each run
    static bool randomize;
    
    /// Static variable to state if the cloud will be animated.
    /// @discussion If set to true, the positions of all atoms in the cloud will be saved to disk after each (?) iteration.
    /// @warning The time intervals will not be regular.
    static bool animate;
    
    /// Static variable that contains the time of each frame
    /// @discussion Is set by the fps parameter of the animate option
    static double frametime;
    
    /// Static variable that is used to decice if already present output at the same location should be overwritten.
    static bool overwrite;

    /// Static variable that is used to continue running even if the phase space density has reached 2.6124
    static bool continueForHighPSD;
    static constexpr double PSD_Limit = 1.001;

    /// Static variable that determines whether TAS should be used or not
    static bool use_subCells;
    
    
    /// Static variable containing the file path to the configuration file.
    static fs::path configFile;
    /// Static variable containing the folder path to the output directory.
    static fs::path outputFolder;
    /// path variable containing the folder path to the animation output
    fs::path animationFolder;
    
// ---------------------------------------------------
    
    
    
    /// Enum class for different initial density distribution.
    enum class DensDistr {
        Harmonic, /**< enum value for a cloud starting in a harmonic distribution */
        Uniform   /**< enum value for a cloud starting in a homogeneous distribution */
    };
    
    /// Enum class for different initial velocity distribution.
    enum class VelDistr {
        MaxwellBoltzmann, ///< enum value for a cloud starting in a regular MB distribution
        Delta             ///< enum value for a cloud starting with only the most probable speed
    };
    
    /// Instance of logging class defined in log.h
    Log log;
    
    
    
private: // member variables
    
    

    /// The wanted number of collision pairs tested in each step per thread.
    /// @discussion The value is chosen to avoid errors stemming from rounding.
    /// E.g. if 50.4 pairs should get chosen, the rounding error is 0.8%.
    /// With 100 pairs per thread, the rounding error has an upper bound of 0.5%.
    static constexpr int _minPairsPerCell = 100;
    
    /// Number of particles per subcell, will only be used is the use of subcells is enabled
    static constexpr int _nPerSubCells = 2;
    
    /// Ratio of dt/tau where dt is the time step and tau the mean expected collision time.
    /// @discussion In DSMC, the time step should be << than the mean coll. time.
    static constexpr double _ratio_timeStep_meanCollisionTime = 0.02;
    
    /// Threshold for the doubling procedure (1 macro-particle -> 2 macro-particles)
    static constexpr int _nDoublingThreshold = 5000;

    /// Maximum relative change in one parameter
    static constexpr double _maxParameterChange = 0.0001;
    
// ---------------------------------------------------
//                  FIXED PARAMETERS
//
    /// Mass of the atomic species.
    double _mass;
    
    /// Maximum simulated duration (in seconds).
    double _Duration;
// ---------------------------------------------------

    
    
// ---------------------------------------------------
//               STORAGE CONTAINERS
//
    /// Main storage for all atoms in the simulation
    Cloud cloud;
    
    /// the number of cells in x, y, z direction
    std::array<unsigned int, 3> _cellStructure;
    /// basis vector of the cell grid (diagonal of one cell)
    vector3d _cellBasis;
    /// minimum values for the coordinates
    vector3d _cellOrigin;
    
    /// Average number of particles per cell
    std::vector<double> _averageParticlesPerCell;
    
    /// Storage for atom pointers for the collisions
    std::vector<Cell> cells;
// ---------------------------------------------------

    
    
// ---------------------------------------------------
//               INITIAL PARAMETERS
// (most of these are not used while running
//  the simulation, only for the setup)
//
    /// Initial number of particles (not atoms)
    uint64_t _nPart;
    /// Initial multiplicity factor
    int _initMultiplicity;
    
    /// Initial cloud shape
    DensDistr _cloudShape;
    /// Side length/radius of the initial cloud
    double _cloudLengthDimension;
    
    /// Initial velocity distribution
    VelDistr _velocityDist;
    /// Maxwell-Boltzmann defined temperature of the cloud
    double _temperature;
//
// ---------------------------------------------------
    
    
    
// ---------------------------------------------------
//                OTHER PARAMETERS
//
//
    /// unique identifier for the atoms
    uint64_t _atomcounter;
    
    /// maximum value of (cross section * relative velocity) over all possible atom pairs.
    double _maxProbability;
    /// mean value of (cross section * relative velocity) over all possible atom pairs
    double _meanProbability;
    
    /// Volume of the cloud.
    double _domainVolume;
    /// Volume of one cell.
    double _cellVolume;
    
    /// RNGs for the collision procedure
    std::vector<MyRandomGenerator> _collisionRNGs;
    
    /// Closest subcells (euclidean distance)
    std::vector<std::array<int,3>> _3d_coords_sorted;
    
    
    /// Time that has passed since the beginning of the simulation (in units of timeUnit)
    int _TIME;
    /// Current step size (in units of timeUnit)
    uint64_t _timeStep;
    
    /// Current multiplicity factor / statistic weight.
    /// @discussion Describes how many real atoms are represented by one particle in the simulation at any given time.
    int _Multiplicity;
    
    /// If the animate option is set, this will count the number of frames already output
    unsigned int _frameCounter;
// ---------------------------------------------------
    
    
    
// ---------------------------------------------------
//                     POTENTIALS
//
    /// Dictionary that stores the IDs of potentials being used and their corresponding indices in other, number indexed containers
    dict<int> _potentialIDs;
    
    /// Storage for all potentials in the simulation.
    /// @discussion Unique pointers are used due to the factory design pattern for the potentials themselves.
    std::vector<std::unique_ptr<Potential>> _potentials;
    
    /// List of dictionaries for the time evolution of each potential parameter.
    /// @discussion Each dictionary corresponds to one potential.
    /// @warning using a list is crucial if cross-compatibility is desired. MSVC doesn't compile vectors of maps of non-copyable objects. *rolls eyes aggressively*
    std::list<dict<Trajectory>> _potentialTrajectories;
    
    /// Time evolution of whether potentials are active or not.
    std::vector<ActiveSwitch> _potentialSwitches;
// ---------------------------------------------------
    
    
    
// ---------------------------------------------------
//                 SCATTERING LENGTH
//
    /// Current value of the scattering length (in SI units)
    double _scatteringLength;
    /// Time evolution of the scattering length
    Trajectory _scatteringLengthTrajectory;
// ---------------------------------------------------

    
    
// ---------------------------------------------------
//                       LOSSES
//
    /// Loss coefficients for 1-body (background), 2-body and 3-body losses.
    /// @discussion Units: 1/s; m^3/s; m^6/s
    std::array<double, 3> _lossCoefficients;
    
    /// Time evolution of loss coefficients.
    /// @warning cannot use an array here because there is no default constructor for Trajectory objects.
    std::vector<Trajectory> _lossCoefficientTrajectories;
// ---------------------------------------------------
    
    
    
    
// ---------------------------------------------------
//                 LOGGED STATISTICS
// (These don't play a role in the simulation
//  itself, but they are logged anyway)
    
    /// Total number of collisions that happen within the simulation.
    unsigned long long _totalCollisions;
    
    /// Total number of particle pairs that get tested for collisions
    unsigned long long int _nPairs;
    /// number of pairs selected in each cell
    std::vector<double> _nPairsCells;
    /// number of collisions in each cell
    std::vector<int> _nCollCells;
// ---------------------------------------------------
    
    
    

private: // methods
    
// ---------------------------------------------------
//                   INITIALISATION
//
    /// Reads configuration from a json file.
    /// @discussion Initialises most member variables.
    /// @param config const std::filesystem::path reference that contains the path to the configuration file
    void readConfig(const fs::path& config);
    
    /// Set the starting value for _maxProbability and _meanProbability.
    /// @discussion Loops over all possible particle pairs and determines the maximum and mean value of the product of cross section and relative velocity.
    void set_initialProbabilities();
    
    /// Prepares the simulation for running.
    void init();
// -------------- END INITIALISATION -----------------
    
    

// ---------------------------------------------------
//                 HELPER FUNCTIONS
//
    // sorted alphabetically
    
    /// Gives the integer index of the cell in which the provided position lies.
    /// @param r a const vector3d reference of the position of the atom whose index shall be calculated.
    /// @returns the integer valued index of the cell.
    unsigned int cellIndex(const vector3d &r) const;
    
    /// Processes one collision event.
    /// @param a1 a pointer to the first particle
    /// @param a2 a pointer to the second particle
    /// @param local_maxProbability the value by which to divide the actual product of cross section * rel. velocity
    /// @param this_Probability a double reference to the cross section / relative velocity product of the particle pair
    /// @param generator a reference to the random number generator currently used
    /// @returns True if the collision happened, False otherwise
    bool collisionEvent(atom* a1, atom* a2, const double& local_maxProbability, double& this_Probability, MyRandomGenerator& generator) const;
    
    /// Calculates the cross section of a collision with a specific relative velocity.
    /// @discussion The cross section returned is not dependent on the angular momentum, it only includes s-wave scattering.
    /// @param vRel a const double reference, the relative velocity in m/s
    /// @returns a double containing the cross section of the collision in m^2
    double crossSection(const double &vRel) const;
    
    /// Calculates the mean collision time of the particles in the cloud.
    /// @returns a double containing the mean collision time in seconds
    double getMeanCollisionTime() const;
    
    /// Calculate only the mean velocity of all particles
    /// @returns a double containing the mean velocity in m/s
    double get_Mean_Velocity() const;
    
    /// Calculates an estimate for the cloud temperature
    /// @discussion Uses the translational temperature definition
    /// @returns a double value containing the temperature estimate in K.
    double get_temperatureEst() const;
    
    /// Calculates the sum of kinetic and potential energy of all atoms in the cloud
    /// @returns a double containing the total energy of all atoms in the cloud in J
    double get_totalEnergy() const;
    
    /// Convenience function for the number of atoms in the cloud
    /// @returns the number of atoms currently in the cloud
    unsigned int NAtoms() const { return (unsigned int)(_Multiplicity * NParticles()); }
    
    /// Convenience function for the number of particles in the cloud
    /// @returns the number of particles currently in the cloud
    unsigned int NParticles() const { return (unsigned int)(cloud.size()); }
    
    /// Function template to output the complete list of particles.
    /// @discussion writes the current time in seconds and the full position and velocity (optional) of all atoms to the stream.
    /// @param stream a reference to either a regular output stream or a output file stream.
    /// @param withVelocity include the velocities or not
    template<typename T>
    void print_Cloud(T& stream, bool withVelocity = false)
    {
        stream << "t=" << _TIME*timeUnit << "\n";
        stream << "M=" << _Multiplicity << "\n";
        stream << "rx,ry,rz";
        if (withVelocity)
            stream << ",vx,vy,vz\n";
        else
            stream << "\n";
        for (const auto &a : cloud) {
            stream << a.r.x << "," << a.r.y << "," << a.r.z;
            if (withVelocity)
                stream << "," << a.v.x << "," << a.v.y << "," << a.v.z << "\n";
            else
                stream << "\n";
        }
    }
    
    /// Generates a unit vector with random direction
    /// @param generator a reference to the random number generator currently used
    /// @returns a vector3d struct that has length 1. and a random orientation
    vector3d randomDirection(MyRandomGenerator& generator) const;
    
    /// Updates the cloud volume
    /// @discussion Defined as the particle number divided by the average density
    void update_domainVolume();
    
    /// Updates all potential parameters, the scattering length and loss coefficients.
    /// @param currentTime a const double reference to the value of time (in s)
    void update_parameters(const double& currentTime);

// ------------------- END HELPERS -------------------
    

    
// ---------------------------------------------------
//               ESSENTIAL DSMC STEPS
//
    /// Doubling procedure.
    /// @discussion Doubles the number of particles present in the cloud while halving the multiplicity.
    /// Every atom is mirrored in position and velocity along the z axis. The z values are not mirrored in case gravity is active, which would violate energy conservation.
    /// @warning This procedure is only energy conserving if all potentials themselves are symmetric under x and y reversal.
    /// @returns True if a doubling has happened, False otherwise (meaning the multiplicity is already at 1)
    bool doubling();
    
    /// Move all particles according to newtonian physics
    /// @discussion Uses the Runge Kutta 4 method to propagate all particles using Newtonian mechanics and the forces presented by the active potentials over the timespan of the current time step.
    void moveStep();
    
    /// Apply boundary conditions
    /// @discussion Check for all particles if they moved outside of the boundary presented by the active potentials and remove them from the cloud.
    void boundaries();
    
    /// Apply loss mechanisms
    /// @discussion Use loss coefficients to determine a number of atoms that should be lost within the timespan of the current timestep and remove these randomly.
    /// @warning considers only the average density. Density dependent losses are not implemented for different parts of the cloud.
    void losses();
    
    /// Fill the cells with pointers to atoms from the cloud.
    /// @discussion Find the cell number for each atom and add pointers to the respective cell
    void fillCells();
    
    /// Collide particles with each other.
    /// @discussion Uses the DSMC collision rules in each cell.
    /// Determines the timestep for the next iteration based on the new value for _maxProbability.
    void collisionStep();
    
    /// A type that carries the return values of the single cell collision
    typedef std::tuple<int, int, double, double> singleCell_return;
    
    /// Apply collision mechanics within a single cell.
    /// @discussion Uses probabilistic rules derived from classical collision rate equations to collide particles with eachother, using only s-wave collisions (relies on low temperatures)
    /// @param cell const reference to a Cell object, containing pointers to atoms that will be tested for collisions
    /// @param RNG a reference to the random number generator for this cell.
    /// @param returnTuple a reference to a tuple with the number of pairs selected in the cell, the number of collisions that happened in the cell and the maximum and mean value of cross section * relative velocity
    void collideSingleCell(const Cell &cell, std::promise<singleCell_return>& returnTuple, MyRandomGenerator& RNG);
    
// ------------ END ESSENTIAL DSMC STEPS -------------
    
public:
    /// Constructor for a DSMC object
    /// @discussion Reads the configuration file and runs the init function.
    /// This way, after creating a DSMC object, the simulation is immediately able to run.
    /// @param config a const std::filesystem::path reference to the configuration file
    DSMC(const fs::path& config = DSMC::configFile);
    
    /// Main loop of the simulation
    void run();
    
    
/*
// --------------- UNUSED -----------------
    
    /// Calculate only the mean position of all particles
    /// @returns a vector3d containing the mean position or centre of mass of all particles
    vector3d get_Mean_Position() const;
    
    /// Calculate only the RMS position of all particles
    /// @returns a vector3d containing the RMS radii of the cloud
    vector3d get_RMS_Position() const;
    
    /// Calculates both the mean and root-mean-square position of all particles
    /// @returns a tuple of vector3d structs, the first being the center of mass, the second a vector of RMS radii in all directions in m
    std::tuple<vector3d, vector3d> get_Mean_RMS_Position() const;
    
    /// Calculate only the RMS velocity of all particles
    /// @returns a double containing the RMS velocity in m/s
    double get_RMS_Velocity() const;
    
    /// Calculates both the mean and root-mean-square velocity of all particles
    /// @returns a tuple of double values, the first being the mean, the second the RMS value in m/s
    std::tuple<double, double> get_Mean_RMS_Velocity() const;
 
*/
};


#endif /* dsmc_h */
