//
//  potential.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 27/02/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef potential_h
#define potential_h

#include <cmath>
#include <exception>            // exception safety
#include <map>                  // dictionary functionality
#include <memory>               // smart pointers
#include <optional>             // optional return values
#include <stdexcept>            // standard exceptions
#include <string>               // standard library string implementation
#include <vector>               // vector container

#include "vector3d.h"           // own simple 3-vector struct with operators

/// A template alias for a dictionary-like map
/// @discussion Uses the standard functionality of a std::map but binds the first type to a string
template <typename T>
using dict = std::map<std::string, T>;

class Potential;
/// A typedef for a function pointer to create a Potential
typedef std::unique_ptr<Potential> (*create_potential_fn)(bool, const dict<double>&);


/// Factory to create and use potentials
/// @discussion This is not the base class for the potential itself, just the class that provides a 'virtual' constructor.
/// This is also intrinsically a singleton, meaning the class is just a wrapper and the object cannot be instantiated more than once.
class Potentials
{
private:
    /// Empty constructor
    Potentials() {}
    /// Singleton: No copy constructor
    Potentials(const Potentials&) = delete;
    /// Assignment operator: Will only provide a reference to the singular instance
    Potentials& operator =(const Potentials&) { return *this; }
    
    /// Map of strings to function pointers
    /// @discussion Is used to provide access to the registered potential child classes
    dict<create_potential_fn> m_map;
    
    /// Internal implementation of the register function.
    /// @discussion Adds the creation function for a specific potential to the map
    /// @param ID a const string reference to the public name that this potential will be accessible by.
    /// @param fct_ptr a create_potential_fn that contains the address of the create() function of the respective potential
    /// @see Register()
    void I_Register(const std::string& ID, create_potential_fn fct_ptr) { m_map[ID] = fct_ptr; }
    
    /// Internal implementation of the IsRegistered() function
    /// @discussion Checks if a specific name is registered
    /// @param ID a const string reference to the ID searched for
    /// @returns True, if the ID exists in the map, False otherwise
    /// @see IsRegistered()
    bool I_IsRegistered(const std::string& ID) const { return m_map.find(ID) != m_map.end(); }
    
    /// Internal implementation of the use function.
    /// @discussion Creates a new potential from the given public name and parameters
    /// @param ID a const string reference to the public name of the potential that shall be created
    /// @param active a bool value that sets the new potential to active or inactive
    /// @param pars a const reference to a dictionary that contains the parameters for the new potential
    /// @see Use()
    /// @warning throws an out_of_range exception if the requested ID is not registered.
    /// Does NOT catch any exceptions thrown by the creation function
    /// @returns a unique pointer to the created potential
    std::unique_ptr<Potential> I_Use(const std::string& ID, bool active, const dict<double>& pars) const {
        // Check if the ID exists in the m_map dictionary
        if (I_IsRegistered(ID)) {
            // Create a potential of type ID with the parameters passed
            return m_map.at(ID)(active, pars);
        } else {
            // Throw an exception if the requested ID was not found
            std::string msg = std::string("'") + ID + "' is not a known potential type.\n";
            throw std::out_of_range(msg);
        }
    }
public:
    /// Destructor
    ~Potentials() { m_map.clear(); }
    
    /// A static method to create and use the singular instance of this class
    /// @discussion As a singleton, the only instance of this factory is constructed statically.
    /// @returns a reference to the static instance
    static Potentials& Get() {
        static Potentials instance;
        return instance;
    }
    
    /// Public implementation of the register function
    /// @discussion Adds the creation function for a specific potential to the map
    /// @param ID a const string reference to the public name that this potential will be accessible by.
    /// @param fct_ptr a create_potential_fn that contains the address of the create() function of the respective potential
    /// @see I_Register()
    static void Register(const std::string& ID, create_potential_fn fct_ptr) { return Get().I_Register(ID, fct_ptr); }
    
    /// Public implementation of the IsRegistered() function
    /// @discussion Checks if a specific name is registered
    /// @param ID a const string reference to the ID searched for
    /// @returns True, if the ID exists in the map, False otherwise
    /// @see I_IsRegistered()
    static bool IsRegistered(const std::string& ID) { return Get().I_IsRegistered(ID); }
    
    /// Public implementation of the use function.
    /// @discussion Creates a new potential from the given public name and parameters
    /// @param ID a const string reference to the public name of the potential that shall be created
    /// @param active a bool value that sets the new potential to active or inactive
    /// @param pars a const dict<double> reference to a dictionary that contains the parameters for the new potential
    /// @see I_Use()
    /// @warning throws an out_of_range exception if the requested ID is not registered.
    /// Does NOT catch any exceptions thrown by the creation function
    /// @returns a unique pointer to the created potential
    static std::unique_ptr<Potential> Use(const std::string& ID, bool active, const dict<double>& pars) { return Get().I_Use(ID, active, pars); }
};


/// Base class for potentials
/// @discussion Provides boiler plate code and basic functionality.
/// Any child class has to implement functions for the potential itself, the force it gives as well as a function that checks if a point is within the boundary.
class Potential
{
protected:
    /// Potential is active or inactive
    bool m_active;
    /// Map from string to double reference
    /// @discussion Is used to address parameters of the child classes
    dict<double&> m_parameterMap;
    
    /// Set all parameters in constructor
    /// @discussion In the constructor of a potential, all parameters have to be set once
    /// @param pars a const reference to a dictionary of doubles which contains the set of parameters
    /// @warning throws an out of range exception if a parameter is missing
    void initialisePars(const dict<double>& pars) {
        std::string msg = "Missing value while creating potential: No value was provided for the parameter '";
        for (const auto&[name,reference] : m_parameterMap) {
            try { reference = pars.at(name); }
            catch(const std::out_of_range& oor) {
                throw std::out_of_range(msg + name + "'.\n");
            }
        }
    }
    
    /// Base Constructor for a potential
    /// @param active a boolean, the state of the potential at initialisation, default inactive
    /// @param parameterMap a const reference to a dictionary of double references that contains the names and addresses of the parameters
    /// @param pars a const reference to a dictionary of doubles that contains the values of the parameters
    Potential(bool active = false, const dict<double&>& parameterMap = {}, const dict<double>& pars = {}) :
    m_active(active), m_parameterMap(parameterMap)
    { initialisePars(pars); }
    
    /// Potentials are non copy-able
    Potential(const Potential&) = delete;
public:
    /// Child classes can have their own parameters so the destructor has to be virtual
    virtual ~Potential() {}
    
    /// Setter for the active value
    void setActive(bool active) { m_active = active; }
    
    /// Getter for the active value
    bool isActive() const { return m_active; }
    
    /// Setter for the numeric parameters
    /// @param pars a const reference to a dictionary of doubles containing the new parameter values
    /// @warning throws an invalid_argument exceptions if the passed dictionary contains unknown parameters that are not applicable to the potential
    void setParameters(const dict<double>& pars) {
        // counter for the number of valid parameters
        int changesApplied = 0;
        
        // Loop through the known parameters and check pars if there are any changes
        for (const auto&[name,reference] : m_parameterMap) {
            if (pars.find(name) != pars.end()) {
                reference = pars.at(name);
                ++changesApplied;
            }
        }
        
        // If there are keys in 'pars' that don't fit any of the known parameters,
        // throw an exception containing the names of the unknown parameters
        if (changesApplied != pars.size()) {
            std::string wrong_keys;
            for (const auto&[name,_] : pars) {
                if (m_parameterMap.find(name) == m_parameterMap.end())
                    wrong_keys += wrong_keys.size() == 0 ? name : std::string(", ") + name;
            }
            throw std::invalid_argument(wrong_keys);
        }
    }
    
    /// Getter for the numeric parameters
    /// @returns a dictionary of doubles with the names and values of the parameters
    dict<double> getParameters() const {
        dict<double> pars;
        for (const auto&[name,value] : m_parameterMap)
            pars[name] = value;
        return pars;
    }
    
    /// Pure virtual function for the potential at a certain point
    /// @param r a const reference to a vector3d, the position to evaluate
    /// @returns a double containing the value of the potential at the position r in J
    virtual double potential(const vector3d& r) const = 0;
    /// Pure virtual function for the force at a certain point
    /// @param r a const reference to a vector3d, the position to evaluate
    /// @returns a vector3d struct containing the force at the position r in N
    virtual vector3d force(const vector3d& r) const = 0;
    /// Pure virtual function for the boundary condition of the potential
    /// @param r a const reference to a vector3d, the position to evaluate
    /// @returns a boolean value, True if r is within the bounding box of the potential, False otherwise
    virtual bool withinBoundary(const vector3d& r) const = 0;
};

#endif /* potential_h */

