#ifndef model_factory_h
#define model_factory_h

#include <cmath>            // mathematics library
#include <exception>        // exception safety
#include <map>              // dictionary functionality
#include <memory>           // smart pointers
#include <string>           // standard library string implementation


/// A template 'typedef' for a dictionary-like map
/// @discussion Uses the standard functionality of a std::map but binds the first type to a string
template <typename T>
using dict = std::map<std::string, T>;

class Model;
/// A typedef for a function pointer to create a Model
typedef std::unique_ptr<Model> (*create_fn)(double, double, double, const dict<double>&);


/// Base class for models
/// @discussion Provides boiler plate code and basic functionality.
class Model
{
protected:
    /// y value at the start of the interval
    double y0;
    /// y value at the end of the interval
    double y1;
    /// length of the interval in seconds
    double T;
    /// Map from string to double reference
    /// @discussion Is used to address parameters of the child classes
    dict<double&> map = {};
    
    /// Set all parameters in constructor
    /// @discussion In the constructor of a model, all parameters have to be set once
    /// @param pars a const reference to a dictionary of doubles which contains the set of parameters
    /// @warning throws an out of range exception if a parameter is missing
    void init(const dict<double>& pars) {
        std::string msg = "No value was provided for the initialisation of '";
        for (const auto&[name,reference] : map) {
            try { reference = pars.at(name); }
            catch(const std::out_of_range& oor) {
                throw std::invalid_argument(msg + name + "'.\n");
            }
        }
    }
    
    /// Base constructor for a model
    /// @param y0 a double that contains the y value at the beginning of the interval
    /// @param y1 a double that contains the y value at the end of the interval
    /// @param T a double that contains the length of the interval in seconds
    /// @param map a const reference to a map of double references that contains names and addresses of parameters
    /// @param pars a const reference to a map that contains any additional parameters
    Model(double y0, double y1 = NAN, double T = 0, const dict<double&>& map = {}, const dict<double>& pars = {}) :
    y0(y0), y1(y1), T(T), map(map)
    {
        init(pars);
    }
    /// Potentials are non copy-able
    Model(const Model&) = delete;
public:
    /// Child classes can have their own parameters so the destructor has to be virtual
    virtual ~Model() {}
    /// Pure virtual function to evaluate the model at a certain time
    virtual double eval(double t) const = 0;
};

/// Exponential model
class Exp : public Model
{
private:
    /// Internal parameter: time constant
    double tc;
    
public:
    /// Constructor for an exponential model
    /// @param y0 a double that contains the y value at the start of the interval
    /// @param y1 a double that contains the y value at the end of the interval
    /// @param T a double that contains the length of the interval in seconds
    /// @param pars a const reference to a dictionary that contains the value of the time constant
    Exp(double y0, double y1, double T, const dict<double>& pars) : Model(y0, y1, T, {{"tc", tc}}, pars) {}
    
    /// Virtual constructor for this model
    /// @returns a unique pointer to the created model
    static std::unique_ptr<Model> Create(double y0, double y1, double T, const dict<double>& pars) { return std::make_unique<Exp>(y0, y1, T, pars); }
    
    double eval(double t) const {
        // Exponential with a time constant tc
        double x = (y1 - y0)/(exp(-T/tc) - 1.);
        return x*(exp(-t/tc) - 1.) + y0;
    }
};

/// Linear model
class Lin : public Model
{
public:
    /// Constructor for a linear model
    /// @param y0 a double that contains the y value at the start of the interval
    /// @param y1 a double that contains the y value at the end of the interval
    /// @param T a double that contains the length of the interval in seconds
    Lin(double y0, double y1, double T) : Model(y0, y1, T) {}
    
    /// Virtual constructor for this model
    /// @returns a unique pointer to the created model
    static std::unique_ptr<Model> Create(double y0, double y1, double T, const dict<double>&) { return std::make_unique<Lin>(y0, y1, T); }
    
    double eval(double t) const {
        return (y1 - y0)*(t/T) + y0;
    }
};


/// Switch model
class Switch : public Model
{
public:
    /// Constructor for a switch model
    /// @param y0 a double that contains the y value during the interval
    /// @param y1 a double that contains the y value after switching at the end of the interval
    /// @param T a double that contains the length of the interval in seconds
    Switch(double y0, double y1, double T) : Model(y0, y1, T) {}
    
    /// Virtual constructor for this model
    /// @returns a unique pointer to the created model
    static std::unique_ptr<Model> Create(double y0, double y1, double T, const dict<double>&) { return std::make_unique<Switch>(y0, y1, T); }
    
    double eval(double t) const { return t >= T ? y1 : y0; }
};

/// Factory to create and use models
/// @discussion This is not the base class for the model itself, just the class that provides a 'virtual' constructor.
/// This is also intrinsically a singleton, meaning the class is just a wrapper and the object cannot be instantiated more than once.
class Models
{
private:
    /// Map of strings to function pointers
    /// @discussion Is used to provide access to the registered model child classes

    dict<create_fn> m_map;
    /// Constructor registers some basic models immediately */
    Models() {
        m_map["lin"] = &Lin::Create;
        m_map["exp"] = &Exp::Create;
        m_map["switch"] = &Switch::Create;
    }
    /// Singleton: No copy constructor */
    Models(const Models&) = delete;
    /// Assignment operator: Will only provide a reference to the singular instance */
    Models& operator = (const Models&) { return *this; }
    
    
    /// Internal implementation of the register function.
    /// @discussion Adds the creation function for a specific model to the map
    /// @param ID a const string reference to the public name that this model will be accessible by.
    /// @param fn a create_fn that contains the address of the create() function of the respective model
    /// @see Register()
    void I_Register(const std::string& ID, create_fn fn) { m_map[ID] = fn; }


    /// Internal implementation of the use function.
    /// @discussion Creates a new model from the given public name and parameters
    /// @param ID a const string reference to the public name of the model that shall be created
    /// @param y0 a double that contains the y value at the beginning of the interval
    /// @param y1 a double that contains the y value at the end of the interval
    /// @param T a double that contains the length of the interval in seconds
    /// @param pars a const reference to a  map that contains any additional parameters
    /// @see Use()
    /// @warning throws an out_of_range exception if the requested ID is not registered.
    /// Does NOT catch any exceptions thrown by the creation function
    /// @returns a unique pointer to the created model
    std::unique_ptr<Model> I_Use(const std::string& ID, double y0, double y1, double T, const dict<double>& pars) const {
        // Check if the ID exists in the m_map dictionary
        if (m_map.find(ID) != m_map.end()) {
            // Create a model of type ID with the parameters passed
            return m_map.at(ID)(y0, y1, T, pars);
        } else {
            // Throw an exception if the requested ID was not found
            throw std::out_of_range(std::string("'") + ID + "' is not a known model type.\n");
        }
    }
    
public:
    /// Destructor
    ~Models() { m_map.clear(); }
    
    /// A static method to create and use the singular instance of this class
    /// @discussion As a singleton, the only instance of this factory is constructed statically.
    /// @returns a reference to the static instance
    static Models& Get() { static Models m; return m; }
    
    /// Public implementation of the register function.
    /// @discussion Adds the creation function for a specific model to the map
    /// @param ID a const string reference to the public name that this model will be accessible by.
    /// @param fn a create_fn that contains the address of the create() function of the respective model
    /// @see I_Register()
    static void Register(const std::string& ID, create_fn fn) { return Get().I_Register(ID, fn); }


    /// Public implementation of the use function.
    /// @discussion Creates a new model from the given public name and parameters
    /// @param ID a const string reference to the public name of the model that shall be created
    /// @param y0 a double that contains the y value at the beginning of the interval
    /// @param y1 a double that contains the y value at the end of the interval
    /// @param T a double that contains the length of the interval in seconds
    /// @param pars a const reference to a map that contains any additional parameters
    /// @see I_Use()
    /// @warning throws an out_of_range exception if the requested ID is not registered.
    /// Does NOT catch any exceptions thrown by the creation function
    /// @returns a unique pointer to the created model
    static std::unique_ptr<Model> Use(const std::string& ID, double y0, double y1, double T, const dict<double>& pars = {}) { return Get().I_Use(ID, y0, y1, T, pars); }
};

#endif // model_factory_h
