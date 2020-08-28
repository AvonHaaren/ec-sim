#ifndef trajectory_h
#define trajectory_h

#include "model-factory.hpp"        // Interval based functions


/// A trajectory contains multiple models
/// @discussion This offers the possibility of functions defined over arbitrary intervals
class Trajectory
{
private:
    /// A vector of doubles that contains the end times of the intervals (in seconds)
    std::vector<double> m_tEnds;
    /// A vector on unique pointers to the models over the different intervals
    std::vector<std::unique_ptr<Model>> m_models;
public:
    /// Constructor for a trajectory object
    /// @discussion Even if the trajectory does not contain any models, it will still evaluate to a constant value
    /// @param y starting value of the trajectory
    Trajectory(double y) {
        m_tEnds.push_back(0);
        m_models.push_back(Models::Use("switch", y, y, 0));
    }
    
    /// Any trajectory is not copy-able
    Trajectory(const Trajectory&) = delete;
    /// Default move constructor
    Trajectory(Trajectory&&) = default;
    /// No copy assignment
    Trajectory& operator=(const Trajectory&) = delete;
    /// Default move assignment
    Trajectory& operator=(Trajectory&&) = default;
    
    /// Adding another interval
    /// @discussion Adds another interval with its function definition to the trajectory.
    /// If no parameter for the y value at the start of the interval is supplied, the last available y value in the trajectory is used, making the function piecewise continuos.
    /// @param model a const string reference with the ID of the model to add
    /// @param t1 a double containing the time where the new interval ends (in seconds)
    /// @param y1 a double containing the y value at the end of the interval
    /// @param pars a const reference to a map that contains optional parameters
    /// @param y0 a double (optional) containing the y value at the start of the interval
    /// @warning throws an invalid argument exception if the new interval is not at the end of the sequence
    void addModel(const std::string& model, double t1, double y1, const dict<double>& pars = {}, double y0 = std::nan(0)) {
        // start at same level where last model ended, if no value is provided
        if (isnan(y0)) y0 = eval(t1);
        // Starting time of the interval
        double t0 = m_tEnds[m_tEnds.size() - 1];
        
        // Check if the new function keeps the order correct
        if (t0 <= t1) {
            // Add the new model and the new end time to the vectors
            m_models.push_back(Models::Use(model, y0, y1, t1 - t0, pars));
            m_tEnds.push_back(t1);
        } else {
            throw std::invalid_argument(std::to_string(t1) + " is smaller than the latest timestamp in trajectory (" + std::to_string(t0) + ").\n");
        }
    }
    
    /// Evaluate the Trajectory at a specific point in time
    /// @param t a double that contains the time where the trajectory will be evaluated (in seconds)
    /// @returns a double containing the value of the trajectory at time t
    double eval(double t) const {
        // Find the interval in which t lies
        int i = 0;
        for (; t > m_tEnds[i]; ++i)
            if (i == m_tEnds.size())
                break;
        
        
        if (i == m_tEnds.size()) {
            // if t is not within any interval, return the last value of the last interval
            return m_models[i - 1]->eval(m_tEnds[i - 1] - m_tEnds[i - 2]);
        } else {
            // if t is in the first interval, the starting time of that interval is 0,
            // otherwise the starting time is the end time of the last interval
            double t0 = i == 0 ? 0 : m_tEnds[i - 1];
            // evaluate the model corresponding to the interval i
            return m_models[i]->eval(t - t0);
        }
    }
};



#endif //trajectory_h
