#ifndef active_inactive_h
#define active_inactive_h

#include <exception>            // exception safety
#include <string>               // to_string()
#include <vector>               // vector container


/// Similar to Trajectory but for the active value of potentials
class ActiveSwitch {
private:
    /// Vector of doubles that contains the times where the potential switches between active/inactive
    std::vector<double> switchingTimes;
    /// Boolean that selects whether the potential starts as active or inactive
    bool start;
    /// Current position in the switchingTimes vector
    int index;
public:
    /// Constructor for an ActiveSwitch
    /// @param start the value of active at time 0
    ActiveSwitch(bool start) : switchingTimes(), start(start), index(0) {}

    /// Function to add a time when to switch
    /// @param t a double containing the new switching time (in seconds)
    /// @warning throws an invalid argument exception if the new time value is in the wrong order
    void switchAt(double t) {
        // if the switchingTimes vector is empty, add the new value immediately
        if (not switchingTimes.size()) {
            switchingTimes.push_back(t);
            return;
        }

        if (switchingTimes[switchingTimes.size() - 1] > t) {
            // if the last value in the switchingTimes vector is bigger than the new one, throw an exception
            throw std::invalid_argument(std::to_string(t) + " is smaller than the latest switching time. Please add the switching times in the correct order.");
        } else {
            switchingTimes.push_back(t);
            return;
        }
         
    }
    
    /// Evaluation of the Switcher
    /// @param t a double containing the time to check
    /// @returns a boolean containing the value of active/inactive at time t
    bool check(double t) {
        // if the switchingTimes vector is empty, return the start value
        if (not switchingTimes.size()) 
            return start;
        
        // case that t is after the last switch
        if (index == switchingTimes.size()) {
            return index % 2 ? !start : start;
        }
        // propagate index until it's in the correct place in the switchingTimes vector
        while (switchingTimes[index] < t) {
            ++index;
            if (index == switchingTimes.size())
                return index % 2 ? !start : start;
        }
        return index % 2 ? !start : start;
    }
};


#endif // active_inactive_h
