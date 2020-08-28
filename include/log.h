//
//  log.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 26/02/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef log_h
#define log_h

#include <iostream>         // input + output


/// Uninterruptible status bar.
/// @discussion Used to show the percentual propagation of some computation.
/// @warning If anything is printed to std::cout while the status bar is active, it will be broken.
struct statusBar
{
private:
    // Private variables to store fixed information and the current percentage value
    double _min, _max;
    int _percentage;
    
public:
    /// Constructor for a statusBar object.
    /// @discussion Assigns the internal minimum and maximum values.
    /// The percentage is set to -1 to indicate that no progress has been made.
    /// @param min the starting value of the variable whose progress is shown
    /// @param max the expected final value of the variable at the end of the calculation
    statusBar(double min, double max) : _min(min), _max(max), _percentage(-1) {}
    
    /// Destructor for a statusBar object.
    /// @discussion Ensures that a new line is printed when the status bar has finished.
    ~statusBar() { std::cout << "\n"; }
    
    /// Public member function to update the status bar
    /// @param currentValue the new value of the variable whose progress is shown
    void update(double currentValue)
    {
        // Restart the same line (carriage return)
        std::cout << "\r";
        
        // Calculate the new percentage with the current value
        int percentage = (int)(std::round(100.*(currentValue - _min)/(double)(_max - _min)));
        
        // To avoid unneccessary output, only update if the percentage value changes
        if (percentage == _percentage)
            return;
        
        // If the percentage has changed, update the internal storage
        _percentage = percentage;
        
        // Print the percentage value with the appropriate amount of spaces (single/double digit)
        if (_percentage / 100 == 0) {
            std::cout << " ";
            if (_percentage / 10 == 0) {
                std::cout << " ";
            }
        }
        
        // Status bar output in the form [|||||||||     ]
        std::cout << _percentage << " %  [";
        for (int i = 0; i < _percentage; ++i)
            std::cout << "|";
        for (int i = _percentage; i < 100; ++i)
            std::cout << " ";
        std::cout << "]";
        
        // Depending on compiler optimisation, the new output might not be flushed automatically
        std::cout.flush();
    }
};


/// A simple logging class.
class Log
{
public:
    
    /// An enum containing the different logging levels
    enum Level {
        LevelError, LevelWarning, LevelInfo
    };

private:
    /// The logging level
    Level _level;
    
public:
    /// Unelegant solution to print various variable to the output stream.
    /// @discussion Provides operator << functionality.
    std::ostream* logger;
    
    /// Constructor for Log class.
    /// @param level the logging level to be used. Default is the highest level (all output)
    /// @param ostream pointer to an output stream
    Log(Level level = LevelInfo, std::ostream *ostream = &std::cout) noexcept :
    _level(level), logger(ostream) {}
    
    /// Public method to change the logging level if desired.
    /// @param level the new logging level
    void SetLevel(Level level) { _level = level; }
    
    /// Template member function to print an error message
    /// @param message the error message to be printed
    /// @warning Giving a message type that cannot be printed to the output stream using << will result in an uncaught exception.
    template <typename T>
    void Error(T message) const {
        if (_level >= LevelError)
            *logger << "[ERROR]: " << message;
    }

    /// Template member function to print a warning message
    /// @param message the warning message to be printed
    /// @warning Giving a message type that cannot be printed to the output stream using << will result in an uncaught exception.
    template <typename T>
    void Warn(T message) const {
        if (_level >= LevelWarning)
            *logger << "[WARNING]: " << message;
    }

    /// Template member function to print an info message
    /// @param message the info message to be printed
    /// @warning Giving a message type that cannot be printed to the output stream using << will result in an uncaught exception.
    template <typename T>
    void Info(T message) const {
        if (_level >= LevelInfo)
            *logger << "[INFO]: " << message;
    }
};

#endif /* log_h */
