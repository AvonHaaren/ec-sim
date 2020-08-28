//
//  Timer.h
//  EvaporationSim
//
//  Created by Andreas von Haaren on 28/02/2020.
//  Copyright © 2020 Andreas von Haaren. All rights reserved.
//

#ifndef Timer_h
#define Timer_h

#include <chrono>
#include <iostream>

template<typename T>
struct Timer {
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    std::chrono::duration<float> duration;
    
    T m_name;
    
    Timer(T name = "Timer") : m_name(name) {
        start = std::chrono::high_resolution_clock::now();
        duration = std::chrono::seconds(0);
    }
    
    double currentSeconds() {
        duration = std::chrono::high_resolution_clock::now() - start;
        return duration.count();
    }
    
    ~Timer() {
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        
        float ms = duration.count() * 1000.f;
        
        std::cout << m_name << " took " << ms << " ms\n";
    }
};

#endif /* Timer_h */
