//
//  main.cpp
//  FluxoniumGTC
//
//  Created by Jon on 4/19/22.
//

#include <iostream>
#include "threshold_runner.hpp"
#include <Eigen/Dense>
#include <random>
#include <chrono>

using namespace std;
using namespace Eigen;

int main(int argc, const char * argv[]) {
    ///////////////////////////////PARAMS
    string filename = ""; //filename to output, "" for no output file
    if (argc > 1) {
        filename = argv[1]; //input filename as command line argument
    }
    
    int n_trials = 1000; //number of trials
    float p_start = .04; //start value of total idle error p=p_X+p_Y+p_Z
    float p_stop = .07; //stop value
    int p_steps = 11; //number steps, including start and stop
    Eigen::Vector3f idle_errors_normalized {0.5, 0, 0.5}; //normalized p_X+p_Y+p_Z
    float faulty_meas_factor = 0.5; //probability of faulty measurement = faulty_meas_factor*p
    
    std::vector<Eigen::Vector2i> sizes {Eigen::Vector2i {8,8},Eigen::Vector2i {16,16},Eigen::Vector2i {24,24}}; //system sizes, as vector of w,h
    ///////////////////////////////
    
    XZZX_threshold_params params {
        filename,
        n_trials,
        p_start,
        p_stop,
        p_steps,
        idle_errors_normalized,
        faulty_meas_factor,
        sizes
    };
    mt19937 engine;
    engine.seed(std::random_device{}());
    XZZX_threshold_runner runner(params, engine);
    
    auto t1 = chrono::high_resolution_clock::now();
    
    runner.run();

    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> dur = std::chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << dur.count() << endl;
    
    return 0;
}
