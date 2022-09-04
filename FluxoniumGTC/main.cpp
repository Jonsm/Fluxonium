//
//  main.cpp
//  FluxoniumGTC
//
//  Created by Jon on 4/19/22.
//

#include <iostream>
#include "GTC_threshold_runner.hpp"
#include <Eigen/Dense>
#include <random>
#include <chrono>

using namespace std;


int main(int argc, const char * argv[]) {
    ///////////////////////////////PARAMS
    string filename = ""; //filename to output, "" for no output file
    if (argc > 1) {
        filename = argv[1]; //input filename as command line argument
    }
    
    int n_trials = 2000; //number of trials
    float p_start = .1; //start value of total idle error p=p_X+p_Y+p_Z
    float p_stop = .3; //stop value
    int p_steps = 11; //number steps, including start and stop
    Eigen::Vector3f idle_errors_normalized {0, 0, 1}; //normalized p_X,p_Y,p_Z
    float faulty_meas_factor = 1.0; //probability of faulty measurement = faulty_meas_factor*p
    
    //periodicity of the lattice is given by a set of two vectors, l1, l2.
    std::vector<Eigen::Vector2i> l1s {Eigen::Vector2i {7,0},Eigen::Vector2i {11,0},Eigen::Vector2i {13,0}};
    std::vector<Eigen::Vector2i> l2s {Eigen::Vector2i {0,11},Eigen::Vector2i {0,13},Eigen::Vector2i {0,17}};
    std::vector<int> n_rounds {10,10,10}; //number of rounds to do error correction. One number for each set of l1,l2. *Might change later.
    ///////////////////////////////
    
    GTC_threshold_params params {
        filename,
        n_trials,
        p_start,
        p_stop,
        p_steps,
        idle_errors_normalized,
        faulty_meas_factor,
        l1s,
        l2s,
        n_rounds
    };
    mt19937 engine;
    engine.seed(std::random_device{}());
    GTC_threshold_runner runner(params, engine);
    
    auto t1 = chrono::high_resolution_clock::now();
    
    runner.run();

    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> dur = std::chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << dur.count() << endl;
    
    return 0;
}
