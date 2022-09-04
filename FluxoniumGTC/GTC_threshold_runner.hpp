//
//  GTC_threshold_runner.hpp
//  FluxoniumGTC
//
//  Created by Jon on 8/19/22.
//

#ifndef GTC_threshold_runner_hpp
#define GTC_threshold_runner_hpp

#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <Eigen/Dense>

struct GTC_threshold_params {
    std::string filename;
    int n_trials;
    float p_start;
    float p_stop;
    int p_steps;
    Eigen::Vector3f idle_errors_normalized;
    float faulty_meas_factor;
    std::vector<Eigen::Vector2i> l1s;
    std::vector<Eigen::Vector2i> l2s;
    std::vector<int> n_rounds;
};

class GTC_threshold_runner {
public:
    GTC_threshold_runner(GTC_threshold_params& params, std::mt19937& engine);
    void run();
    
private:
    std::mt19937& engine;
    GTC_threshold_params params;
    std::ofstream file;
    
    void run_instance(float p, int j);
    void write_output(std::stringstream& str);
};

#endif /* GTC_threshold_runner_hpp */
