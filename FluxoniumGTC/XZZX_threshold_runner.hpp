//
//  threshold_runner.hpp
//  FluxoniumGTC
//
//  Created by Jon on 7/4/22.
//

#ifndef threshold_runner_hpp
#define threshold_runner_hpp

#include <fstream>
#include <sstream>
#include <random>
#include <vector>
#include <Eigen/Dense>

struct XZZX_threshold_params {
    std::string filename;
    int n_trials;
    float p_start;
    float p_stop;
    int p_steps;
    Eigen::Vector3f idle_errors_normalized;
    float faulty_meas_factor;
    std::vector<Eigen::Vector2i> sizes;
};

class XZZX_threshold_runner {
public:
    XZZX_threshold_runner(XZZX_threshold_params& params, std::mt19937& engine);
    void run();
    
private:
    std::mt19937& engine;
    XZZX_threshold_params params;
    std::ofstream file;
    
    void run_instance(float p, int w, int h);
    void write_output(std::stringstream& str);
};

#endif /* threshold_runner_hpp */
