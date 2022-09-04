//
//  threshold_runner.cpp
//  FluxoniumGTC
//
//  Created by Jon on 7/4/22.
//

#include "XZZX_threshold_runner.hpp"
#include "XZZX_code.hpp"
#include <iostream>

using namespace std;
using namespace Eigen;

XZZX_threshold_runner::XZZX_threshold_runner(XZZX_threshold_params& params, std::mt19937& engine) :
    params(params),
    engine(engine)
{
    
}

void XZZX_threshold_runner::write_output(std::stringstream& strstream) {
    string str = strstream.str();
    if (params.filename != "") {
        file << str;
    }
    cout << str;
}

void XZZX_threshold_runner::run_instance(float p, int w, int h) {
    Vector4f idle_errors;
    idle_errors.tail(3) = params.idle_errors_normalized * p;
    idle_errors[0] = 1.0 - p;
    float faulty_meas = params.faulty_meas_factor * p;
    int n_rounds = 1;
    if (faulty_meas > 0) {
        n_rounds = max(w,h);
    }

    XZZX_params instance_params {
        idle_errors,
        faulty_meas,
        n_rounds,
        w,
        h
    };
    XZZX_code code(instance_params, engine);
    
    int n_success = 0;
    for (int i = 0; i < params.n_trials; i++) {
        n_success += code.run_decode();
    }
    
    stringstream str;
    str << p << "," << w << "," << h << "," << params.n_trials << "," << n_success << endl;
    write_output(str);
}

void XZZX_threshold_runner::run() {
    if (params.filename != "") {
        file.open(params.filename);
    }
    
    stringstream str;
    str << "p,w,h,n_trials,n_success" << endl;
    write_output(str);
    
    float p_step = (params.p_stop - params.p_start) / (params.p_steps - 1);
    for (int i = 0; i < params.p_steps; i++) {
        float p = params.p_start + p_step * i;
        
        for (int j = 0; j < params.sizes.size(); j++) {
            int w = params.sizes[j](0);
            int h = params.sizes[j](1);
            
            run_instance(p, w, h);
        }
    }
    
    if (params.filename != "") {
        file.close();
    }
}
