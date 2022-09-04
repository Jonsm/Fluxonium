//
//  GTC_threshold_runner.cpp
//  FluxoniumGTC
//
//  Created by Jon on 8/19/22.
//

#include "GTC_threshold_runner.hpp"
#include "GTC_code.hpp"
#include <iostream>

using namespace std;
using namespace Eigen;

GTC_threshold_runner::GTC_threshold_runner(GTC_threshold_params& params, std::mt19937& engine) :
    params(params),
    engine(engine)
{
    
}

void GTC_threshold_runner::write_output(std::stringstream& strstream) {
    string str = strstream.str();
    if (params.filename != "") {
        file << str;
    }
    cout << str;
}

void GTC_threshold_runner::run_instance(float p, int j) {
    Vector4f idle_errors;
    idle_errors.tail(3) = params.idle_errors_normalized * p;
    idle_errors(0) = 1.0 - p;
    float faulty_meas = params.faulty_meas_factor * p;
    int n_rounds = params.n_rounds[j];
    Vector2i l1 = params.l1s[j];
    Vector2i l2 = params.l2s[j];
    
    GTC_params instance_params {
        idle_errors,
        faulty_meas,
        n_rounds,
        l1,
        l2
    };
    GTC_code code(instance_params, engine);
    
    int n_success = 0;
    for (int i = 0; i < params.n_trials; i++) {
        n_success += code.run_decode();
    }
    
    stringstream str;
    str << p << "," << l1.x() << "," << l1.y() << "," << l2.x() << "," << l2.y() << "," << params.n_trials << "," << n_success << endl;
    write_output(str);
}

void GTC_threshold_runner::run() {
    if (params.filename != "") {
        file.open(params.filename);
    }
    
    stringstream str;
    str << "p,l1_x,l1_y,l2_x,l2_y,n_trials,n_success" << endl;
    write_output(str);
    
    float p_step = (params.p_stop - params.p_start) / (params.p_steps - 1);
    for (int i = 0; i < params.p_steps; i++) {
        float p = params.p_start + p_step * i;
        
        for (int j = 0; j < params.l1s.size(); j++) {
            run_instance(p, j);
        }
    }
    
    if (params.filename != "") {
        file.close();
    }
}
