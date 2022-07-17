//
//  XZZX_code.hpp
//  FluxoniumGTC
//
//  Created by Jon on 4/21/22.
//

//DONT USE THIS YET

#include <vector>
#include "boost/multi_array.hpp"
#include <Eigen/Dense>
#include "paulis.hpp"

#ifndef XZZX_GTC_hpp
#define XZZX_GTC_hpp

struct XZZX_GTC_params {
public:
    std::vector<float> idle_errors;
    std::vector<std::vector<float>> meas_induced_errors;
    float faulty_meas;
    float t_weight_MWPM;
    std::vector<float> xz_weight_MWPM;
    int n_rounds;
    Eigen::Vector2i l1;
    Eigen::Vector2i l2;
};

class XZZX_GTC {
public:
    XZZX_GTC(XZZX_GTC_params params);
    bool run_decode();
    
private:
    bool is_identified(int x1, int y1, int x2, int y2);
    int to_rep(int x, int y);
    void init_qubits();
    void init_stabs();
    void init_logicals();
    void init_graph();
    void idle_errors_round();
    void meas_errors_round();
    void single_round();
    void error_rounds();
    void decode();
    
    
    XZZX_GTC_params params;
    int w;
    int h;
    int n=0;
    boost::multi_array<int,2> qubits_pos;
    std::vector<std::vector<int>> stabs;
    std::vector<std::vector<Pauli>> logicals;
    boost::multi_array<int, 2> shortest_paths;
};

#endif /* XZZX_code_hpp */
