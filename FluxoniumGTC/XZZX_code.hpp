//
//  xzzx_code.hpp
//  FluxoniumGTC
//
//  Created by Jon on 6/6/22.
//

#ifndef xzzx_code_hpp
#define xzzx_code_hpp

#include <Eigen/Dense>
#include <vector>
#include <random>
#include "blossom5-v2.05.src/PerfectMatching.h"

struct XZZX_params {
    Eigen::Vector4f idle_errors;
    float faulty_meas;
    int n_rounds;
    int w;
    int h;
};

class XZZX_code {
public:
    XZZX_code(XZZX_params& params, std::mt19937& engine);
    bool run_decode();
    void print_error();
    void print_syndrome(int round);
    void print_correction();
    
private:
    std::mt19937& engine;
    std::uniform_real_distribution<float> dis;
    XZZX_params params;
    
    int n;
    Eigen::Vector4f cumulative_errors;
    Eigen::MatrixXi stabs;
    std::vector<Eigen::VectorXi> logicals_X;
    std::vector<Eigen::VectorXi> logicals_Z;
    Eigen::Vector2d xz_weight_MWPM;
    double t_weight_MWPM;
    
    int round;
    Eigen::Vector2i n_syndromes_v;
    Eigen::VectorXi errs_X;
    Eigen::VectorXi errs_Z;
    std::vector<Eigen::VectorXi> syndrome;
    Eigen::VectorXi correction_X;
    Eigen::VectorXi correction_Z;
    
    std::vector<Eigen::Vector2i> match_vertices;
    
    void init_stabs();
    void init_logicals();
    void init_weights();
    
    int stab_parity(int stab, bool with_correction);
    void reset();
    void add_errs();
    void meas_stabs();
    void single_round();
    
    bool check_infinity(int x_dist, int y_dist, int t_dist);
    double stab_dist(Eigen::Vector2i& s1, Eigen::Vector2i& s2);
    void get_matching(PerfectMatching& matching, int s_type);
    void make_string(int i1, int i2, int s_type);
    void make_strings(PerfectMatching& matching, int s_type);
    void get_correction();
    bool codespace_check();
    bool logical_check();
    bool check_correction();
    
    void print_helper(Eigen::VectorXi& X_type, Eigen::VectorXi& Z_type);
    std::string stab_to_str(int i);
    std::string qubit_to_str(int i);
};

#endif /* xzzx_code_hpp */
