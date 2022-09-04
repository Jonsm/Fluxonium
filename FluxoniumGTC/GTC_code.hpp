//
//  GTC_code.hpp
//  FluxoniumGTC
//
//  Created by Jon on 7/23/22.
//

#ifndef GTC_code_hpp
#define GTC_code_hpp

#include <stdio.h>
#include <Eigen/Dense>
#include <random>
#include <vector>

struct GTC_params {
public:
    Eigen::Vector4f idle_errors;
    float faulty_meas;
    int n_rounds;
    Eigen::Vector2i l1;
    Eigen::Vector2i l2;
};

class GTC_code {
public:
    GTC_code(GTC_params& params, std::mt19937& engine);
    bool run_decode();
    void print_error();
    void print_syndrome(int round);
    void print_identification();
    void print_stabs();
    
private:
    std::mt19937& engine;
    std::uniform_real_distribution<float> dis;
    GTC_params params;
    
    int w;
    int h;
    int n;
    Eigen::MatrixXi qubits_pos;
    Eigen::MatrixXi stabs;
    Eigen::Vector2d xz_weight_MWPM;
    double t_weight_MWPM;
    std::vector<Eigen::VectorXi> logicals; //Z then X
    Eigen::MatrixXi connecting_qubits;
    
    Eigen::Vector4f cumulative_errors;
    Eigen::MatrixXd edge_weights;
    std::vector<Eigen::MatrixXi> logical_crossings;
    
    int round;
    int n_syndromes;
    Eigen::VectorXi errs_X;
    Eigen::VectorXi errs_Z;
    std::vector<Eigen::VectorXi> syndrome;
    Eigen::Vector4i logical_parities;
    
    std::vector<Eigen::Vector2i> verts_to_match;
    Eigen::VectorXi matches;
    
    bool is_identified(int x1, int y1, int x2, int y2);
    int to_rep(int x, int y);
    void init_qubits();
    void init_stabs();
    
    void init_logicals_2x();
    void init_4x_first2();
    void init_4x_second2();
    void init_logicals_4x();
    void init_logicals();
    
    void init_weights();
    void init_edges(std::vector<std::pair<int,int>>& edges, std::vector<double>& weights);
    void reconstruct_path(std::vector<int>& p, std::vector<double>& d, int src, int dest);
    void get_shortest_paths(std::vector<std::pair<int, int>> &edges, std::vector<double> &weights);
    void init_graph();
    
    void reset();
    void flip(int i, char type);
    void add_errs();
    int stab_parity(int stab);
    void meas_stabs();
    void single_round();
    
    void lemon_match();
    double stab_dist(Eigen::Vector2i& s1, Eigen::Vector2i& s2);
    void get_correction();
    bool check_correction();
};

#endif /* GTC_code_hpp */
