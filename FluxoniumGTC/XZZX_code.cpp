//
//  xzzx_code.cpp
//  FluxoniumGTC
//
//  Created by Jon on 6/6/22.
//

#include "XZZX_code.hpp"
#include <iostream>
#include "math.h"
#include <sstream>
#include <limits>

using namespace Eigen;
using namespace std;

void XZZX_code::init_stabs() {
    int np = params.w*params.h;
    for (int x = 0; x < params.w; x++) {
        for (int y = 0; y < params.h; y++) {
            int stab_ind = 2*(y*params.w+x);
            
            stabs(stab_ind,0) = stab_ind/2;
            stabs(stab_ind,1) = stab_ind/2 + np;
            stabs(stab_ind,2) = y*params.w + (x+1)%params.w;
            stabs(stab_ind,3) = ((y+1)%params.h) * params.w + x + np;
            
            stabs(stab_ind+1,0) = stab_ind/2 + np;
            stabs(stab_ind+1,1) = ((y-1+params.h) % params.h) * params.w + (x+1) % params.w;
            stabs(stab_ind+1,2) = y*params.w + (x+1) % params.w + np;
            stabs(stab_ind+1,3) = y*params.w + (x+1) % params.w;
        }
    }
}

void XZZX_code::init_logicals() {
    for (int i = 0; i < 2; i++) {
        logicals_X[i].resize(params.h);
        logicals_Z[i].resize(params.w);
        
        int x = 0;
        for (int y = 0; y < params.h; y++) {
            logicals_X[i](y) = params.w*y+x + i*params.h*params.w;
        }
        int y = 0;
        for (int x =0; x < params.w; x++) {
            logicals_Z[i](x) = params.w*y+x + i*params.h*params.w;
        }
    }
}

void XZZX_code::init_weights() {
    double px = params.idle_errors(1) + params.idle_errors(2);
    double pz = params.idle_errors(2) + params.idle_errors(3);
    double p_id = params.idle_errors(0);
    double pt = params.faulty_meas;
    xz_weight_MWPM(0) = log(px/p_id)*-1.0;
    xz_weight_MWPM(1) = log(pz/p_id)*-1.0;
    t_weight_MWPM = log(pt/(1-pt))*-1.0;
    
    if (px == 0) {
        xz_weight_MWPM(0) = 0;
    }
    if (pz == 0) {
        xz_weight_MWPM(1) = 0;
    }
    if (pt == 0) {
        t_weight_MWPM = 0;
    }
}

XZZX_code::XZZX_code(XZZX_params& params, std::mt19937& engine) :
params(params),
stabs(params.w*params.h*2,4),
errs_X(VectorXi::Constant(params.w*params.h*2, 0)),
errs_Z(VectorXi::Constant(params.w*params.h*2, 0)),
correction_X(VectorXi::Constant(params.w*params.h*2, 0)),
correction_Z(VectorXi::Constant(params.w*params.h*2, 0)),
syndrome(params.n_rounds),
logicals_X(2),
logicals_Z(2),
engine(engine),
dis(0.0,1.0),
n(2*params.w*params.h),
match_vertices(2*params.w*params.h*params.n_rounds)
{
    for (int i = 0; i < params.n_rounds; i++) {
        syndrome[i] = VectorXi::Constant(n, 0);
    }
    cumulative_errors[0] = params.idle_errors[0];
    for (int i = 1; i < 4; i++) {
        cumulative_errors[i] = cumulative_errors[i-1] + params.idle_errors[i];
    }
    
    init_stabs();
    init_logicals();
    init_weights();
}

int XZZX_code::stab_parity(int stab, bool with_correction) {
    int parity = 0;
    for (int j = 0; j < 4; j++) {
        if (j % 2 == 0) {
            parity ^= errs_Z(stabs(stab,j));
            parity ^= (correction_Z(stabs(stab,j)) & with_correction);
        } else {
            parity ^= errs_X(stabs(stab,j));
            parity ^= (correction_X(stabs(stab,j)) & with_correction);
        }
    }
    return parity;
}

void XZZX_code::reset() {
    round = 0;
    n_syndromes_v = {0,0};
    errs_X = VectorXi::Constant(n, 0);
    errs_Z = VectorXi::Constant(n, 0);
    for (int i = 0; i < params.n_rounds; i++) {
        syndrome[i] = VectorXi::Constant(n, 0);
    }
    correction_X = VectorXi::Constant(n, 0);
    correction_Z = VectorXi::Constant(n, 0);
}

void XZZX_code::add_errs() {
    for (int i = 0; i < n; i++) {
        float p = dis(engine);
        if (p > cumulative_errors[0] && p <= cumulative_errors[1]) {
            errs_X(i) ^= 1;
        } else if (p > cumulative_errors[1] && p <= cumulative_errors[2]) {
            errs_X(i) ^= 1;
            errs_Z(i) ^= 1;
        } else if (p > cumulative_errors[2]) {
            errs_Z(i) ^= 1;
        }
    }
}

//do more stuff here, if we want a more complicated error model
void XZZX_code::meas_stabs() {
    for (int i = 0; i < n; i++) {
        int parity = stab_parity(i, false);
        
        if (round != params.n_rounds - 1) {
            float p = dis(engine);
            if (p < params.faulty_meas) {
                parity ^= 1;
            }
        }
        
        syndrome[round](i) = parity;
        
        if (round == 0) {
            n_syndromes_v[i%2] += parity;
        } else {
            n_syndromes_v[i%2] += (parity ^ syndrome[round-1](i));
        }
    }
}

void XZZX_code::single_round() {
    add_errs();
    meas_stabs();
    round++;
}

bool XZZX_code::check_infinity(int x_dist, int y_dist, int t_dist) {
    double px = params.idle_errors(1) + params.idle_errors(2);
    double pz = params.idle_errors(2) + params.idle_errors(3);
    double pt = params.faulty_meas;
    
    if (x_dist != 0 && pz == 0) {
        return true;
    } else if (y_dist != 0 && px == 0) {
        return true;
    } else if (t_dist != 0 && pt == 0) {
        return true;
    }
    return false;
}

double XZZX_code::stab_dist(Eigen::Vector2i& s1, Eigen::Vector2i& s2) {
    int stab_loc1 = s1.x()/2;
    int stab_loc2 = s2.x()/2;
    int t1 = s1.y();
    int t2 = s2.y();
    int x1 = stab_loc1 % params.w;
    int y1 = stab_loc1 / params.w;
    int x2 = stab_loc2 % params.w;
    int y2 = stab_loc2 / params.w;
    
    int x_dist = min(abs(x1-x2),params.w-abs(x1-x2));
    int y_dist = min(abs(y1-y2),params.h-abs(y1-y2));
    int t_dist = abs(t1-t2);
    
    if (check_infinity(x_dist, y_dist, t_dist)) {
        return -1;
    }
    
    return xz_weight_MWPM(0)*y_dist + xz_weight_MWPM(1)*x_dist + t_weight_MWPM*t_dist;
}

void XZZX_code::get_matching(PerfectMatching& matching, int s_type) {
    int ct = 0;
    for (int r = 0; r < params.n_rounds; r++) {
        for (int i = s_type; i < n; i+=2) {
            if (r == 0) {
                if (syndrome[r](i)) {
                    match_vertices[ct] = {i,r};
                    ct++;
                }
            } else {
                if (syndrome[r](i) != syndrome[r-1](i)) {
                    match_vertices[ct] = {i,r};
                    ct++;
                }
            }
        }
    }
    
    int n_syndromes = n_syndromes_v[s_type];
    for (int i = 0; i < n_syndromes; i++) {
        for (int j = i+1; j < n_syndromes; j++) {
            double dist = stab_dist(match_vertices[i], match_vertices[j]);
            if (dist != -1) {
                matching.AddEdge(i, j, dist);
            }
        }
    }
    
    matching.options.verbose = false;
    matching.Solve();
}

void XZZX_code::make_string(int i1, int i2, int s_type) {
    int stab_loc1 = match_vertices[i1].x()/2;
    int stab_loc2 = match_vertices[i2].x()/2;
    int x1 = stab_loc1 % params.w;
    int y1 = stab_loc1 / params.w;
    int x2 = stab_loc2 % params.w;
    int y2 = stab_loc2 / params.w;
    
    int x_dir = 1;
    int y_dir = 1;
    if (x2 - x1 > params.w / 2 || (x1 - x2 >=0 && x1 - x2 < params.w / 2)) {
        x_dir = 0; //going left
    }
    if (y2 - y1 > params.h / 2 || (y1 - y2 >=0 && y1 - y2 < params.h / 2)) {
        y_dir = 0; //going up
    }
    
    int x = x1;
    while (x != x2) {
        int Z_ind = y1*params.w + (x + x_dir) % params.w;
        if (s_type == 1) {
            Z_ind += params.h*params.w;
        }
        correction_Z[Z_ind] ^= 1;
        x = (x-1+2*x_dir+params.w) % params.w;
    }
    
    int y = y1;
    while (y != y2) {
        int X_ind = ((y+y_dir) % params.h)*params.w + x2 + params.h*params.w;
        if (s_type == 1) {
            X_ind = ((y-1+y_dir + params.h) % params.h)*params.w + (x2+1) % params.w;
        }
        correction_X[X_ind] ^= 1;
        y = (y-1+2*y_dir+params.h) % params.h;
    }
}

void XZZX_code::make_strings(PerfectMatching& matching, int s_type) {
    int n_syndromes = n_syndromes_v[s_type];
    for (int i = 0; i < n_syndromes; i++) {
        if (match_vertices[i].x() == -1) {
            continue;
        }
        
        int match = matching.GetMatch(i);
        make_string(i, match, s_type);
        match_vertices[match].x() = -1;
    }
}

void XZZX_code::get_correction() {
    for (int s_type = 0; s_type <= 1; s_type++) {
        int n_syndromes = n_syndromes_v[s_type];
        int n_edges = (n_syndromes * (n_syndromes - 1)) / 2;
        PerfectMatching matching(n_syndromes, n_edges);
    
        get_matching(matching,s_type);
        make_strings(matching,s_type);
    }
}

bool XZZX_code::codespace_check() {
    for (int i = 0; i < n; i++) {
        if (stab_parity(i,true)) {
            cout << "VIOLATED STABILIZER AT " << i << endl;
            return false;
        }
    }
    
    return true;
}

bool XZZX_code::logical_check() {
    for (int i = 0; i < 2; i++) {
        int parity = 0;
        for (int j = 0; j < logicals_X[i].size(); j++) {
            int ind = logicals_X[i](j);
            parity ^= (errs_Z[ind] ^ correction_Z[ind]);
        }
        
        if (parity) {
            return false;
        }
    }
    
    for (int i = 0; i < 2; i++) {
        int parity = 0;
        for (int j = 0; j < logicals_Z[i].size(); j++) {
            int ind = logicals_Z[i](j);
            parity ^= (errs_X[ind] ^ correction_X[ind]);
        }
        
        if (parity) {
            return false;
        }
    }
    
    return true;
}

bool XZZX_code::check_correction() {
    return codespace_check() && logical_check();
}

bool XZZX_code::run_decode() {
    reset();
    
    for (int i = 0; i < params.n_rounds; i++) {
        single_round();
    }

    
    get_correction();
    return check_correction();
}

void XZZX_code::print_helper(Eigen::VectorXi& X_type, Eigen::VectorXi& Z_type) {
    for (int y = 0; y < params.h; y++) {
        for (int x = 0; x < params.w; x++) {
            int ind_horiz = params.w*y+x+params.w*params.h;
            char err_type = '-';
            if (X_type(ind_horiz) && Z_type(ind_horiz)) {
                err_type = 'Y';
            } else if (X_type(ind_horiz)) {
                err_type = 'X';
            } else if (Z_type(ind_horiz)) {
                err_type = 'Z';
            }
            cout << "+" << err_type;
        }
        cout << "-" << endl;
        
        for (int x = 0; x < params.w; x++) {
            int ind_vert = params.w*y+x;
            char err_type = '|';
            if (X_type(ind_vert) && Z_type(ind_vert)) {
                err_type = 'Y';
            } else if (X_type(ind_vert)) {
                err_type = 'X';
            } else if (Z_type(ind_vert)) {
                err_type = 'Z';
            }
            cout << err_type << " ";
        }
        cout << endl;
    }
    
    for (int x = 0; x < params.w; x++) {
        cout << "| ";
    }
    cout << endl;
}

void XZZX_code::print_error() {
    print_helper(errs_X, errs_Z);
}

void XZZX_code::print_syndrome(int round) {
    for (int y = 0; y < params.h; y++) {
        char syndrome_char = ' ';
        int last_right =2*(params.w*y + params.w-1)+1;
        if (syndrome[round](last_right)) {
            syndrome_char = 'x';
        }
        cout << syndrome_char;
        
        for (int x = 0; x < params.w-1; x++) {
            int ind_v = 2*(y*params.w+x)+1;
            char syndrome_char = ' ';
            if (syndrome[round](ind_v)) {
                syndrome_char = 'x';
            }
            cout << "-" << syndrome_char;
        }
        cout << "-" << endl;
        
        for (int x = 0; x < params.w; x++) {
            int ind_p = 2*(y*params.w+x);
            char syndrome_char = ' ';
            if (syndrome[round](ind_p)) {
                syndrome_char = 'x';
            }
            cout << "|" << syndrome_char;
        }
        cout << endl;
    }
}

void XZZX_code::print_correction() {
    print_helper(correction_X, correction_Z);
}

string XZZX_code::stab_to_str(int i) {
    int x = (i / 2) % params.w;
    int y = (i / 2) / params.w;
    char type = 'P';
    if (i % 2 == 1) {
        type = 'V';
    }
    stringstream str;
    str << type << "(" << x << "," << y << ")";
    return str.str();
}

string XZZX_code::qubit_to_str(int i) {
    char dir = 'y';
    if (i >= params.w*params.h) {
        i -= params.w*params.h;
        dir = 'x';
    }
    int x = i % params.w;
    int y = i / params.w;
    stringstream str;
    str << dir << "(" << x << "," << y << ")";
    return str.str();
}
