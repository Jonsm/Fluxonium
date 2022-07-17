//
//  XZZX_code.cpp
//  FluxoniumGTC
//
//  Created by Jon on 4/21/22.
//

//DONT USE THIS YET

#include "XZZX_GTC.hpp"
#include <iostream>
#include "math.h"

using namespace std;
using namespace Eigen;

XZZX_GTC::XZZX_GTC(XZZX_GTC_params params) :
params(params)
{
    init_qubits();
    init_stabs();
    init_graph();
    init_logicals();
}

//warning not efficient
bool XZZX_GTC::is_identified(int x1, int y1, int x2, int y2) {
    Matrix2f transform;
    transform.col(0) = params.l1.cast<float>();
    transform.col(1) = params.l2.cast<float>();
    Matrix2i transform_int = transform.cast<int>();
    Matrix2f transform_i = transform.inverse();
    
    Vector2f diff (x2-x1,y2-y1);
    Vector2f diff_tr = transform_i*diff;
    Vector2i diff_tr_round (rint(diff_tr.x()),rint(diff_tr.y()));
    Vector2i diff_closest = transform_int * diff_tr_round;
    
    return diff_closest.x() == x2 -x1 && diff_closest.y() == y2 - y1;
}

//not efficient
int XZZX_GTC::to_rep(int x, int y) {
    for (int x2 = 0; x2 < w; x2++) {
        for (int y2 = 0; y2 < h; y2++) {
            if (qubits_pos[x2][y2]==-1) {
                continue;
            }
            
            if (is_identified(x, y, x2, y2)) {
                return qubits_pos[x2][y2];
            }
        }
    }

    return -1;
}

void XZZX_GTC::init_qubits() {
    w = params.l1.x() + params.l2.x()+1;
    h = params.l1.y() + params.l2.y()+1;
    qubits_pos.resize(boost::extents[w][h]);
    fill(qubits_pos.data(), qubits_pos.data()+qubits_pos.num_elements(), -1);
    
    int x_max = 0;
    int y_max = 0;
    
    
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            bool added = false;
            for (int y2 = 0; y2 <= y; y2++) {
                int x_lim = w;
                if (y2 == y) {
                    x_lim = x;
                }
                for (int x2 = 0; x2 < x_lim; x2++) {
                    added = added || is_identified(x, y, x2, y2);
                }
            }
            
            if (!added) {
                qubits_pos[x][y] = n;
                n++;
                if (x+1 > x_max) {
                    x_max = x+1;
                }
                if (y+1 > y_max) {
                    y_max = y+1;
                }
            }
        }
    }
    
    w = x_max;
    h = y_max;
//
//    for (int y = 0; y < h; y++) {
//        for (int x = 0; x < w; x++) {
//            cout << qubits[x][y] << " ";
//        }
//        cout << endl;
//    }
}

void XZZX_GTC::init_stabs() {
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            if (qubits_pos[x][y]==-1) {
                continue;
            }
            
            int tl = qubits_pos[x][y];
            int tr = to_rep(x+1, y);
            int bl = to_rep(x, y+1);
            int br = to_rep(x+1, y+1);
            
            stabs.push_back(vector<int> {tl,tr,br,bl});
        }
    }
    
//    cout << endl;
//    for (int i = 0; i < stabs.size(); i++) {
//        for (int j = 0; j < 4; j++) {
//            cout << stabs[i][j] << ",";
//        }
//        cout << endl;
//    }
}

void XZZX_GTC::init_logicals() {
    Pauli paulis[] {Pauli::X, Pauli::Z};
    
    bool even_norm = !((params.l1.x() + params.l1.y()) % 2) &&
    !((params.l2.x() + params.l2.y()) % 2);
    
    for (Pauli p : paulis) {
        if (even_norm) {
            for (int ofs = 0; ofs <= 1; ofs++) {
                logicals.push_back(vector<Pauli>(n));
                fill(logicals.back().begin(), logicals.back().end(), Pauli::I);
                
                int curr = ofs, curr_x = ofs, curr_y = 0;
                do {
                    logicals.back()[curr] = p;
                    curr_y++;
                    if (p == Pauli::X) {
                        curr_x--;
                    } else {
                        curr_x++;
                    }
                    
                    curr = to_rep(curr_x, curr_y);
                } while (curr != ofs);
            }
        } else {
            logicals.push_back(vector<Pauli>(n));
            fill(logicals.back().begin(), logicals.back().end(), p);
        }
    }
    
    for (vector<Pauli>& logical: logicals) {
        for (int i = 0; i < n; i++) {
            cout << pauli_to_s(logical[i]) << " ";
        }
        cout << endl;
    }
}

void XZZX_GTC::init_graph() {
    boost::multi_array<pair<Pauli,int>, 2> adj_matrix(boost::extents[n][n]);
    fill(adj_matrix.data(), adj_matrix.data() + adj_matrix.num_elements(), pair(Pauli::NA,-1));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 4; j++) {
            Pauli p = Pauli::Z;
            if (j % 2 == 1) {
                p = Pauli::X;
            }
            
            int j_flip = (j+2) % 4;
            for (int k = 0; k < n; k++) {
                if (stabs[i][j] == stabs[k][j_flip]) {
                    adj_matrix[i][k] = pair(p,stabs[i][j]);
                    break;
                }
            }
        }
    }
    
    shortest_paths.resize(boost::extents[n][n]);
    
//    cout << endl;
//    for (int j = 0; j < n; j++) {
//        for (int i = 0; i < n; i++) {
//            char p_c = '-';
//            if (adj_matrix[i][j] == X) {
//                p_c = 'X';
//            } else if (adj_matrix[i][j] == Z) {
//                p_c = 'Z';
//            }
//            cout << p_c;
//        }
//        cout << endl;
//    }
}
