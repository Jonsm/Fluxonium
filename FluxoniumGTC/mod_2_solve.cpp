//
//  mod_2_solve.cpp
//  FluxoniumGTC
//
//  Created by Jon on 5/16/22.
//

#include "mod_2_solve.hpp"

using namespace Eigen;
using namespace std;

VectorXi mod_2_solve(VectorXi L, MatrixXi M) {
    int l = (int)L.size();
    int h = (int)M.rows();
    int n_pivots = 0;
    VectorXi soln(h);
    vector<int> pivots(l);
    
    for (int col = 0; col < l; col++) {
        for (int row = 0; row < h; row++) {
            if (M(row,col) == 1 && pivots[row] != 1) {
                pivots[row] = 1;
                n_pivots++;
                
                for (int row2 = 0; row < l; row2++) {
                    if (M(row2,col) == 1) {
                        M.row(row2) += M.row(row); //mod 2
                        L[row2] += L[row];
                    }
                }
                break;
            }
        }
        
        if (n_pivots == h) {
            break;
        }
    }
    
    return soln;
}
