//
//  mod_2_solve.hpp
//  FluxoniumGTC
//
//  Created by Jon on 5/16/22.
//

#ifndef mod_2_solve_hpp
#define mod_2_solve_hpp

#include <Eigen/Dense>
#include <vector>

//return solution of M*x = L
Eigen::VectorXi mod_2_solve(Eigen::VectorXi L, Eigen::MatrixXi M);

#endif /* mod_2_solve_hpp */
