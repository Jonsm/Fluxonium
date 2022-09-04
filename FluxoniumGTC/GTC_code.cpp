//
//  GTC_code.cpp
//  FluxoniumGTC
//
//  Created by Jon on 7/23/22.
//

#include "GTC_code.hpp"
#include <iostream>
#include <vector>
#include <limits>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include "blossom5-v2.05.src/PerfectMatching.h"
extern "C" {
    #include "spasm.src/spasm.h"
}
#include <lemon/matching.h>
#include <lemon/list_graph.h>

using namespace std;
using namespace Eigen;
using namespace boost;

bool GTC_code::is_identified(int x1, int y1, int x2, int y2) {
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

int GTC_code::to_rep(int x, int y) {
    for (int x2 = 0; x2 < w; x2++) {
        for (int y2 = 0; y2 < h; y2++) {
            if (qubits_pos(x2,y2)==-1) {
                continue;
            }
            
            if (is_identified(x, y, x2, y2)) {
                return qubits_pos(x2,y2);
            }
        }
    }

    return -1;
}

void GTC_code::init_qubits() {
    w = params.l1.x() + params.l2.x()+1;
    h = params.l1.y() + params.l2.y()+1;
    n=0;
    qubits_pos = MatrixXi::Constant(w, h, -1);
    assert (!is_identified(0, 0, 2, 0));
    
    for (int y = 0; y <h; y++) {
        for (int x = 0; x < w; x++) {
            if (to_rep(x,y) == -1) {
                qubits_pos(x,y)=n;
                n++;
            }
        }
    }
}

void GTC_code::init_stabs() {
    stabs.resize(n, 4);
    int stab_ct = 0;
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            if (qubits_pos(x,y)==-1) {
                continue;
            }
            
            stabs(stab_ct, 0) = qubits_pos(x,y); //tl
            stabs(stab_ct, 1) = to_rep(x+1, y); //tr
            stabs(stab_ct, 2) = to_rep(x+1, y+1); //br
            stabs(stab_ct, 3) = to_rep(x, y+1); //bl
            stab_ct++;
        }
    }
}

void GTC_code::init_weights() {
    double px = params.idle_errors(1) + params.idle_errors(2);
    double pz = params.idle_errors(2) + params.idle_errors(3);
    double p_id = params.idle_errors(0);
    double pt = params.faulty_meas;
    xz_weight_MWPM(0) = log(px/p_id)*-1.0;
    xz_weight_MWPM(1) = log(pz/p_id)*-1.0;
    t_weight_MWPM = log(pt/(1-pt))*-1.0;
    
    if (px == 0) {
        xz_weight_MWPM(0) = -1;
    }
    if (pz == 0) {
        xz_weight_MWPM(1) = -1;
    }
    if (pt == 0) {
        t_weight_MWPM = -1;
    }
}

void GTC_code::init_logicals_2x() {
    logicals.push_back(VectorXi::Constant(2*n,0));
    logicals.back().head(n) = VectorXi::Constant(n,1);
    logicals.push_back(VectorXi::Constant(2*n,0));
    logicals.back().tail(n) = VectorXi::Constant(n,1);
}

void GTC_code::init_4x_first2() {
    logicals.push_back(Eigen::VectorXi::Constant(2*n, 0));
    logicals.push_back(Eigen::VectorXi::Constant(2*n, 0));
    
    for (int i = 0; i < 2; i++) {
        int x = i;
        int y = 0;
        int current = to_rep(x,y);
        int start = current;
        
        do {
            logicals[i](current) = 1;
            x--;
            y++;
            current = to_rep(x, y);
        } while (current != start);
    }
}

void GTC_code::init_4x_second2() {
    vector<int> sparse_cols;
    vector<int> sparse_rows;
    vector<int> sparse_sol;
    
    for (int i = 0; i < stabs.rows(); i++) {
        sparse_cols.push_back(i+1);
        sparse_rows.push_back(stabs(i,0));//X
        sparse_cols.push_back(i+1);
        sparse_rows.push_back(stabs(i,2));//X
        sparse_cols.push_back(i+1);
        sparse_rows.push_back(n+stabs(i,1));//Z
        sparse_cols.push_back(i+1);
        sparse_rows.push_back(n+stabs(i,3));//Z
        sparse_sol.push_back(0);
    }
    
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < logicals[j].size(); i++) {
            if (logicals[j](i)) {
                sparse_cols.push_back(n+j);
                sparse_rows.push_back(i);
            }
        }
    }
    sparse_sol.push_back(1);
    sparse_sol.push_back(0);
    
    int nzmax = (int)sparse_rows.size();
    int nz = nzmax;
    int rows = 2*n;
    int cols = n+2;
    int prime = 2;
    vector<int> entries(sparse_rows.size());
    fill(entries.begin(), entries.end(), 1);
    spasm_triplet T {nzmax,nz,rows,cols,sparse_rows.data(),sparse_cols.data(),entries.data(),prime};
    spasm* A = spasm_compress(&T);
    
    VectorXi swapped(2*n);
    spasm_LU_solve(A, sparse_sol.data(), swapped.data());
    logicals.push_back(VectorXi(2*n));
    logicals.back().head(n) = swapped.tail(n);
    logicals.back().tail(n) = swapped.head(n);
    
    sparse_sol[n] = 0;
    sparse_sol[n+1] = 1;
    spasm_LU_solve(A, sparse_sol.data(), swapped.data());
    logicals.push_back(VectorXi(2*n));
    logicals.back().head(n) = swapped.tail(n);
    logicals.back().tail(n) = swapped.head(n);
}

void GTC_code::init_logicals_4x() {
    init_4x_first2();
    init_4x_second2();
}

void GTC_code::init_logicals() {
    if ((params.l1.x()+params.l1.y()) % 2 == 0 && (params.l2.x() % 2 == 0+params.l2.y()) % 2 == 0) {
        init_logicals_4x();
    } else {
        init_logicals_2x();
    }
}

void GTC_code::init_edges(std::vector<std::pair<int, int>> &edges, std::vector<double> &weights) {
    int tl_stab_i = -1;
    for (int i = 0; i < n; i++) {
        if (xz_weight_MWPM(1) != -1) {
            for (int j = 0; j < n; j++) {
                if (stabs(j,2)==stabs(i,0)) {
                    tl_stab_i = j;
                }
            }
            
            if (connecting_qubits(i,tl_stab_i) == -1) {
                edges.push_back(pair<int,int>{i, tl_stab_i});
                weights.push_back(xz_weight_MWPM(1));
                connecting_qubits(i,tl_stab_i) = stabs(i,0);
                connecting_qubits(tl_stab_i,i) = stabs(i,0);
            }
        }
        
        int tr_stab_i = -1;
        if (xz_weight_MWPM(0) != -1) {
            for (int j = 0; j < n; j++) {
                if (stabs(j,3) == stabs(i,1)) {
                    tr_stab_i = j;
                }
            }
            
            if (connecting_qubits(i,tr_stab_i)) {
                edges.push_back(pair<int,int>{i, tr_stab_i});
                weights.push_back(xz_weight_MWPM(0));
                connecting_qubits(i,tr_stab_i) = stabs(i,1);
                connecting_qubits(tr_stab_i,i) = stabs(i,1);
            }
        }
    }
}

void GTC_code::reconstruct_path(std::vector<int> &p, std::vector<double> &d, int src, int dest) {
    if (d[dest] == numeric_limits<double>::max()) {
        edge_weights(src,dest)=-1;
        edge_weights(dest,src)=-1;
        return;
    }
    
    edge_weights(src,dest)=d[dest];
    edge_weights(dest,src)=d[dest];
    
    int curr = dest;
    while (curr != src) {
        int next = p[curr];
        int connecting_qubit = connecting_qubits(curr,next);
        if (connecting_qubit == stabs(curr,1) || connecting_qubit == stabs(curr,3)) {
            connecting_qubit += n;
        }
        
        for (int i = 0; i < logicals.size(); i++) {
            (logical_crossings[i])(src,dest) ^= (logicals[i])(connecting_qubit);
            (logical_crossings[i])(dest,src) ^= (logicals[i])(connecting_qubit);
        }
        
        curr = next;
    }
}

void GTC_code::get_shortest_paths(std::vector<std::pair<int, int>> &edges, std::vector<double> &weights) {
    typedef adjacency_list < listS, vecS, undirectedS,
      no_property, property < edge_weight_t, double > > graph_t;
    graph_t graph(edges.begin(), edges.end(), weights.begin(), n);
    
    vector<int> p(n);
    vector<double> d(n);
    for (int src = 0; src < n; src++) {
        int s = (int)vertex(src, graph);
        dijkstra_shortest_paths(graph, s, predecessor_map(&p[0]).distance_map(&d[0]));
        
        for (int dest = src; dest < n; dest++) {
            reconstruct_path(p, d, src, dest);
        }
    }
}

void GTC_code::init_graph() {
    vector<pair<int,int>> edges;
    vector<double> weights;
    connecting_qubits = MatrixXi::Constant(n, n, -1);
    init_edges(edges, weights);

    edge_weights = MatrixXd::Constant(n,n,0);
    
    for (int i = 0; i < logicals.size(); i++) {
        logical_crossings.push_back(MatrixXi::Constant(n, n, 0));
    }
    get_shortest_paths(edges, weights);
}

GTC_code::GTC_code(GTC_params& params, std::mt19937& engine) :
params(params),
engine(engine),
dis(0.0,1.0),
syndrome(params.n_rounds)
{
    init_qubits();
    init_stabs();
    init_logicals();
    init_weights();
    init_graph();
    
    for (int i = 0; i < params.n_rounds; i++) {
        syndrome[i] = VectorXi::Constant(n, 0);
    }
    cumulative_errors[0] = params.idle_errors[0];
    for (int i = 1; i < 4; i++) {
        cumulative_errors[i] = cumulative_errors[i-1] + params.idle_errors[i];
    }
    verts_to_match.resize(n*params.n_rounds);
    matches.resize(n*params.n_rounds);
}

void GTC_code::reset() {
    round = 0;
    n_syndromes = 0;
    errs_X = VectorXi::Constant(n, 0);
    errs_Z = VectorXi::Constant(n, 0);
    for (int i = 0; i < params.n_rounds; i++) {
        syndrome[i] = VectorXi::Constant(n, 0);
    }
    logical_parities = Vector4i::Constant(0);
}

void GTC_code::flip(int i, char type) {
    if (type == 'X') {
        errs_X(i) ^= 1;
        for (int j = 0; j < logicals.size(); j++) {
            logical_parities(j) ^= logicals[j](i+n);
        }
    } else {
        errs_Z(i) ^= 1;
        for (int j = 0; j < logicals.size(); j++) {
            logical_parities(j) ^= logicals[j](i);
        }
    }
}

void GTC_code::add_errs() {
    for (int i = 0; i < n; i++) {
        float p = dis(engine);
        if (p > cumulative_errors[0] && p <= cumulative_errors[1]) {
            flip(i, 'X');
        } else if (p > cumulative_errors[1] && p <= cumulative_errors[2]) {
            flip(i, 'X');
            flip(i, 'Z');
        } else if (p > cumulative_errors[2]) {
            flip(i, 'Z');
        }
    }
}

int GTC_code::stab_parity(int stab) {
    int parity = 0;
    for (int j = 0; j < 4; j++) {
        if (j % 2 == 0) {
            parity ^= errs_Z(stabs(stab,j));
        } else {
            parity ^= errs_X(stabs(stab,j));
        }
    }
    return parity;
}

//do more stuff here, if we want a more complicated error model
void GTC_code::meas_stabs() {
    for (int i = 0; i < n; i++) {
        int parity = stab_parity(i);
        
        if (round != params.n_rounds - 1) {
            float p = dis(engine);
            if (p < params.faulty_meas) {
                parity ^= 1;
            }
        }
        
        syndrome[round](i) = parity;
        
        if (round == 0) {
            n_syndromes += parity;
        } else {
            n_syndromes += (parity ^ syndrome[round-1](i));
        }
    }
}

void GTC_code::single_round() {
    add_errs();
    meas_stabs();
    round++;
}

double GTC_code::stab_dist(Eigen::Vector2i& s1, Eigen::Vector2i& s2) {
    int t1 = s1.y();
    int t2 = s2.y();
    int x1 = s1.x();
    int x2 = s2.x();
    double d = edge_weights(x1,x2);
    
    if (d == -1) {
        return -1;
    }
    
    int t_dist = abs(t1-t2);
    if (t_weight_MWPM == -1 && t_dist != 0) {
        return -1;
    }
    
    return d + t_weight_MWPM*t_dist;
}

//alternate graph matching library, unused because slower
void GTC_code::lemon_match() {
    lemon::ListGraph graph;
    lemon::ListGraph::EdgeMap<double> weight(graph);
    
    for (int i = 0; i < n_syndromes; i++) {
        graph.addNode();
    }
    
    for (lemon::ListGraph::NodeIt i(graph); i!=lemon::INVALID; ++i) {
        for (lemon::ListGraph::NodeIt j = i; j!=lemon::INVALID; ++j) {
            if (j == i) {
                continue;
            }
            int i_ind = graph.id(i);
            int j_ind = graph.id(j);
            double dist = stab_dist(verts_to_match[i_ind], verts_to_match[j_ind]);
            
            if (dist != -1) {
                lemon::ListGraph::Edge edge = graph.addEdge(i,j);
                weight.set(edge, dist*-1);
            }
        }
    }
    
    lemon::MaxWeightedPerfectMatching<lemon::ListGraph, lemon::ListGraph::EdgeMap<double> > mwpm(graph,weight);
    mwpm.run();
    
    for (lemon::ListGraph::NodeIt it(graph); it != lemon::INVALID; ++it) {
        int i = graph.id(it);
        matches[i] = graph.id(mwpm.mate(it));
    }
}

void GTC_code::get_correction() {
    int ct = 0;
    for (int r = 0; r < params.n_rounds; r++) {
        for (int i = 0; i < n; i++) {
            if (r == 0) {
                if (syndrome[r](i)) {
                    verts_to_match[ct] = {i,r};
                    ct++;
                }
            } else {
                if (syndrome[r](i) != syndrome[r-1](i)) {
                    verts_to_match[ct] = {i,r};
                    ct++;
                }
            }
        }
    }
    
    int n_edges = (n_syndromes * (n_syndromes - 1)) / 2;
    PerfectMatching matching(n_syndromes, n_edges);

    for (int i = 0; i < n_syndromes; i++) {
        for (int j = i+1; j < n_syndromes; j++) {
            double dist = stab_dist(verts_to_match[i], verts_to_match[j]);
            if (dist != -1) {
                matching.AddEdge(i, j, dist);
            }
        }
    }

    matching.options.verbose = false;
    matching.Solve();

    for (int i = 0; i < n_syndromes; i++) {
        matches[i] = matching.GetMatch(i);
    }
    
//    lemon_match();
}

bool GTC_code::check_correction() {
    for (int i = 0; i < n_syndromes; i++) {
        if (verts_to_match[i].x() == -1) {
            continue;
        }
        
        int src = verts_to_match[i].x();
        int dest = verts_to_match[matches[i]].x();
        
        for (int i = 0; i < logicals.size(); i++) {
            logical_parities(i) ^= (logical_crossings[i](src,dest));
        }
        
        verts_to_match[matches[i]].x() = -1;
    }
    
    for (int i = 0; i < 4; i++) {
        if (logical_parities(i)) {
            return false;
        }
    }
    return true;
}

bool GTC_code::run_decode() {
    reset();
    
    for (int i = 0; i < params.n_rounds; i++) {
        single_round();
    }

    get_correction();
    return check_correction();
}

void GTC_code::print_error() {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            if (qubits_pos(x,y) != -1) {
                int i = qubits_pos(x,y);
                char err = 'I';
                if (errs_X(i) && errs_Z(i)) {
                    err = 'Y';
                } else if (errs_X(i)) {
                    err = 'X';
                } else if (errs_Z(i)) {
                    err = 'Z';
                }
                cout << err << "-";
            } else {
                cout << "+-";
            }
        }
        cout << endl;
        
        for (int x = 0; x < w; x++) {
            cout << "| ";
        }
        cout << endl;
    }
}

void GTC_code::print_syndrome(int round) {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            cout << "+-";
        }
        cout << endl;
        
        for (int x = 0; x < w; x++) {
            int syndrome_ind = 0;
            int id = to_rep(x, y);
            while(stabs(syndrome_ind,0) != id) {
                syndrome_ind++;
            }
            
            char syndrome_char = ' ';
            if (syndrome[round](syndrome_ind)) {
                syndrome_char = 'x';
            }
            cout << "|" << syndrome_char;
        }
        cout << endl;
    }
}

void GTC_code::print_identification() {
    int digits = int(log10(n));
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int identification = to_rep(x, y);
            int rep_digits = int(log10(identification));
            if (identification == 0) {
                rep_digits = 0;
            }
            
            cout << identification;
            for (int i = 0; i < digits - rep_digits; i++) {
                cout << "-";
            }
            cout << "-";
        }
        cout << endl;
        
        for (int x = 0; x < w; x++) {
            cout << "|";
            for (int i = 0; i < digits + 1; i++) {
                cout << " ";
            }
        }
        cout << endl;
    }
}

void GTC_code::print_stabs() {
    int digits = int(log10(n));
    
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            cout << "+";
            for (int i = 0; i < digits + 1; i++) {
                cout << "-";
            }
        }
        cout << endl;
        
        for (int x = 0; x < w; x++) {
            cout << "|";
            int syndrome_ind = 0;
            int id = to_rep(x, y);
            while(stabs(syndrome_ind,0) != id) {
                syndrome_ind++;
            }
            
            int syndrome_digits = int(log10(syndrome_ind));
            if (syndrome_ind == 0) {
                syndrome_digits = 0;
            }
            for (int i = 0; i < digits - syndrome_digits; i++) {
                cout << " ";
            }
            
            cout << syndrome_ind;
        }
        cout << endl;
    }
}
