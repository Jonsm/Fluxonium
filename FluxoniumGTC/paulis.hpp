//
//  paulis.hpp
//  FluxoniumGTC
//
//  Created by Jon on 4/23/22.
//

#ifndef paulis_hpp
#define paulis_hpp
#include <string>

enum class Pauli {X, Y, Z, I, NA};

std::string pauli_to_s(Pauli p);

#endif /* Paulis_h */
