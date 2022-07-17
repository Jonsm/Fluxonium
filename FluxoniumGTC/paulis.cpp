//
//  paulis.cpp
//  FluxoniumGTC
//
//  Created by Jon on 5/6/22.
//

#include <stdio.h>
#include "paulis.hpp"

using namespace std;

string pauli_to_s(Pauli p) {
   switch(p) {
       case Pauli::X:
           return "X";
       case Pauli::Y:
           return "Y";
       case Pauli::Z:
           return "Z";
       case Pauli::I:
           return "I";
       default:
           return "NA";
   }
}
