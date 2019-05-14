%module GeantNPData
%{
#include "GeantNPData.hh"
%}

%include "std_vector.i"

namespace std {
    %template(vectord) vector<double>;
    %template(vectordd) vector<vector<double>>;
}

%include "GeantNPData.hh"
