%module GeantPPData
%{
#include "GeantPPData.hh"
%}

%include "std_vector.i"

namespace std {
    %template(vectord) vector<double>;
    %template(vectordd) vector<vector<double>>;
}

%include "GeantPPData.hh"
