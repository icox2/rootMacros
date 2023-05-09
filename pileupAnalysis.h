#include "TGraph.h"
#include "TROOT.h"
#include "TF1.h"
#include "functions.hpp"
#include <vector>

int pileupAnalysis(vector<double> trace){
    vector<double> derivative, normalized;
    int numPile=0;
    if(trace.size()==0){
        return -1;
    }

    normalized = normalization(trace);
    derivative = der(normalized);
    for(int i=1;i<derivative.size()-1;i++){
        if(derivative.at(i)>0.1 && derivative.at(i)>derivative.at(i-1) && derivative.at(i)>derivative.at(i+1))
            numPile++;
            if(numPile==5) break;
    }
    return numPile;
}