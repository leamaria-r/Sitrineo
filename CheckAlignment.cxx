// Loop to find the best alignment for the planes based on the same procedure that is written in MRaw.cxx
// all distances are computed wrt plane 1

#include <iostream>
#include<cstdio>
using namespace std;

void CheckAlignment()
{float plane2x = 1000; // initialization with value higher than threshold for x-axis
float plane3x = 1000;
float plane4x = 1000;

float plane2y = 1000; // initialization with value higher than threshold for y-axis
float plane3y = 1000;
float plane4y = 1000;

float thresholdP2 = 100; // when value lower than threshold, break the while loop and return the values
float thresholdP3 = 500;
float thresholdP4 = 500;

int i; // number of iterations, in order to avoid infinite loop 

for(i=0;i<30;i++){
    while(plane2x>thresholdP2 || plane3x>thresholdP3 || plane4x>thresholdP4)
    {plane2=meandiffx[0];
    plane3=meandiffx[1];
    plane4=meandiffx[2];

    MimosaAnalysis *gTAF = new MimosaAnalysis();
    gTAF->InitSession(runNumber,planeNumber,-1,configFile.Data());

    }
}

}
