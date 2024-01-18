//This software has been authored by Shaun Cooke with The University of British Columbia - Vancouver. 
//Research was supervised by Dr. Chad W. Sinclair and Dr. Daan M. Maijer with support by the National Sciences and Engineering Research Council of Canada (NSERC)
/*Copyright 2024, Shaun Cooke, All rights reserved.
*
* All Rights Reserved
*
* Authors: Shaun Cooke <src@student.ubc.ca>
*
* Redistribtion and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Redistributions of source code must retain the above copyright notice,
*	 this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. Neither the name of the project/files nor the names of its
*    contributors may be used to endorse or promote products derived from
*    this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*/

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <cmath>

using namespace std;

double getpmult(double flux, double powy){

    double pmult = ((powy - flux) /powy);
    if (pmult < 0){
        pmult = 0;
    }
    else if (pmult > 1){
        pmult = 1;
    }

    return pmult;
}

void getnewT(double pmult, double baseT, vector <double> &ogT, vector <double> &T){
    for (int i = 0; i < ogT.size(); i++) {
        T[i] = (ogT[i] - baseT)*pmult + baseT;
    }
}

double calculate_rad(vector <double> &T, double x_node_spacing, double y_node_spacing, double emissivity) {

    static double boltz = 5.67e-8; //Boltzmann constant
    double ambient_temp = 30.0;
    double Aelem = x_node_spacing*y_node_spacing; 

    double sum = 0.0;
    int size = T.size();

    for (int i = 0; i < size; i++) {
        double x = T[i];
        sum += ((pow(x + 273.15, 4) - pow(ambient_temp + 273.15, 4)));
    }

    double rad_pow =  Aelem * emissivity * boltz * sum / pow(1000.0, 2);
    return rad_pow;
}

double iter_radiation(vector<double> tempT, double xnodes, double ynodes,double emis, double baseT, double eta, double power ){
    double oldflux;
    double newflux;
    double flux, pmult;
    double rad = 0;

    int max_iter = 100; //max number of iterations to convergence
    double converg_thresh = 0.5; //convergence threshold
                
    vector<double> newT(tempT.size(), 0);
    vector<double> newT2(tempT.size(), 0);

    double scaler = 0.95;
    flux = calculate_rad(tempT, xnodes, ynodes,emis);
    oldflux = flux;
    newflux = flux;

    double startflux = flux;

    if (newflux >= (power/1000)*eta*(0.67) && power/1000 != 0){
        while (newflux >= (power/1000)*eta*(0.67)){
            pmult = scaler;//scales temperature field down
            getnewT(pmult, baseT, tempT, newT);
            newflux = calculate_rad(newT, xnodes, ynodes,emis);
            scaler -=0.05;
        }
    }

    else{
        pmult = getpmult(flux*0.67, power/1000 * eta);
        getnewT(pmult, baseT, tempT, newT);
        newflux = calculate_rad(newT, xnodes, ynodes,emis);
    }

    bool track = true;
    while (track){
        for (int iter = 0; iter < max_iter; ++iter){
            double pmultnew = getpmult(newflux, power/1000 * eta);
            getnewT(pmultnew, baseT, newT, newT2);
            newflux = calculate_rad(newT2, xnodes, ynodes,emis);
                            
            if (abs(oldflux - newflux) < converg_thresh) {
                //cout << "Converged after " << iter + 1 << " iterations - " << "flux = " << newflux << endl;
                rad = newflux;
                track = false;
                break;
            }
            oldflux = newflux;
        }
        // if not converged within max iterations, try again with smaller starting field
        if (track == true){
            pmult = scaler;
            getnewT(pmult, baseT, tempT, newT);
            newflux = calculate_rad(newT, xnodes, ynodes,emis);
            scaler -=0.05;
        }
    }
    return rad;
}