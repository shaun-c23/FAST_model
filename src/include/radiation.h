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

double calculate_rad(vector <double> &T, double x_node_spacing, double y_node_spacing, double emissivity, double ambient_temp) {

    static double boltz = 5.67e-8; //Boltzmann constant
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

double iter_radiation(vector<double> tempT, double xnodes, double ynodes, double emis, double ambient_temp, double baseT, double eta, double power, double radloss) {
    double oldflux, newflux, flux, pmult;
    double rad = 0;

    int max_iter = 100; // max number of iterations to convergence
    double converg_thresh = 0.5; // convergence threshold
                
    vector<double> newT(tempT.size(), 0);
    vector<double> newT2(tempT.size(), 0);

    double scaler = 0.95;
    flux = calculate_rad(tempT, xnodes, ynodes, emis,ambient_temp);
    oldflux = flux;
    newflux = flux;

    double target_flux = (power / 1000) * eta * 0.67;

    // Adjust based on previous radiation loss
    if (radloss > 0) {
        newflux = (flux + radloss) / 2.0; // Smoothing between flux and previous radloss
    }

    // Scaling down if flux exceeds the target value
    if (newflux >= target_flux && power / 1000 != 0) {
        while (newflux >= target_flux) {
            pmult = scaler;  // scales temperature field down
            getnewT(pmult, baseT, tempT, newT);
            newflux = calculate_rad(newT, xnodes, ynodes, emis,ambient_temp);
            scaler -= 0.05;
            if (scaler < 0.1) { // Prevents over-scaling
                scaler = 0.1;
                break;
            }
        }

    } else {
        // Scaling using flux and previous radloss
        pmult = getpmult((flux * 0.67 + radloss) / 2.0, target_flux);
        getnewT(pmult, baseT, tempT, newT);
        newflux = calculate_rad(newT, xnodes, ynodes, emis,ambient_temp);
    }



    bool track = true;
    while (track) {
        for (int iter = 0; iter < max_iter; ++iter) {
            double pmultnew = getpmult(newflux, target_flux);
            getnewT(pmultnew, baseT, newT, newT2);
            newflux = calculate_rad(newT2, xnodes, ynodes, emis,ambient_temp);
                            
            if (abs(oldflux - newflux) < converg_thresh) {
                // Convergence reached
                rad = newflux;
                track = false;
                break;
            }
            oldflux = newflux;
        }
        
        // If not converged within max iterations, retry with smaller starting field
        if (track) {
            pmult = scaler;
            getnewT(pmult, baseT, tempT, newT);
            newflux = calculate_rad(newT, xnodes, ynodes, emis,ambient_temp);
            scaler -= 0.05;
            if (scaler < 0.1) {
                scaler = 0.1;
                break;
            }
        }
    }
    
    return rad;
}
