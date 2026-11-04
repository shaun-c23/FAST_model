// This software was authored by Shaun Cooke at The University of British Columbia – Vancouver.
// The research was supervised by Dr. Chad W. Sinclair and Dr. Daan M. Maijer,
// with funding support from the Natural Sciences and Engineering Research Council of Canada (NSERC).

/*
 * Copyright (c) 2024 Shaun Cooke. All rights reserved.
 * 
 * Author: Shaun Cooke <src@student.ubc.ca>
 * 
 * This software is provided for academic and non-commercial research purposes only.
 * Commercial use is strictly prohibited without prior written permission from the author.
 * Redistribution and modification are permitted for non-commercial use,
 * provided that proper credit is given and the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions, and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions, and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of the project nor the names of its contributors may be used
 *    to endorse or promote products derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY CLAIM,
 * DAMAGES, OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT, OR OTHERWISE,
 * ARISING FROM, OUT OF, OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
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

    return Aelem*boltz* sum / pow(1000.0, 2);
}

double iter_radiation(vector<double> tempT, double xnodes, double ynodes,double emis, double ambient_temp,double baseT, double eta, double power,double radloss){
    double oldflux;
    double newflux;
    double flux, pmult;
    double rad = 0;

    int max_iter = 100; //max number of iterations to convergence
    double converg_thresh = 0.5; //convergence threshold
                
    vector<double> newT(tempT.size(), 0);
    vector<double> newT2(tempT.size(), 0);

    double scaler = 0.95;
    flux = calculate_rad(tempT, xnodes, ynodes,emis,ambient_temp);
    oldflux = flux;
    newflux = flux;

    double startflux = flux;

    

    if (newflux >= (power/1000)*eta*(0.67) && power/1000 != 0){
        while (newflux >= (power/1000)*eta*(0.67)){
            pmult = scaler;//scales temperature field down
            getnewT(pmult, baseT, tempT, newT);
            newflux = calculate_rad(newT, xnodes, ynodes,emis,ambient_temp);
            scaler -=0.01;
        }
    }

    
    else{
        pmult = getpmult(radloss, power/1000 * eta);
        getnewT(pmult, baseT, tempT, newT);
        newflux = calculate_rad(newT, xnodes, ynodes,emis,ambient_temp);
    }
    

    bool track = true;
    while (track){
        for (int iter = 0; iter < max_iter; ++iter){
            double pmultnew = getpmult(newflux, power/1000 * eta);
            getnewT(pmultnew, baseT, newT, newT2);
            newflux = calculate_rad(newT2, xnodes, ynodes,emis,ambient_temp);
                            
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
            newflux = calculate_rad(newT, xnodes, ynodes,emis,ambient_temp);
            scaler -=0.01;
        }
    }
    return rad*emis;
}
