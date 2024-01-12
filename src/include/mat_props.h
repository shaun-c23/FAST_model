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
#include <cmath>

double interpolate(double x, const std::vector<double>& xData, const std::vector<double>& yData) {
    int size = xData.size();
    if (x <= xData[0]) {
        return yData[0];
    } else if (x >= xData[size - 1]) {
        return yData[size - 1];
    } else {
        int index = 0;
        while (xData[index] < x) {
            index++;
        }
        double t = (x - xData[index - 1]) / (xData[index] - xData[index - 1]);
        return yData[index - 1] + t * (yData[index] - yData[index - 1]);
    }
}

struct MaterialProperties {
    std::vector<double> ktemps, k, ctemps, c;
    double rho, emissivity, meltT;
};

MaterialProperties getprops(const std::string& mat) {
    MaterialProperties props;
    //Ti-6Al-4V
    if (mat == "Ti"){
        props.ktemps = {35.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 995.0, 1005.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 1600,1610,1620,1625};
        props.k = {7.0, 7.45, 8.75, 10.15, 11.35, 12.6, 14.18, 15.47, 17.55, 19.06, 19.41, 19.46, 21, 22.9, 23.7, 24.6, 25.8, 27.06,27.6,30.81,32.8};
        props.ctemps = {35, 400, 500, 700, 800, 900, 995, 1005, 1100, 1200, 1300, 1400, 1500, 1600,1610,1620,1625};
        props.c = {546*1e6, 629*1e6, 651*1e6, 684*1e6, 691.3*1e6, 677.4*1e6, 644.5*1e6, 642.7*1e6, 660*1e6, 678*1e6, 696*1e6, 714*1e6, 732*1e6, 750.7*1e6,758.5*1e6,803.2*1e6,831*1e6};
        props.rho = 4.2e-09;
        props.emissivity = 0.7;
        props.meltT = 1625;
    }
    //4140 Steel
    else if (mat == "Fe"){
        props.ktemps = {25,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400};
        props.k = {46,44,41,39,37.5,35,33,32.5,30,29,28,27.5,30,30,30};
        props.ctemps = {0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400};
        props.c = {475*1e6,475*1e6,525*1e6,540*1e6,560*1e6,640*1e6,730*1e6,825*1e6,825*1e6,600*1e6,600*1e6,600*1e6,600*1e6,600*1e6,600*1e6};
        props.rho = 7.85e-09;
        props.emissivity = 0.55;
        props.meltT = 1460;
    }
    //IN718
    else if (mat == "IN"){
        props.ktemps = {25,100,300,500,700,1350};
        props.k = {11.4,12.5,14.0,15.5,21.5,31.3};
        props.ctemps = {25,100,300,500,700,1350};
        props.c = {427*1e6,441*1e6,481*1e6,521*1e6,601*1e6,691*1e6};
        props.rho = 8.146e-09;
        props.emissivity = 0.5;
        props.meltT = 1255;
    }

    return props;
}

double calculate_weighted_average(double T0, double T_range, vector<double> temps, vector<double> vals) {
    int num_intervals = 100;  // Number of intervals for integration

    // Generate temperature values with logarithmic spacing
    vector<double> T_values;
    for (int i = 0; i < num_intervals; i++) {
        double T = T0 * pow((T0 + T_range) / T0, i / (double)(num_intervals - 1));
        T_values.push_back(T);
    }
    vector<double> values;

    for (double T : T_values) {
            values.push_back(interpolate(T, temps, vals));
    }

    double weighted_avg = 0.0;
    double total_weight = 0.0;

    for (int i = 0; i < num_intervals; i++) {
        weighted_avg += values[i] * T_values[i];
        total_weight += T_values[i];
    }

    weighted_avg/= total_weight;

    return weighted_avg;
}