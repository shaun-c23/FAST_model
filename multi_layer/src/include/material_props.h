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

void readMatProps(const std::string& filename,
                      std::vector<double>& ktemps, std::vector<double>& k,
                      std::vector<double>& ctemps, std::vector<double>& c,
                      std::vector<double>& rho,
                      std::vector<double>& emissivity, std::vector<double>& meltT) {
    std::ifstream inputFile(filename);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file." << std::endl;
        return;
    }

    auto readVector = [&inputFile](std::vector<double>& vec) {
        inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '{');
        while (inputFile >> std::skipws >> std::ws && inputFile.peek() != '}') {
            double value;
            inputFile >> value;
            vec.push_back(value);
            inputFile >> std::ws;
            if (inputFile.peek() == ',')
                inputFile.ignore();
        }
        inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '}');
    };

    readVector(ktemps);
    readVector(k);
    readVector(ctemps);
    readVector(c);
    readVector(rho);
    readVector(emissivity);
    readVector(meltT);

    inputFile.close();
}

struct MaterialProperties {
    std::vector<double> ktemps, k, ctemps, c,rho, emissivity, meltT;
};

MaterialProperties getprops(const std::string mat) {
    MaterialProperties props;

    readMatProps(mat, props.ktemps, props.k, props.ctemps, props.c, props.rho, props.emissivity, props.meltT);
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