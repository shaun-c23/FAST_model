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


vector<double> fitPolynomial(const vector<double>& x, const vector<double>& y, int degree) {
    int numSamples = x.size();

    // Create matrices A and B for the least squares method
    vector<vector<double>> A(degree + 1, vector<double>(degree + 1, 0.0));
    vector<double> B(degree + 1, 0.0);

    // Populate matrices A and B
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k < numSamples; ++k) {
                A[i][j] += pow(x[k], i + j);
            }
        }

        for (int k = 0; k < numSamples; ++k) {
            B[i] += y[k] * pow(x[k], i);
        }
    }

    // Solve the linear system to get coefficients
    vector<double> coefficients(degree + 1, 0.0);

    // Perform Gaussian elimination
    for (int i = 0; i <= degree; ++i) {
        for (int j = i + 1; j <= degree; ++j) {
            double factor = A[j][i] / A[i][i];
            for (int k = 0; k <= degree; ++k) {
                A[j][k] -= factor * A[i][k];
            }
            B[j] -= factor * B[i];
        }
    }

    // Back substitution
    for (int i = degree; i >= 0; --i) {
        coefficients[i] = B[i];
        for (int j = i + 1; j <= degree; ++j) {
            coefficients[i] -= A[i][j] * coefficients[j];
        }
        coefficients[i] /= A[i][i];
    }

    return coefficients;
}

double evaluatePolynomial(const vector<double>& coefficients, double x) {
    int degree = coefficients.size() - 1;

    double result = 0;
    for (int j = 0; j <= degree; ++j) {
        result += coefficients[j] * pow(x, j);
    }

    return result;
}