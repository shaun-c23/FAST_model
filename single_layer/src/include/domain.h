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
#include <cmath>
#include <array>
using namespace std;

void linspace(double a, double b, int num, vector<double> &arr)
{
    num-=1;
    for (int i = 0; i <= num; i++)
    {
        arr.push_back(a + i * ( (b - a) / num ));
    }
    if (arr.back() != b){
        arr.push_back(b);
    }
}

vector< vector<double> > getdomain(int nx, int ny, int nz, double x1, double x2, double y1, double y2, double z1, double z2){
    vector<double> xi{};
    vector<double> yi{};
    vector<double> zi{};
    if (x1 != x2){
        linspace(x1, x2, nx, xi);
    }
    else{
        xi.push_back(x1);
    }
    if (y1 != y2){
        linspace(y1, y2, ny, yi);
    }
    else{
        yi.push_back(y1);
    }
    if (z1 != z2){
        linspace(z1, z2, nz, zi);
    }
    else{
        zi.push_back(z1);
    }

    vector < vector<double> > domain{};
    domain.push_back(xi);
    domain.push_back(yi);
    domain.push_back(zi);

    return domain;
}
