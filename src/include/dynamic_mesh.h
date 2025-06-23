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
#include <filesystem>

using namespace std;


void plotPointsSquareZ(double center_x, double center_y, double Ztop, double depth, double dr, double radius,
                      vector<double>& x_coords, vector<double>& y_coords, vector<double>& z_coords, bool modelMPdepth, double depthpoints) {
    double half_side_length = radius;
    double x_start = center_x - half_side_length;
    double y_start = center_y - half_side_length;
    depth += 0.1;
    double z = Ztop;
    if (modelMPdepth){
        for (double z = Ztop; z > -depth; z-= depthpoints){
            for (double x = x_start; x <= center_x + half_side_length; x += dr) {
                for (double y = y_start; y <= center_y + half_side_length; y += dr) {
                    x_coords.push_back(x);
                    y_coords.push_back(y);
                    z_coords.push_back(z);
                }
            }
        }
    }
    else{
        for (double x = x_start; x <= center_x + half_side_length; x += dr) {
            for (double y = y_start; y <= center_y + half_side_length; y += dr) {
                x_coords.push_back(x);
                y_coords.push_back(y);
                z_coords.push_back(z);
            }
        }
    }
}

void dyn_mesh(vector<double> &xdynam,vector<double> &ydynam,vector<double> &zdynam, double Xpos, double Ypos, double Zpos, double send, vector <vector<double>>  meltednodes, vector <double> meltedtime, double nodelengths,vector<double> &rads,vector <double> MPparams, bool modelMPdepth)
{

    double MPmult = MPparams[0];
    double MPspacing = MPparams[1];
    double depthspacing = MPparams[2];
    static double pi = 3.14159265358979323846;

    double dr = 1;      // Increment of radius
    double radius;
    if (nodelengths == 0){
        radius = 0;
    }
    else{
        if (meltednodes[meltednodes.size()-1][0] > 0.25){
            dr = MPspacing; 
            radius = (2.5*MPmult*sqrt(meltednodes[meltednodes.size()-1][0]/pi));
            int n = ceil(radius / dr);
            dr = radius / n ;
        }
        else{
            radius = 0.75;
            dr = MPspacing;
            int n = ceil(radius / dr);
            dr = radius / n ; 
            }
    }

    rads.push_back(radius);
    rads.push_back(dr);
    rads.push_back(dr);

    plotPointsSquareZ(Xpos, Ypos, Zpos, meltednodes[meltednodes.size()-1][1], dr, radius, xdynam, ydynam, zdynam, modelMPdepth,depthspacing);

}
