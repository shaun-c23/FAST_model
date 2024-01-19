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

    int indMP = 0;
    if (send > meltedtime[meltedtime.size()-1]){
        indMP = meltedtime.size()-1;
    }
	else if (send > meltedtime[indMP+1]){
		while (send> meltedtime[indMP]){
            indMP+=1;
			}
	}
		    else if (send> meltedtime[indMP]){
		        indMP +=1;
		    }



    double dr = 1;      // Increment of radius
    double radius;
    if (nodelengths == 0){
        radius = 0;
    }
    else{
        if (meltednodes[indMP][0] > 0.25){
            dr = MPspacing; 
            radius = (2.5*MPmult*sqrt(meltednodes[indMP][0]/M_PI));
            int n = ceil(radius / dr);
            dr = radius / n ;
        }
        else{
            radius = 0.25;
            dr = 0.05;
            int n = ceil(radius / dr);
            dr = radius / n ; 
            }
    }

    rads.push_back(radius);
    rads.push_back(dr);
    rads.push_back(dr);

    plotPointsSquareZ(Xpos, Ypos, Zpos, meltednodes[meltednodes.size()-1][1], dr, radius, xdynam, ydynam, zdynam, modelMPdepth,depthspacing); //meltednodes[indMP][1]

}
