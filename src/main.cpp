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
#include <cstdlib>

//include libraries
#include "include/read_laserpos.h" 
#include "include/gauss_quad.h" 
#include "include/mesh.h" 
#include "include/domain.h"
#include "include/runcode.h"

namespace fs = std::filesystem;
using namespace std;

int main(int argc, char* argv[]){

    /*
    MESH DETAILS
    1. "static" - static grid of points (cumputationally unefficient and requires fine mesh in melt pool for accruate radiaiton)
    2. "dynamic" - static set of points dynmaically refined around melt pool
    3. "meltpool" - dynamic set of point just in and around the melt pool area
    4. "points" - select points but cannot iclude radiation since the meltpool is not modelled (can combine points with melt pool to include radiaiton effects)
    */

    string type = "dynamic";

    vector < vector<double> > domain{};
    vector <double> MPparams;

    //meshed domain x-y is surface, z is depth
    if (type == "dynamic" or "static"){
        int nx = 75; //number of nodes along x
        int ny = 75; //number of nodes along y
        int nz =  1; //number of nodes along z
        double x1 = -2.5;//.5; //specify x lower bound
        double x2 = 12.5;//12.5; //specify x upper bound
        double y1 = -2.5;//-2.5; //specify y lower bound
        double y2 = 12.5;//12.5; //specify y upper bound
        double z1 = 0; //specify z lower bound
        double z2 = 0; //specify z upper bound
        domain = getdomain(nx,ny,nz,x1,x2,y1,y2,z1,z2);
    }

    //single points 
    else if (type == "points"){ 
        vector<double> xi{};
        vector<double> yi{};
        vector<double> zi{};

        xi = {5,-2.5,-5}; //specify x points
        yi = {0,5,10}; //specify y points
        zi = {-1.5875}; //specify z points
        domain.push_back(xi);
        domain.push_back(yi);
        domain.push_back(zi);
    }

    int no_cores = 6; //number of cores to run on. Generally, more cores will run quicker.

    double baseT = 250; //initial temp
    double eta = 0.45; //laser absorption efficiency
    double sig = 0.145; //Gaussian beam spot size in mm
    double ambient_temp = 30.0; //ambient temperature for radiation

    /*
    Material properties can be found in mat_props.h
    "Ti" - Ti-6Al-4V
    "Fe" - 4140 Steel
    "IN" - IN718
    */

    string mat = "Ti"; //define material based on code above or add data to mat_props.h with a new or existing identifier

    //melt pool calculation points parameters
    double size_mult = 2.5; //multiplier to change size of dynamic meltpool calculation points - default is 2.5
    double point_spacing = 0.075; //spacing of points along x-y within meltpool mesh in mm
    
    bool modelMPdepth = false; //set to true if you want to model the depth of the melt pool. Autmatically switches to to true in dynamic mesh if modelling depth (e.g. z2 != z1)
    double depthspacing = 0.025;//spacing of points in depth of meltpool in mm - smaller spacing with be more accurate at cost of CPU speed
    
    MPparams.push_back(size_mult), MPparams.push_back(point_spacing), MPparams.push_back(depthspacing);

    double time_freq = 0.01; //defines how often (in seconds) you would like to write the the output
    double ex_time = 0; //defines extra time (in seconds) you want to add to analysis (in addition to the amount of time in gcode)
    
    //ensure gcode file has be pre-processed with "gcode_interface.py" before proceeding
    string gcodefile = "gcode/gcodefiles/triangle_laser_pos.txt";//file name and location of processed gcode
    string datapath = "data/"; //location of data files to save for output
    
    fs::path currentPath = fs::current_path();
    fs::path parentPath = currentPath.parent_path();
    fs::path gcode_dir = gcodefile;
    fs::path data_dir = datapath;

    // Construct the path to the target directory
    fs::path gcode_loc = parentPath / gcode_dir;
    fs::path data_loc = parentPath / data_dir;

    //executes code
    runcode(baseT, ambient_temp, sig, eta, no_cores,domain,gcode_loc,data_loc,time_freq,ex_time,type,mat,MPparams,modelMPdepth); 
    
    return 0;
}