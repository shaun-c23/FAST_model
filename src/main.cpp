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
#include <map>
#include <any>
#include <chrono>
#include <regex>

//include libraries
#include "include/read_laserpos.h" 
#include "include/gauss_quad.h" 
#include "include/mesh.h" 
#include "include/domain.h"
#include "include/runcode.h"
#include "include/read_input.h"

namespace fs = std::filesystem;
using namespace std;


int main(int argc, char* argv[]){

    auto start_time = std::chrono::high_resolution_clock::now();

    /*
    MESH DETAILS
    1. "static" - static grid of points (cumputationally unefficient and requires fine mesh in melt pool for accruate radiaiton)
    2. "dynamic" - static set of points dynmaically refined around melt pool
    3. "meltpool" - dynamic set of point just in and around the melt pool area
    4. "points" - select points but cannot iclude radiation since the meltpool is not modelled (can combine points with melt pool to include radiaiton effects)
    */

    fs::path currentPath = fs::current_path();
    fs::path ParentPath = currentPath.parent_path();
    string FilePath = ParentPath.string() + "/";

    ConfigPath configPath;
    ConfigVars setupVars;
    ConfigMPpoints setupMPpoints;
    ConfigDomain setupDomain;
    ConfigBC setupBC;
    vector<double> BCvals;

    filenames(FilePath, "inputs/setup.txt", configPath);
    vars(FilePath,configPath.Config, setupVars);
    MP_vars(FilePath,configPath.Meltpool, setupMPpoints);
    get_domain(FilePath,configPath.Domain, setupDomain);
    BCs(FilePath,configPath.BCPath, setupBC, BCvals);


    vector < vector<double> > domain{};
    domain = getdomain(setupDomain.nx,setupDomain.ny,setupDomain.nz,setupDomain.x1,setupDomain.x2,setupDomain.y1,setupDomain.y2,setupDomain.z1,setupDomain.z2);


    string type = configPath.Meshtype;
    int no_cores = setupVars.MaxThreads; //number of cores to run on. Generally, more cores will run quicker.

    double baseT = setupVars.InitialTemp; //initial temp
    double eta = setupVars.Efficiency; //laser absorption efficiency
    double sig = setupVars.Sigma; //Gaussian beam spot size in mm
    double delta = setupVars.Delta;
    double ambient_temp = setupVars.AmbientTemp; //ambient temperature for radiation
    int everyNframes = setupVars.everyNframes;


    //melt pool calculation points parameters
    double size_mult = setupMPpoints.BoxSizeMult; //multiplier to change size of dynamic meltpool calculation points - default is 1
    double point_spacing = setupMPpoints.SurfSpacing; //spacing of points along x-y within meltpool mesh in mm
    
    bool modelMPdepth = configPath.meltpooldepth; //set to true if you want to model the depth of the melt pool. Autmatically switches to to true in dynamic mesh if modelling depth (e.g. z2 != z1)
    double depthspacing = setupMPpoints.DepthSpacing;//spacing of points in depth of meltpool in mm - smaller spacing with be more accurate at cost of CPU speed
    vector <double> MPparams;
    MPparams.push_back(size_mult), MPparams.push_back(point_spacing), MPparams.push_back(depthspacing);

    double time_freq = setupVars.TimeStep; //defines how often (in seconds) you would like to write the the output
    double ex_time = setupVars.CoolingTime; //defines extra time (in seconds) you want to add to analysis (in addition to the amount of time in gcode)
    
    bool modelBC = configPath.modelBC;

    //ensure gcode file has be pre-processed with "gcode_interface.py" before proceeding
    string gcodefile = FilePath + configPath.BeamPath;//file name and location of processed gcode
    string datapath = "data/"; //location of data files to save for output
    fs::path data_dir = datapath;
    fs::path data_loc = FilePath / data_dir;

    //executes code
    runcode(baseT, ambient_temp, sig, eta, no_cores,domain,gcodefile,data_loc.string(),time_freq,ex_time,type,FilePath + configPath.Material,MPparams,modelMPdepth,modelBC,BCvals,delta,everyNframes);
    
    auto stop_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time);

    std::cout << "Execution time (s): " << duration.count() << "\n";

    return 0;
}
