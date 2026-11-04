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
#include <iomanip>
#include <thread>
#include <boost/math/special_functions/gamma.hpp>

//include libraries
#include "fit_poly.h"
#include "material_props.h"
#include "radiation.h"
#include "dynamic_mesh.h"

namespace fs = std::filesystem;
using namespace std;

void showProgressBar1(double current, double total, int updateInterval) {
    static int lastUpdated = -1; // Track the last progress update
    int percentage = static_cast<int>(static_cast<double>(current) / total * 100);

    // Only update if the percentage has increased by the update interval
    if (percentage / updateInterval > lastUpdated) {
        lastUpdated = percentage / updateInterval;

        int barWidth = 50;
        double progress = static_cast<double>(current) / total;

        if (std::isnan(progress)) {
            progress = 0.0;
        }

        std::cout << "\rProgress: [";
        int pos = static_cast<int>(barWidth * progress);
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "█";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << std::fixed << std::setprecision(1) << progress * 100.0 << "%";
        std::cout.flush();
    }
}

double mirrorXValue(double x, double axis) {
    return axis - (x - axis);
}


double integrand_dim(double prefa, double pwr, double pmlt,double sp, double X,double Y, double Z,double send,double alpha, double sig, double Vx, double Vy, double posX, double posY, double posZ, double slim, bool BC, vector<double> BCvals, double delta)
{

    double phix,phiy,phiz;
    phix = 12*alpha*(send-sp)+pow(sig,2);
    phiy = phix;
    phiz = 12*alpha*(send-sp)+pow(delta,2);

    double spp,Xb,Yb,Zb;
    double Xb1, Xb2, Yb1, Yb2, Zbot;

    spp = send - sp;

    Xb = posX + (Vx)*(sp-slim);
    Yb = posY  + (Vy)*(sp-slim);
    Zb = posZ;

    if (!BC){
        return prefa*pwr*pmlt*(1/(sqrt(phix*phiy*phiz)))*(exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z-Zb),2)/phiz));
    }
    else{
       
        Xb1 = mirrorXValue(posX, BCvals[0]) + (-Vx)*(sp-slim);
        Xb2 = mirrorXValue(posX, BCvals[1]) + (-Vx)*(sp-slim);

        Yb1 = mirrorXValue(posY, BCvals[2]) + (-Vy)*(sp-slim);
        Yb2 = mirrorXValue(posY, BCvals[3]) + (-Vy)*(sp-slim);

        Zbot = BCvals[4];

        return prefa * pwr * pmlt * (1.0 / sqrt(phix * phiy * phiz)) * (
       exp(-3.0 * pow(X - Xb,  2) / phix
         - 3.0 * pow(Y - Yb,  2) / phiy
         - 3.0 * pow(Z - Zb,  2) / phiz)

     + exp(-3.0 * pow(X - Xb,  2) / phix
         - 3.0 * pow(Y - Yb,  2) / phiy
         - 3.0 * pow(Z - (2.0 * Zbot - Zb), 2) / phiz)

     + exp(-3.0 * pow(X - Xb,  2) / phix
         - 3.0 * pow(Y - Yb1, 2) / phiy
         - 3.0 * pow(Z - Zb,  2) / phiz)

     + exp(-3.0 * pow(X - Xb,  2) / phix
         - 3.0 * pow(Y - Yb2, 2) / phiy
         - 3.0 * pow(Z - Zb,  2) / phiz)

     + exp(-3.0 * pow(X - Xb1, 2) / phix
         - 3.0 * pow(Y - Yb,  2) / phiy
         - 3.0 * pow(Z - Zb,  2) / phiz)

     + exp(-3.0 * pow(X - Xb2, 2) / phix
         - 3.0 * pow(Y - Yb,  2) / phiy
         - 3.0 * pow(Z - Zb,  2) / phiz)

     + exp(-3.0 * pow(X - Xb,  2) / phix
         - 3.0 * pow(Y - Yb1, 2) / phiy
         - 3.0 * pow(Z - (2.0 * Zbot - Zb), 2) / phiz)

     + exp(-3.0 * pow(X - Xb,  2) / phix
         - 3.0 * pow(Y - Yb2, 2) / phiy
         - 3.0 * pow(Z - (2.0 * Zbot - Zb), 2) / phiz)

     + exp(-3.0 * pow(X - Xb1, 2) / phix
         - 3.0 * pow(Y - Yb,  2) / phiy
         - 3.0 * pow(Z - (2.0 * Zbot - Zb), 2) / phiz)

     + exp(-3.0 * pow(X - Xb2, 2) / phix
         - 3.0 * pow(Y - Yb,  2) / phiy
         - 3.0 * pow(Z - (2.0 * Zbot - Zb), 2) / phiz));
       
        /*
        return prefa*pwr*pmlt*(1/(sqrt(phix*phiy*phiz)))*(exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z-Zb),2)/phiz)
                                                           +exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z+Zb)-Zbot*2,2)/phiz)
                                                           +exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb1),2)/phiy-3*pow((Z-Zb),2)/phiz)
                                                           +exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb2),2)/phiy-3*pow((Z-Zb),2)/phiz)
                                                           +exp(-3*pow((X-Xb1),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z-Zb),2)/phiz)
                                                           +exp(-3*pow((X-Xb2),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z-Zb),2)/phiz));
        */                                               
    }

}

void linspace_s(double a, double b, int num, vector<double> &arr)
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

template <typename T>
void clearVectors(std::vector<std::vector<T>>& vec) {
    for (auto& innerVec : vec) {
        innerVec.clear();
    }
    vec.clear();
}

double d_first_term_taylor(double alpha, double z, double t, double t_start, double t_end){
    double prefactor = 1 / (pow((12 * alpha),(1.5)));
    
    //Incomplete Gamma function: Gamma(0.5, x)
    double x_start = (t - t_start);
    double x_end = (t - t_end);
    
    //Compute the integral result
    double integral_value = prefactor * (pow(x_start,-1.5) - pow(x_end,-1.5));
    
    return integral_value;
}

double compute_integral_with_incomplete_gamma(double alpha, double z, double t, double t_start, double t_end) {
    double b = std::pow(z, 2) / (4 * alpha);
    double prefactor = 1 / (std::pow(12 * alpha, 1.5) * std::sqrt(b));
    double x_start = b / (t - t_start);
    double x_end = b / (t - t_end);

    double gamma_start = boost::math::gamma_p(0.5, x_start) * boost::math::tgamma(0.5);
    double gamma_end = boost::math::gamma_p(0.5, x_end) * boost::math::tgamma(0.5);

    double integral_value = -1 * (gamma_start - gamma_end);
    return integral_value;
}

bool areAllZValuesSame(const std::vector<double>& zValues) {
    if (zValues.empty()) {
        // If the vector is empty, technically all values are the same (as there are none)
        return true;
    }

    int firstZValue = zValues[0]; // Get the first z value

    // Compare each z value with the first one
    for (size_t i = 1; i < zValues.size(); ++i) {
        if (zValues[i] != firstZValue) {
            return false; // If any z value is different, return false
        }
    }

    return true; // If all z values are the same, return true
}

double average_diff(const std::vector<double>& start_times) {
    if (start_times.size() < 2) {
        // If there are fewer than two elements, there are no differences to calculate.
        return 0.0;
    }

    std::vector<double> diffs;
    for (size_t i = 1; i < start_times.size(); ++i) {
        diffs.push_back(start_times[i] - start_times[i - 1]);
    }

    double sum_diffs = std::accumulate(diffs.begin(), diffs.end(), 0.0);
    return sum_diffs / diffs.size();
}

struct Layer {
    double start_time;
    double end_time;
    double height;
    int layer_index;
    double average_power;
};

void processLayers(const std::vector<double>& z, const std::vector<double>& p, const std::vector<double>& t, std::vector<Layer>& layers) {
    double prev_z = -std::numeric_limits<double>::infinity(); // Initialize to negative infinity
    double current_start = -1.0; // Start time of the current layer
    bool layer_in_progress = false; // Indicates if a valid layer segment is being processed
    int layer_start_index = -1; // To track the index of the current layer start

    double power_sum = 0.0; // Sum of power for the current layer
    int power_count = 0;    // Count of power values for averaging

    for (size_t index = 0; index < z.size(); ++index) {
        double z_ = z[index];
        double p_ = p[index];
        double t_ = t[index];

        // Check if a new layer starts
        if (z_ > prev_z) {
            if (layer_in_progress) {
                // Close the previous layer
                double current_end = t[index - 1];
                double average_power = power_count > 0 ? power_sum / power_count : 0.0;
                layers.push_back({current_start, current_end, prev_z, layer_start_index, average_power});
            }

            // Start a new layer
            prev_z = z_;
            current_start = -1.0; // Reset start time for the new layer
            layer_in_progress = false;
            layer_start_index = index; // Store the index of the new layer
            power_sum = 0.0; // Reset power sum for the new layer
            power_count = 0; // Reset power count for the new layer
        }

        // Update layer progress if power is non-zero
        if (p_ != 0.0) {
            if (current_start < 0) {
                current_start = t_; // Set the start time
            }
            power_sum += p_; // Accumulate power
            power_count++;    // Increment count
            layer_in_progress = true;
        }
    }

    // Handle the last layer if it is still in progress
    if (layer_in_progress) {
        double current_end = t.back(); // End time is the last time in the dataset
        double average_power = power_count > 0 ? power_sum / power_count : 0.0;
        layers.push_back({current_start, current_end, prev_z, layer_start_index, average_power});
    }
}

void runcode(double baseT, double ambient_temp, double sig, double eta, int nthreads, vector< vector<double> > domain,vector< vector<double> > SPs,vector<string> gcodes, string datapath, double time_freq, double ex_time, string type, string mat,vector<double> MPparams, bool modelMPdepth,bool ifSP,bool modelBC, bool secondbeamON, vector<double> BCvals, double delta, int everyNframes, bool LayRemove, bool corr, bool setlayers, int nlayers, double C){

    fs::create_directory(datapath);

    string filename = "MP1";
    string extension =  ".csv";
    string rtn = datapath +'/' + filename + extension;
    ofstream myfile(rtn);

    if (!LayRemove){
        corr = false;
    }

    static double pi = 3.14159265358979323846;

    vector<double> x,vx, vy, y, z, s, p;
    if (filesystem::exists(gcodes[0])) {
        laserpath(gcodes[0], s,x,vx,y,vy,z,p); //read gcode file
    } 
    else {
        cout << "Beam path 1 file not found... Program terminating." << endl;
        std::exit(EXIT_FAILURE);
    }





    // 1️⃣ Compute dz
    vector<double> dz(z.size() - 1);
    for (size_t i = 0; i < z.size() - 1; i++) {
        dz[i] = z[i + 1] - z[i];
    }

    // 2️⃣ Find jump indices (abs(dz) > 0.14)
    vector<size_t> jump_indices;
    for (size_t i = 0; i < dz.size(); i++) {
        if (abs(dz[i]) > 0.14) jump_indices.push_back(i);
    }

    // 3️⃣ Store jump times (time after the jump)
    vector<double> jump_times;
    for (auto idx : jump_indices) {
        jump_times.push_back(s[idx + 1]);
    }

    // 4️⃣ Calculate r^2 = x^2 + y^2
    vector<double> r2(x.size());
    for (size_t i = 0; i < x.size(); i++) {
        r2[i] = x[i]*x[i] + y[i]*y[i];
    }

    // 5️⃣ Find unique z values
    vector<double> unique_z = z;
    sort(unique_z.begin(), unique_z.end());
    auto last = unique(unique_z.begin(), unique_z.end());
    unique_z.erase(last, unique_z.end());

    // 6️⃣ Compute max r2 per z layer
    vector<double> max_r2_per_z;
    for (auto z_val : unique_z) {
        double max_r2 = 0;
        for (size_t i = 0; i < z.size(); i++) {
            if (z[i] == z_val) {
                if (r2[i] > max_r2) max_r2 = r2[i];
            }
        }
        max_r2_per_z.push_back(max_r2);
    }

    // 7️⃣ Main calculation
    vector<size_t> L;
    for (size_t i = 15; i < jump_times.size(); i++) L.push_back(i);

    double Threshold = 0.85;
    vector<size_t> N_;

    for (auto l : L) {
        vector<double> tmp;

        for (size_t i = 1; i < l; i++) {
            double val = exp(-max_r2_per_z[l - i] / (4.0 * 30.0 * (jump_times[l] - jump_times[l - i])));
            tmp.push_back(val);
        }

        // Find first index where tmp > Threshold
        auto it = find_if(tmp.begin(), tmp.end(), [&](double v){ return v > Threshold; });
        if (it != tmp.end()) {
            N_.push_back(distance(tmp.begin(), it));
        } else {
            N_.push_back(l);
        }
    }
    
    vector<int> gcodeinds;
    vector<vector<double>> xm,vxm, vym, ym, zm, sm, pm; 
    xm.resize(5);
    vxm.resize(5);
    vym.resize(5);
    ym.resize(5);
    zm.resize(5);
    sm.resize(5);
    pm.resize(5);
    if (secondbeamON){
        for (int kd = 0; kd < gcodes.size()-1; kd++){

            if (filesystem::exists(gcodes[kd+1])) {
                laserpath(gcodes[kd+1], sm[kd],xm[kd],vxm[kd],ym[kd],vym[kd],zm[kd],pm[kd]);
                gcodeinds.push_back(kd);
            } 
            else {
                cout << "Beam " << kd+2 << " file not found..." << endl;
            }
        }

        if (gcodeinds.size() == 0){
            cout << "Running only Beam 1" << endl;
            secondbeamON = false;
        }
        else{
            cout << "Running:" << endl;
            cout << "Beam 1" << endl;
            for (int uu = 0; uu < gcodeinds.size(); uu++){
                cout << "Beam " << gcodeinds[uu]+2 << endl;
            }
        }
    }
    else{
        cout << "Running Beam 1" << endl;
    }

    vector<double> layer_inds1,end_times1, start_times1, layer_heights1,layer_pows1,layer_nums1;

    std::vector<Layer> layers1; // To store layer information
    processLayers(z, p, s, layers1);

    for (size_t i = 0; i < layers1.size(); ++i) {
        layer_inds1.push_back(layers1[i].layer_index);
        start_times1.push_back(layers1[i].start_time);
        end_times1.push_back(layers1[i].end_time);
        layer_heights1.push_back(layers1[i].height);
        layer_pows1.push_back(layers1[i].average_power);
        layer_nums1.push_back(i+1);
        //cout << layers1[i].start_time << endl;
        //cout << layer_nums1[i] << endl;
    }


    vector<vector<Layer>> layers(5);
    std::vector<std::vector<double>> layer_inds(5);
    std::vector<std::vector<double>> end_times(5);
    std::vector<std::vector<double>> start_times(5);
    std::vector<std::vector<double>> layer_heights(5);
    std::vector<std::vector<double>> layer_pows(5);
    std::vector<std::vector<double>> layer_nums(5);
    if (secondbeamON){
        for (int kd = 0; kd < gcodeinds.size(); kd++){
            int indexval = gcodeinds[kd];

            processLayers(zm[indexval], pm[indexval], sm[indexval], layers[indexval]);
            for (size_t i = 0; i < layers[indexval].size(); ++i) {
                layer_inds[indexval].push_back(layers[indexval][i].layer_index);
                start_times[indexval].push_back(layers[indexval][i].start_time);
                end_times[indexval].push_back(layers[indexval][i].end_time);
                layer_heights[indexval].push_back(layers[indexval][i].height);
                layer_pows[indexval].push_back(layers[indexval][i].average_power);
                layer_nums[indexval].push_back(i+1);
                //cout << layers[indexval][i].start_time << endl;
            }
        }
    }

    bool layerON = true;
    if (areAllZValuesSame(layer_heights1)){
        layerON = false;
    }


    //node spaces
    vector<node> course_nodes;
    vector<node> MP_nodes; 
    vector<int> nodeinds;
    vector<int> nodeindsm;
    nodeindsm.resize(5);

    int nodenum = 1;
    if (type == "static" || type == "dynamic") {
        for (int i = 0; i < domain[0].size(); i++){
            for (int j = 0; j <domain[1].size(); j++){
                for (int k = 0; k < domain[2].size(); k++){
                    double x = domain[0][i];
                    double y = domain[1][j];
                    double z = domain[2][k];
                    course_nodes.push_back(make_node(nodenum,x,y,z,0,0));
                    nodenum+=1;
                }
            }
        }
    }

    if (ifSP || type == "points" || type == "dynamic-points"){
        for (int i = 0; i < SPs[0].size(); i++){
            double x = SPs[0][i];
            double y = SPs[1][i];
            double z = SPs[2][i];
            course_nodes.push_back(make_node(nodenum,x,y,z,0,0));
            nodenum+=1;
        }
    }

    //get material properties
    MaterialProperties props = getprops(mat);

    double rho = props.rho[0];
    double meltT = props.meltT[0];
    double emis = props.emissivity[0];

    double k = calculate_weighted_average(baseT, meltT-baseT, props.ktemps, props.k);
    double cp = calculate_weighted_average(baseT, meltT-baseT, props.ctemps, props.c);
    double a = k / (rho * cp);

    double k2 = calculate_weighted_average(baseT, baseT, props.ktemps, props.k);
    double cp2 = calculate_weighted_average(baseT, baseT, props.ctemps, props.c);
    double a2 = k2 / (rho * cp2);

    double k_T0 = interpolate(baseT, props.ktemps, props.k);
    double cp_T0 = interpolate(baseT, props.ctemps, props.c);

    vector<double> cavgs, kavgs;
    for (int Titer : props.ktemps){
        double Trange = Titer-baseT;
        double val = calculate_weighted_average(baseT, Trange, props.ktemps, props.k);
        kavgs.push_back(val);
    }
    for (int Titer : props.ctemps){
        double Trange = Titer-baseT;
        double val = calculate_weighted_average(baseT, Trange, props.ctemps, props.c);
        cavgs.push_back(val);
    }

    // Perform polynomial curve fit
    vector<double> c_coefficients = fitPolynomial(props.ctemps, cavgs, 4); //fit weighted avg cp values with 4th order polynomial
    vector<double> k_coefficients = fitPolynomial(props.ktemps, kavgs, 1); //fit weighted avg cp values with 1st order polynomial
    
    //prefactor 
    double prefactor = 2*eta/(rho*cp*pow((pi/3),1.5)); //non dimen

    //define output time frames
    vector<double> si;

    double s1 = s[0];
    double s2 = (s[s.size()-1] + ex_time);
    //s1 = 1026.98;//1225.08;//
    //s1 = 1046.73;

    int N = ceil((s2) / time_freq);
    linspace_s(s1, s2, N, si);

    /*make a vector of vecotrs with gauss quad points from n=2 to n=16*/
    vector< vector <double> > roots;
    vector< vector <double> > weights;
    double max_n = 16;
    for (int nval = 1; nval <max_n + 1; nval++){
        vector<double> temp_roots, temp_weights;
        temp_roots.clear();
        temp_weights.clear();
        findRoots(nval,temp_roots,temp_weights);
        roots.push_back(temp_roots);
        weights.push_back(temp_weights);
    }
    
    // intialize vectors
    vector<double> int_points, spoints, vxpoints, vypoints, Ppoints, posxpoints, posypoints, poszpoints, prefpoints, apoints;
    //vector<double> int_pointsm, spointsm, vxpointsm, vypointsm, Ppointsm, posxpointsm, posypointsm, poszpointsm, prefpointsm, p_pointsm, weights_totm;
    vector<double> weights_tot, p_points, melted;
    vector<double> currT, currx, curry, currz, rad_pow;
    vector< vector<double> > Tatdt_course, Tatdt_fine, Tatdt, MPatdt, xatdt, yatdt, zatdt;
    vector <double> probe1;

    vector< vector<double> > int_pointsm, spointsm, vxpointsm, vypointsm, Ppointsm, posxpointsm, posypointsm, poszpointsm, prefpointsm, p_pointsm, weights_totm;

    // Resize the outer vectors to have 5 elements (each element is an inner vector)
    int_pointsm.resize(5);
    spointsm.resize(5);
    vxpointsm.resize(5);
    vypointsm.resize(5);
    Ppointsm.resize(5);
    posxpointsm.resize(5);
    posypointsm.resize(5);
    poszpointsm.resize(5);
    prefpointsm.resize(5);
    p_pointsm.resize(5);
    weights_totm.resize(5);

    vector<double>  MPlen;
    MPlen.push_back(0);
    rad_pow.push_back(0);

    vector<vector<double>> rad_powm;
    rad_powm.resize(5);
    rad_powm[0].push_back(0);
    rad_powm[1].push_back(0);
    rad_powm[2].push_back(0);
    rad_powm[3].push_back(0);
    rad_powm[4].push_back(0);

    vector< vector<double> > meltednodes;
    vector< double > meltedtime;
    vector< double > meltedwidth;
    meltednodes.push_back({0,0});
    meltedtime.push_back(0);
    meltedwidth.push_back(0);

    double nodelengths = 0;
    int updateInterval = 2;

    vector< vector<vector<double>> > meltednodesm;
    vector<double> nodelengthsm(5, 0.0);
    meltednodesm.resize(5);
    meltednodesm[0].push_back({0,0});
    meltednodesm[1].push_back({0,0});
    meltednodesm[2].push_back({0,0});
    meltednodesm[3].push_back({0,0});
    meltednodesm[4].push_back({0,0});



    vector< vector<double> > MPcoords;
    vector< vector<vector<double> > > MPcoordsm;
    MPcoordsm.resize(5);
    MPcoords.push_back({0,0,0});
    MPcoordsm[0].push_back({0,0,0});
    MPcoordsm[1].push_back({0,0,0});
    MPcoordsm[2].push_back({0,0,0});
    MPcoordsm[3].push_back({0,0,0});
    MPcoordsm[4].push_back({0,0,0});

    // get G1 feed rate
    double check = 0;
    double V;
    double Vd;
    double P;
    while (check < p.size()){
        if (p[check]!=0){
            V = sqrt(pow(vx[check],2)+pow(vy[check],2))*(sig/a);
            Vd = round(sqrt(pow(vx[check],2)+pow(vy[check],2)));
            P = p[check]/1000;
            break;
        }
        check+=1;
    }

    int numFrames = s2/time_freq + 1;
    double startValue = 0.0;
    double endValue = s2;

    int startidd = 0;
    double newspp = s1;

    int layerNumber = 0;
    int startLayer = 0;
    int currentIndex = 1;
    bool hasIncremented = false; // Flag to track whether layerNumber has been incremented

    double zheight = layer_heights1[0];
    double averagerad = 0;

    vector <double> inds;
    vector <double> avgrads1;
    vector <vector<double>> avgrads;
    avgrads.resize(5);

    double prevspp = 0;
    double prevT = 0;
    double prevlay = 0;
    int idd = 0;

    double send = s1;
    int cnt = 0;

    for (int j = 0; j < numFrames; j++){
    
        int_points.clear();
        weights_tot.clear();
        spoints.clear();
        posxpoints.clear();
        posypoints.clear();
        poszpoints.clear();
        vxpoints.clear();
        vypoints.clear();
        Ppoints.clear();
        p_points.clear();

        // Clear contents of inner vectors but keep the outer size as 5
        for (int i = 0; i < 5; i++) {
            int_pointsm[i].clear();
            spointsm[i].clear();
            vxpointsm[i].clear();
            vypointsm[i].clear();
            Ppointsm[i].clear();
            posxpointsm[i].clear();
            posypointsm[i].clear();
            poszpointsm[i].clear();
            prefpointsm[i].clear();
            p_pointsm[i].clear();
            weights_totm[i].clear();
        }


        if (layerON){

            int currentIdx = layer_inds1[currentIndex];
	        int lastind = layer_inds1.size();

            // Compare the time value 'tn' with the corresponding value in 's' at the current index
            if (send < s[currentIdx]) {
                // Increment the layer number only once until the next 's' value
                if (!hasIncremented) {
                    layerNumber++;
                    hasIncremented = true;
                }
            } 

            else if (currentIndex >= lastind){
                if (!hasIncremented) {
                    layerNumber++;
                    hasIncremented = true;
                    //cout << layerNumber << endl;
                }
            }
            else {
                hasIncremented = false;
                currentIndex++;
                //cout << layerNumber << endl;
            }
            
            
            double tstart = start_times1[idd];
            double tend = end_times1[idd];
            zheight = layer_heights1[layerNumber-1];
            double pwrr = layer_pows1[idd];
            double preflayC = 2*pwrr*eta/(rho*cp2*pow((pi/3),1.5));

            if (LayRemove){
            if (!setlayers){
                //if (send > tend && corr && layerNumber > prevlay){
                if (send > tend && layerNumber > prevlay){
                    double slope = -1 * preflayC * d_first_term_taylor(a2, 0, send, tstart,tend);
                        if (slope <= C){
                            newspp = start_times1[idd+1];
                            cout << "layer : " << layerNumber <<" layer removed: " << idd+1 << endl;
                            idd++;
                        }
                }
                prevlay = layerNumber;

            }
            
            
           
            else{
                //if (send > tend && corr && layerNumber > prevlay){
                if (layerNumber > 16){

                    int nlayer = N_[idd];
                    
                    double startlayer = layerNumber - nlayer;
                    //cout << startlayer << endl;
                        if (startlayer > 0 && send > tend && layerNumber > prevlay){
                            newspp = start_times1[idd+1];
                            //cout << newspp << endl;
                            cout << "layer : " << layerNumber <<" layer removed: " << idd+1 << endl;
                            idd++;
                        }
                        prevlay = layerNumber;
                }

            }
            

            }
            
            averagerad = 0;
            if (!hasIncremented){
                double sum = std::accumulate(rad_pow.begin(), rad_pow.end(), 0.0,
                [](double a, double b) { return b != 0.0 ? a + b : a; });

                // Count non-zero elements
                size_t non_zero_count = std::count_if(rad_pow.begin(), rad_pow.end(),
                [](double value) { return value != 0.0; });

                // Calculate the average, ensuring we don't divide by zero
                averagerad = (non_zero_count > 0) ? sum / non_zero_count : 0.0;
                avgrads1.push_back(averagerad*1000);

                
                if (secondbeamON){
                    for (int kd = 0; kd < gcodeinds.size(); kd++){
                        int indexval = gcodeinds[kd];
                        double sum = std::accumulate(rad_powm[indexval].begin(), rad_powm[indexval].end(), 0.0,
                        [](double a, double b) { return b != 0.0 ? a + b : a; });

                        // Count non-zero elements
                        size_t non_zero_count = std::count_if(rad_powm[indexval].begin(), rad_powm[indexval].end(),
                        [](double value) { return value != 0.0; });

                        // Calculate the average, ensuring we don't divide by zero
                        averagerad = (non_zero_count > 0) ? sum / non_zero_count : 0.0;
                        avgrads[indexval].push_back(averagerad*1000);

                    }
                }
            }
        }

        double spp = s1; //0
        double dspp;
	    int ind = 1;

        if (layerON && idd > 0){
            spp = newspp;
        }

        int indMP = 0;
        double spp_;

        //print the model progress
        showProgressBar1(send, s2, updateInterval);

        double mult = a / pow(sig,2);

        while (spp < send){

            spp_ = send-spp;

            double ux;
            double lx;
            double uxm;
            double lxm;
            vector<int> indms;
            indms.resize(5);

            double ds_spot = pow(2,floor(sqrt(log2(1+spp_*mult))));
            int gauss_order = max(2, int(max_n / (ds_spot)));
            double pmult;
            vector<double> pmultm(5, 0.0);

            if (V < 1){
                V = 1;
            }

            dspp = (sqrt((12*spp_*mult+1)*log(sqrt(2)))/(5*V)); // integration segment modified from Stump et al. (Applied Mathematical Modelling, 2019)
            dspp = dspp / mult;
            //dspp = 0.01;

            if ((spp + dspp) > send){
                dspp = send-spp;
            }
            
            ux = spp+dspp;
            lx = spp;

            uxm = spp+dspp;
            lxm = spp;

            if (ux > s.back()) {
            ind = s.size() - 1;
            } 
            
            else {
                int lower_bound = 0;
                int upper_bound = s.size() - 1;
                while (lower_bound < upper_bound) {
                    int mid = (lower_bound + upper_bound) / 2;
                    if (s[mid] < ux) {
                        lower_bound = mid + 1;
                    } 
                    else {
                        upper_bound = mid;
                    }
                }
                ind = upper_bound;
            }

            if (secondbeamON){
                for (int kd = 0; kd < gcodeinds.size(); kd++){
                    int indm = 1;
                    int indexval = gcodeinds[kd];

                    if (uxm > sm[indexval].back()) {
                        indm = sm[indexval].size() - 1;
                    } 
                    else {
                        int lower_bound = 0;
                        int upper_bound = sm[indexval].size() - 1;
                        while (lower_bound < upper_bound) {
                            int mid = (lower_bound + upper_bound) / 2;
                            if (sm[indexval][mid] < uxm) {
                                lower_bound = mid + 1;
                            } 
                            else {
                                upper_bound = mid;
                            }
                        }
                        indm = upper_bound;
                    }
                    indms[indexval] = indm;
                }
            }
            

            if (!(type == "static" && secondbeamON) || !(type == "points" && secondbeamON)){

                if (ux > meltedtime[meltedtime.size()-1]){
                    indMP = meltedtime.size()-1;
                }

		        else if (ux > meltedtime[indMP+1]){
			        while (ux > meltedtime[indMP]){
                    indMP+=1;
			        }
		        }

		        else if (ux > meltedtime[indMP]){
		            indMP +=1;
		        }

                //power reducer for radiation from previous timestep
                pmult = ( ((p[ind])/1000) - (rad_pow[indMP]))/ (p[ind]/1000); //
                //ensure pmult is within 0 and 1
                pmult = (pmult > 1) ? 1 : (pmult < 0 || isnan(pmult)) ? 0 : pmult;
                if (secondbeamON){
                    for (int kd = 0; kd < gcodeinds.size(); kd++){
                        int indexval = gcodeinds[kd];
                        int indm = indms[indexval];
                        double pmult_ = ( ((pm[indexval][indm])/1000) - (rad_powm[indexval][indMP]))/ (pm[indexval][indm]/1000);
                        pmult_ = (pmult_ > 1) ? 1 : (pmult_ < 0 || isnan(pmult_)) ? 0 : pmult_;
                        pmultm[indexval] = pmult_;
                    }
                }
            }
            else {
                pmult = 1;
            }

            spp+=dspp;

            // G-L integration parameters
            for (int i = 0; i < roots[gauss_order-1].size(); i ++){
                int_points.push_back(((ux-lx)/2)*roots[gauss_order-1][i] + (ux+lx)/2); //integration points for G-L integration
                weights_tot.push_back(((ux-lx)/2)*weights[gauss_order-1][i]); //weights for G-L integration
                spoints.push_back(s[ind-1]); //time history
                posxpoints.push_back(x[ind-1]); //x location history
                posypoints.push_back(y[ind-1]); //y location history
                poszpoints.push_back(z[ind-1]); //z location history
                vxpoints.push_back(vx[ind]); //x velocity history
                vypoints.push_back(vy[ind]); //y velocity history
                Ppoints.push_back(p[ind]); //beam power history
                p_points.push_back(pmult); //radiation loss history
            
                
                if (secondbeamON){
                    for (int kd = 0; kd < gcodeinds.size(); kd++){
                        int indexval = gcodeinds[kd];
                        int indm = indms[indexval];

                        int_pointsm[indexval].push_back(((uxm-lxm)/2)*roots[gauss_order-1][i] + (uxm+lxm)/2);
                        weights_totm[indexval].push_back(((uxm-lxm)/2)*weights[gauss_order-1][i]);
                        spointsm[indexval].push_back(sm[indexval][indm-1]); //time history
                        posxpointsm[indexval].push_back(xm[indexval][indm-1]); //x location history
                        posypointsm[indexval].push_back(ym[indexval][indm-1]); //y location history
                        poszpointsm[indexval].push_back(zm[indexval][indm-1]); //z location history
                        vxpointsm[indexval].push_back(vxm[indexval][indm]); //x velocity history
                        vypointsm[indexval].push_back(vym[indexval][indm]); //y velocity history
                        Ppointsm[indexval].push_back(pm[indexval][indm]); //beam power history
                        p_pointsm[indexval].push_back(pmultm[indexval]); //radiation loss history
                    }
                }
            }
        }

        // get current beam position 
        vector<int> indms;
        indms.resize(5);
        vector<vector<double>> radsm;
        radsm.resize(5);

        vector<double> rads;
        
        double Xpos,Ypos,Zpos;
        ind = 1; 
        if (send > s[s.size()-1]){
            ind = s.size()-1;
        }

		else if (send > s[ind+1]){
			while (send > s[ind]){
                ind+=1;
			}
		}

	    else if (send > s[ind]){
		    ind +=1;
		}

        Xpos = (x[ind-1]) +(vx[ind])*((send-s[ind-1]));
        Ypos = (y[ind-1]) + (vy[ind])*((send-s[ind-1]));
        Zpos = (z[ind-1]);
        //Zpos = fabs(z[ind-1]-BCvals[4]*2);
        MPcoords.push_back({Xpos,Ypos,Zpos});

        if (type == "meltpool" or type == "dynamic" || type == "dynamic-points") {

            vector <double> xdynam;
            vector <double> ydynam;
            vector<double> zdynam;
            MP_nodes.clear();
            nodeinds.clear();

            if (type != "meltpool"){
                if (domain[2][0] != domain[2][domain[2].size()-1]){
                    modelMPdepth = true;
                }
            }

            dyn_mesh(xdynam,ydynam,zdynam,Xpos,Ypos,Zpos,send,meltednodes,meltedtime,nodelengths,rads,MPparams,modelMPdepth);

            int nodenum = 1;
            for (int w = 0; w < xdynam.size(); w++){
                double x = xdynam[w];
                double y = ydynam[w];
                double z = zdynam[w];
                MP_nodes.push_back(make_node(nodenum,x,y,z,baseT,0));
                nodenum+=1;
            }

            nodeinds.push_back(xdynam.size());

            xdynam.clear();
            ydynam.clear();
            zdynam.clear();

            if (secondbeamON){
                nodeindsm.clear();
                radsm.clear();
                int nodesum = nodeinds[0];
                for (int kd = 0; kd < gcodeinds.size(); kd++){
                    int indexval = gcodeinds[kd];

                    double Xposm,Yposm,Zposm;
                    ind = 1; 

                    if (send > sm[indexval][sm[indexval].size()-1]){
                        ind = sm[indexval].size()-1;
                    }

		            else if (send > sm[indexval][ind+1]){
			            while (send > sm[indexval][ind]){
                            ind+=1;
			            }
		            }

	                else if (send > sm[indexval][ind]){
		                ind +=1;
		            }

                    Xposm = (xm[indexval][ind-1]) +(vxm[indexval][ind])*((send-sm[indexval][ind-1]));
                    Yposm = (ym[indexval][ind-1]) + (vym[indexval][ind])*((send-sm[indexval][ind-1]));
                    Zposm = (zm[indexval][ind-1]);
                    indms[indexval] = ind;

                    xdynam.clear();
                    ydynam.clear();
                    zdynam.clear();
                    
                    dyn_mesh(xdynam,ydynam,zdynam,Xposm,Yposm,Zposm,send,meltednodesm[indexval],meltedtime,nodelengthsm[indexval],radsm[indexval],MPparams,modelMPdepth);
                    
                    nodesum+=xdynam.size();

                    nodeindsm[indexval] = xdynam.size();

                    for (int w = 0; w < xdynam.size(); w++){
                        double x = xdynam[w];
                        double y = ydynam[w];
                        double z = zdynam[w];
                        MP_nodes.push_back(make_node(nodenum,x,y,z,baseT,0));
                        nodenum+=1;
                    }


                    MPcoordsm[indexval].push_back({Xposm,Yposm,Zposm});
                    
                }
            }

        }
        else{
            rads.push_back(0);
            rads.push_back(fabs(domain[0][domain[0].size()-1]- domain[0][0])/(domain[0].size()-1));
            rads.push_back(fabs(domain[1][domain[1].size()-1]- domain[1][0])/(domain[1].size()-1));

        }

        //melt pool nodes integration (constant properties above liquidus)
        if (type == "meltpool" or type == "dynamic" || type == "dynamic-points") {
            #pragma omp parallel for num_threads(nthreads)
            for (int iter = 0; iter < MP_nodes.size(); iter++)
            {
                double T = 0;
                double Tn = 0;
	            double X, Y, Z;
	            X  = MP_nodes[iter].x;
                Y = MP_nodes[iter].y;
	            Z = MP_nodes[iter].z;

                for (int i = 0; i < weights_tot.size(); i++)
                {
 		            T+=weights_tot[i]*integrand_dim(prefactor,Ppoints[i],p_points[i],int_points[i],X,Y,Z,send,a,sig,vxpoints[i],vypoints[i], posxpoints[i], posypoints[i], poszpoints[i], spoints[i],modelBC,BCvals,delta);
                    if (secondbeamON){
                        //#pragma omp parallel for num_threads(2)
                        for (int kd = 0; kd < gcodeinds.size(); kd++){
                            int indexval = gcodeinds[kd];
 		                    T+=weights_totm[indexval][i]*integrand_dim(prefactor,Ppointsm[indexval][i],p_pointsm[indexval][i],int_pointsm[indexval][i],X,Y,Z,send,a,sig,vxpointsm[indexval][i],vypointsm[indexval][i], posxpointsm[indexval][i], posypointsm[indexval][i], poszpointsm[indexval][i], spointsm[indexval][i],modelBC,BCvals,delta);
                        }
                    }
                }

                if (corr){
                    double incTs = 0;
                    #pragma omp parallel for num_threads(2)
                    for (int ind = startidd; ind < idd; ind++){
                        double t_start = start_times1[ind];
                        double t_end = end_times1[ind];
                        double zh = Z - layer_heights1[ind];
                                
                        if (zh == 0){
                            zh = 1e-6;
                        }
                                
                        double b = std::pow(zh, 2) / (4 * a);
                        //double preff = 1 / (std::pow(12 * a, 1.5) * std::sqrt(b));

                        double preff = -eta/(2*pow(pi,3/2)*rho*cp*a*zh);


                        incTs += preff*(layer_pows1[ind]-avgrads1[ind])*compute_integral_with_incomplete_gamma(a,zh,send,t_start,t_end);
                        
                        if (secondbeamON){

                            for (int kd = 0; kd < gcodeinds.size(); kd++){
                                int indexval = gcodeinds[kd];
                                #pragma omp parallel for num_threads(2)
                                for (int ind = startidd; ind < idd; ind++){
                                    double t_start = start_times[indexval][ind];
                                    double t_end = end_times[indexval][ind];
                                    double zh = Z - layer_heights[indexval][ind];
                                
                                    if (zh == 0){
                                        zh = 1e-6;
                                    }
                                
                                    incTs += preff*(layer_pows[indexval][ind]-avgrads[indexval][ind])*compute_integral_with_incomplete_gamma(a,zh,send,t_start,t_end);
                                }
                            }
                        }
                        

                        if (modelBC){
                            double Zbot = BCvals[4];
                            zh = Z + layer_heights1[ind] - Zbot*2;
                                
                            if (zh == 0){
                                zh = 1e-6;
                            }
                                
                            b = std::pow(zh, 2) / (4 * a);

                            preff = 1 / (std::pow(12 * a, 1.5) * std::sqrt(b));

                            incTs += preff*(layer_pows1[ind]-avgrads1[ind])*compute_integral_with_incomplete_gamma(a,zh,send,t_start,t_end);

                            
                            if (secondbeamON){
                                for (int kd = 0; kd < gcodeinds.size(); kd++){
                                    int indexval = gcodeinds[kd];
                                    double Zbot = BCvals[4];
                                    zh = Z + layer_heights[indexval][ind] - Zbot*2;
                                
                                    if (zh == 0){
                                        zh = 1e-6;
                                    }
                                
                                    b = std::pow(zh, 2) / (4 * a);

                                    preff = 1 / (std::pow(12 * a, 1.5) * std::sqrt(b));

                                    incTs += preff*(layer_pows[indexval][ind]-avgrads[indexval][ind])*compute_integral_with_incomplete_gamma(a,zh,send,t_start,t_end);
                                }
            
                            }
                            
                        }                                
                    } 

                    Tn = T + prefactor*incTs + baseT;
                    MP_nodes[iter].T = Tn;
                }
                else{
                    Tn = T + baseT;
                    MP_nodes[iter].T = Tn;
                }
            
                if (MP_nodes[iter].MP != 1 && Tn>= meltT){
                    MP_nodes[iter].MP = 1;
                }
            }

        }

        //longer range nodes integration (variable properties below liquidus)
        if (type != "meltpool" || ifSP) {
            #pragma omp parallel for num_threads(nthreads)
            for (int iter = 0; iter < course_nodes.size(); iter++)
            {
                double T = 0;
                double Tn = 0;
	            double X  = course_nodes[iter].x;
                double Y = course_nodes[iter].y;
	            double Z = course_nodes[iter].z;

                //initialize properties and prefactor for variable properties
                double prefd = prefactor;
                double ad = a;
                double cpd = cp;
                double kd = k;
            
                //check if prior temperature history
                if (Tatdt.size() > 0){
                    double Tcurr = Tatdt_course[Tatdt_course.size()-1][iter];

                    // get node property value based on previous temperature
                    kd = 0;
                    for (int i = 0; i < k_coefficients.size(); i++){
                        kd += k_coefficients[i]*pow(Tcurr,i);
                    }
                    cpd = 0;
                    for (int i = 0; i < c_coefficients.size(); i++){
                        cpd += c_coefficients[i]*pow(Tcurr,i);
                    
                    }

                    // ensure properties are bound at edges of polynomial
                    kd = (isnan(kd) || kd > k) ? k : (kd < k_T0) ? k_T0 : kd;
                    cpd = (isnan(cpd) || cpd > cp) ? cp : (cpd < cp_T0) ? cp_T0 : cpd;

                    prefd = 2*eta/(rho*cpd*pow((pi/3),1.5));
                    ad = kd/(rho * cpd);

                }

                // start with weighted average properties when no prior temperature history
                else{
                    prefd = prefactor;
                    ad = a;
                    cpd = cp;
                    kd = k;
                }

                // perform integral
                #pragma omp parallel for num_threads(nthreads)
                for (int i = 0; i < weights_tot.size(); i++)
                {
 		            T+=weights_tot[i]*integrand_dim(prefd,Ppoints[i],p_points[i],int_points[i],X,Y,Z,send,ad,sig,vxpoints[i],vypoints[i], posxpoints[i], posypoints[i], poszpoints[i], spoints[i],modelBC,BCvals,delta);
                    
                    if (secondbeamON){
                        for (int kd = 0; kd < gcodeinds.size(); kd++){
                            int indexval = gcodeinds[kd];
                            T+=weights_totm[indexval][i]*integrand_dim(prefd,Ppointsm[indexval][i],p_pointsm[indexval][i],int_pointsm[indexval][i],X,Y,Z,send,ad,sig,vxpointsm[indexval][i],vypointsm[indexval][i], posxpointsm[indexval][i], posypointsm[indexval][i], poszpointsm[indexval][i], spointsm[indexval][i],modelBC,BCvals,delta);
                        }
                    }
                }

                if (corr){
                    double incTs = 0;
                    #pragma omp parallel for num_threads(2)
                    for (int ind = startidd; ind < idd; ind++){
                        double t_start = start_times1[ind];
                        double t_end = end_times1[ind];
                        double zh = Z - layer_heights1[ind];
                                
                        if (zh == 0){
                            zh = 1e-6;
                        }
                                
                        double b = std::pow(zh, 2) / (4 * ad);
                        //double preff = 1 / (std::pow(12 * ad, 1.5) * std::sqrt(b));
                        double preff = -eta/(2*pow(pi,3/2)*rho*cpd*ad*zh);

                        incTs += preff*(layer_pows1[ind]-avgrads1[ind])*compute_integral_with_incomplete_gamma(ad,zh,send,t_start,t_end);

                        //cout << incTs << endl;
                        
                        if (secondbeamON){

                            for (int kd = 0; kd < gcodeinds.size(); kd++){
                                int indexval = gcodeinds[kd];
                                #pragma omp parallel for num_threads(2)
                                for (int ind = startidd; ind < idd; ind++){
                                    double t_start = start_times[indexval][ind];
                                    double t_end = end_times[indexval][ind];
                                    double zh = Z - layer_heights[indexval][ind];
                                
                                    if (zh == 0){
                                        zh = 1e-6;
                                    }

                                    incTs += preff*(layer_pows[indexval][ind]-avgrads[indexval][ind])*compute_integral_with_incomplete_gamma(ad,zh,send,t_start,t_end);
                                }
                            }
                        }
                        
                        if (modelBC){
                            double Zbot = BCvals[4];
                            zh = Z + layer_heights1[ind] - Zbot*2;
                                
                            if (zh == 0){
                                zh = 1e-6;
                            }
                                
                            b = std::pow(zh, 2) / (4 * ad);

                            preff = 1 / (std::pow(12 * ad, 1.5) * std::sqrt(b));

                            incTs += preff*(layer_pows1[ind]-avgrads1[ind])*compute_integral_with_incomplete_gamma(ad,zh,send,t_start,t_end);

                            
                            if (secondbeamON){
                                for (int kd = 0; kd < gcodeinds.size(); kd++){
                                    int indexval = gcodeinds[kd];
                                    double Zbot = BCvals[4];
                                    zh = Z + layer_heights[indexval][ind] - Zbot*2;
                                
                                    if (zh == 0){
                                        zh = 1e-6;
                                    }
                                
                                    b = std::pow(zh, 2) / (4 * ad);

                                    preff = 1 / (std::pow(12 * ad, 1.5) * std::sqrt(b));

                                    incTs += preff*(layer_pows[indexval][ind]-avgrads[indexval][ind])*compute_integral_with_incomplete_gamma(ad,zh,send,t_start,t_end);
                                }
            
                            }
                            
                        }                                
                    } 

                    Tn = T + prefd*incTs + baseT;
                    course_nodes[iter].T = Tn;
                }
                else{
                    Tn = T + baseT;
                    course_nodes[iter].T = Tn;
                }
            
                if (course_nodes[iter].MP != 1 && Tn>= meltT){
                    course_nodes[iter].MP = 1;
                }
            } 
        }      

        vector <double> tempx, tempy, tempz, tempT;
        bool myBool = false;
        vector <vector <double>> tempxm, tempym, tempzm, tempTm;
        tempxm.resize(5);
        tempym.resize(5);
        tempzm.resize(5);
        tempTm.resize(5);

        if (type == "meltpool" or type == "dynamic" || type == "dynamic-points")  
        {
            // Process MP_nodes (before adding nodesm)
            for (int iter = 0; iter < nodeinds[0]; iter++) {  // Everything before first start index

                if (MP_nodes[iter].T >= meltT) {
                    currT.push_back(MP_nodes[iter].T);
                    melted.push_back(MP_nodes[iter].MP);
                    currx.push_back(MP_nodes[iter].x);
                    curry.push_back(MP_nodes[iter].y);
                    currz.push_back(MP_nodes[iter].z);
                    tempz.push_back(MP_nodes[iter].z);
                    if (MP_nodes[iter].z == Zpos) {
                        tempT.push_back(MP_nodes[iter].T);
                        tempx.push_back(MP_nodes[iter].x);
                        tempy.push_back(MP_nodes[iter].y);
                    }
                }
            }

            if (secondbeamON) {
                int start_idx = nodeinds[0];
                for (int kd = 0; kd < gcodeinds.size(); kd++){
                    int indexval = gcodeinds[kd];
                    //cout << indexval << endl;
                    int end_idx = start_idx+nodeindsm[indexval];

                    for (int iter = start_idx; iter <= end_idx; iter++) {
                        if (MP_nodes[iter].T >= meltT) {
                            tempzm[indexval].push_back(MP_nodes[iter].z);
                            currT.push_back(MP_nodes[iter].T);
                            melted.push_back(MP_nodes[iter].MP);
                            currx.push_back(MP_nodes[iter].x);
                            curry.push_back(MP_nodes[iter].y);
                            currz.push_back(MP_nodes[iter].z);

                            if (MP_nodes[iter].z == MPcoordsm[indexval][MPcoordsm[indexval].size()-1][2]) { // Use last z in section
                                tempTm[indexval].push_back(MP_nodes[iter].T);
                                tempxm[indexval].push_back(MP_nodes[iter].x);
                                tempym[indexval].push_back(MP_nodes[iter].y);
                            }
                        }
                    }
                    start_idx+=nodeindsm[indexval];
                }
            }
        }

        else if (type == "static")  {
            for (int iter = 0; iter < course_nodes.size(); iter++ ){
                if (course_nodes[iter].T >= meltT){
                    if (course_nodes[iter].z == Zpos){
                        tempT.push_back(course_nodes[iter].T);
                        tempx.push_back(course_nodes[iter].x);
                        tempy.push_back(course_nodes[iter].y);
                        tempz.push_back(course_nodes[iter].z);
                    }
                }   
            }

        } 

        if (send != 0){
            vector <double >tempmelted;
            tempmelted.reserve(2);
            nodelengths = tempT.size();

            if (nodelengths > 3){
                tempmelted.push_back(nodelengths * rads[1] * rads[2]); //calculate melt pool area
                double maxz = *max_element(tempz.begin(), tempz.end());
                double minz = *min_element(tempz.begin(), tempz.end());
                double maxy = *max_element(tempy.begin(), tempy.end());
                double miny = *min_element(tempy.begin(), tempy.end());
                tempmelted.push_back(fabs(maxz-minz)); //meltpool depth
                meltedwidth.push_back(fabs(maxy-miny));

                double effpow =  p[ind];
                double radloss = 0;
                //calculate radiation iteratively
                if (rad_pow.size() > 1){
                    effpow = p[ind] -rad_pow[rad_pow.size()-1];
                    radloss = rad_pow[rad_pow.size()-1];
                }
                if (effpow < 0){
                    effpow = 0;
                }

                double rad = iter_radiation(tempT,rads[1],rads[2],emis,ambient_temp,baseT,eta,p[ind], radloss);
                //rad = 0;
                //cout << endl << "beam 1 rad: " << rad << endl;
                rad_pow.push_back(rad);
            }
            else{
                tempmelted.push_back(0);
                tempmelted.push_back(0);
                meltedwidth.push_back(0);
                rad_pow.push_back(0);
            }

            meltednodes.push_back(tempmelted);

            if (secondbeamON){
                vector <double>tempmeltedm;
                tempmeltedm.reserve(2);
                for (int kd = 0; kd < gcodeinds.size(); kd++){
                    int indexval = gcodeinds[kd];
                    tempmeltedm.clear();
                    nodelengthsm[indexval] = tempTm[indexval].size();
                    if (nodelengthsm[indexval] > 3){
                        tempmeltedm.push_back(nodelengthsm[indexval] * radsm[indexval][1] * radsm[indexval][2]); //calculate melt pool area
                        double maxz = *max_element(tempzm[indexval].begin(), tempzm[indexval].end());
                        double minz = *min_element(tempzm[indexval].begin(), tempzm[indexval].end());
                        double maxy = *max_element(tempym[indexval].begin(), tempym[indexval].end());
                        double miny = *min_element(tempym[indexval].begin(), tempym[indexval].end());
                        tempmeltedm.push_back(fabs(maxz-minz)); //meltpool depth
                        //meltedwidthm[indexval].push_back(fabs(maxy-miny));

                        double effpow =  pm[indexval][indms[indexval]];
                        double radloss = 0;
                        //calculate radiation iteratively
                        if (rad_powm[indexval].size() > 1){
                            effpow = pm[indexval][indms[indexval]] -rad_powm[indexval][rad_powm[indexval].size()-1];
                            radloss = rad_powm[indexval][rad_powm[indexval].size()-1];
                        }
                        if (effpow < 0){
                            effpow = 0;
                        }

                        double rad = iter_radiation(tempTm[indexval],radsm[indexval][1],radsm[indexval][2],emis,ambient_temp,baseT,eta,pm[indexval][indms[indexval]], radloss);
                        //cout << "beam "<< indexval +2 <<" rad: " << rad << endl;
                        rad_powm[indexval].push_back(rad);
                    }
                    else{
                        tempmeltedm.push_back(0);
                        tempmeltedm.push_back(0);
                        //meltedwidth.push_back(0);
                        rad_powm[indexval].push_back(0);
                    }
                    meltednodesm[indexval].push_back(tempmeltedm);
                }  
            }
        }

        else{
            meltednodes.push_back({0,0});
        }
        if (send != 0){
            meltedtime.push_back(send);
        }

        vector <double> courseT = {};

        if (type != "meltpool" || ifSP){
            for (int iter = 0; iter < course_nodes.size(); iter++ ){
                courseT.push_back(course_nodes[iter].T);
                if (type == "dynamic" || type == "dynamic-points"){
                    if (course_nodes[iter].T <= meltT){
                        currT.push_back(course_nodes[iter].T);
                        melted.push_back(course_nodes[iter].MP);
                        currx.push_back(course_nodes[iter].x);
                        curry.push_back(course_nodes[iter].y);
                        currz.push_back(course_nodes[iter].z);
                    }
                }
                else if (type == "static" or type == "points" || ifSP){
                    currT.push_back(course_nodes[iter].T);
                    melted.push_back(course_nodes[iter].MP);
                    currx.push_back(course_nodes[iter].x);
                    curry.push_back(course_nodes[iter].y);
                    currz.push_back(course_nodes[iter].z);

                    if (course_nodes[iter].x == SPs[0][0] && course_nodes[iter].y == SPs[1][0] && course_nodes[iter].z == SPs[2][0]){
                        if (!myBool){
                            myBool = true;
                            probe1.push_back(course_nodes[iter].T);
                            //cout << course_nodes[iter].T << endl;
                        }
                    }
                }
	        }
        }

        Tatdt_course.push_back(courseT);
        Tatdt.push_back(currT);
	    MPatdt.push_back(melted);
        xatdt.push_back(currx);
        yatdt.push_back(curry);
        zatdt.push_back(currz);
 
        currT.clear();
        melted.clear();
        currx.clear();
        curry.clear();
        currz.clear();

        int i = Tatdt.size()-1;
        
        myfile << std::fixed << std::setprecision(4) << meltedtime[i] << ',' << meltednodes[i][0] << ',' << meltednodes[i][1] << ',' << MPcoords[i][0] << ',' << MPcoords[i][1] << ',' << MPcoords[i][2] << endl;

        send += time_freq; //increment timestep
    }
    myfile.close();

    for (int kd = 0; kd < gcodeinds.size(); kd++){
        int indexval = gcodeinds[kd];
        stringstream ss;
        ss<<indexval+2;  
        string s;  
        ss>>s;  
        string filename = "MP";
        string filenum = s;
        string extension =  ".csv";
        string rtn1 = datapath +'/' + filename + filenum + extension;
        ofstream myfile1(rtn1);
        for (int i = 0 ; i < meltedtime.size() ; i ++){
            myfile1 << std::fixed << std::setprecision(4) << meltedtime[i] << ',' << meltednodesm[indexval][i][0] << ',' << meltednodesm[indexval][i][1] << ',' << MPcoordsm[indexval][i][0] << ',' << MPcoordsm[indexval][i][1] << ',' << MPcoordsm[indexval][i][2] << endl;
        }
        myfile1.close();
    }

    for (int i = 0 ; i < Tatdt.size() ; i ++){

        if (i%everyNframes == 0){
            stringstream ss;
            ss<<i;  
            string s;  
            ss>>s;  
            string filename = "snap.";
            string filenum = s;
            string extension =  ".csv";
            string rtn1 = datapath +'/' + filename + filenum + extension;
            ofstream myfile1(rtn1);
            myfile1 << meltedtime[i] << endl;
            int vsize1 = Tatdt[i].size();
            for (int n=0; n<vsize1; n++)
            {   
                myfile1 << Tatdt[i].at(n) <<','<< MPatdt[i].at(n) << ','<< xatdt[i].at(n) << ','<< yatdt[i].at(n)<< ','<< zatdt[i].at(n)<<endl;

            }
            myfile1.close();

        }

        
    }



}
