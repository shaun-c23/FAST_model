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
    Yb = posY + (Vy)*(sp-slim);
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
    
        return prefa*pwr*pmlt*(1/(sqrt(phix*phiy*phiz)))*(exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z-Zb),2)/phiz)
                                                           +exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z+Zb)-Zbot*2,2)/phiz)
                                                           +exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb1),2)/phiy-3*pow((Z-Zb),2)/phiz)
                                                           +exp(-3*pow((X-Xb),2)/phix-3*pow((Y-Yb2),2)/phiy-3*pow((Z-Zb),2)/phiz)
                                                           +exp(-3*pow((X-Xb1),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z-Zb),2)/phiz)
                                                           +exp(-3*pow((X-Xb2),2)/phix-3*pow((Y-Yb),2)/phiy-3*pow((Z-Zb),2)/phiz));

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


void runcode(double baseT, double ambient_temp, double sig, double eta, int nthreads, vector< vector<double> > domain,vector< vector<double> > SPs,vector<string> gcodes, string datapath, double time_freq, double ex_time, string type, string mat,vector<double> MPparams, bool modelMPdepth,bool ifSP,bool modelBC, bool secondbeamON, vector<double> BCvals, double delta, int everyNframes){

    fs::create_directory(datapath);

    string filename = "MP1";
    string extension =  ".csv";
    string rtn = datapath +'/' + filename + extension;
    ofstream myfile(rtn);

    static double pi = 3.14159265358979323846;

    vector<double> x,vx, vy, y, z, s, p;
    if (filesystem::exists(gcodes[0])) {
        laserpath(gcodes[0], s,x,vx,y,vy,z,p); //read gcode file
    } 
    else {
        cout << "Beam path 1 file not found... Program terminating." << endl;
        std::exit(EXIT_FAILURE);
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


    //node spaces
    vector<node> course_nodes;
    vector<node> nodes;
    vector<node> MP_nodes; 
    vector<std::pair<int, int>> index_ranges;
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

    double s1 = 0;
    double s2 = (s[s.size()-1] + ex_time);

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

    double send = 0;
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

        double spp = 0;
        double dspp;
	    int ind = 1;

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
        nodeinds.clear();

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
        MPcoords.push_back({Xpos,Ypos,Zpos});

        if (type == "meltpool" || type == "dynamic" || type == "dynamic-points") {

            vector <double> xdynam;
            vector <double> ydynam;
            vector<double> zdynam;
            nodes.clear();
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
	            double X, Y, Z;
	            X  = MP_nodes[iter].x;
                Y = MP_nodes[iter].y;
	            Z = MP_nodes[iter].z;

                for (int i = 0; i < weights_tot.size(); i++)
                {
 		            T+=weights_tot[i]*integrand_dim(prefactor,Ppoints[i],p_points[i],int_points[i],X,Y,Z,send,a,sig,vxpoints[i],vypoints[i], posxpoints[i], posypoints[i], poszpoints[i], spoints[i],modelBC,BCvals,delta);
                    if (secondbeamON){
                        #pragma omp parallel for num_threads(2)
                        for (int kd = 0; kd < gcodeinds.size(); kd++){
                            int indexval = gcodeinds[kd];
 		                    T+=weights_totm[indexval][i]*integrand_dim(prefactor,Ppointsm[indexval][i],p_pointsm[indexval][i],int_pointsm[indexval][i],X,Y,Z,send,a,sig,vxpointsm[indexval][i],vypointsm[indexval][i], posxpointsm[indexval][i], posypointsm[indexval][i], poszpointsm[indexval][i], spointsm[indexval][i],modelBC,BCvals,delta);
                        }
                    }
                }

                double Tn = T + baseT;
                MP_nodes[iter].T = Tn;
            
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

                //save temperatures
                double Tn = T + baseT;
                course_nodes[iter].T = Tn;
            
                if (course_nodes[iter].MP != 1 && Tn>= meltT){
                    course_nodes[iter].MP = 1;
                }
            } 
        }      

        vector <double> tempx, tempy, tempz, tempT;
        vector <vector <double>> tempxm, tempym, tempzm, tempTm;
        tempxm.resize(5);
        tempym.resize(5);
        tempzm.resize(5);
        tempTm.resize(5);

        

        if (type == "meltpool" || type == "dynamic" || type == "dynamic-points")  
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


        else if (type == "static" || type == "points")  {
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

            if (nodelengths > 3 && !(type == "static" && secondbeamON)){
            //if (nodelengths > 3){
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
        
        myfile << meltedtime[i] << ',' << meltednodes[i][0] << ',' << meltednodes[i][1] << ',' << MPcoords[i][0] << ',' << MPcoords[i][1] << ',' << MPcoords[i][2] << endl;

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
            myfile1 << meltedtime[i] << ',' << meltednodesm[indexval][i][0] << ',' << meltednodesm[indexval][i][1] << ',' << MPcoordsm[indexval][i][0] << ',' << MPcoordsm[indexval][i][1] << ',' << MPcoordsm[indexval][i][2] << endl;
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
