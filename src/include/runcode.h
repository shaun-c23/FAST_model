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

    string filename = "MParea";
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

    if (ifSP || type == "points"){
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

    vector< vector<double> > meltednodes;
    vector< double > meltedtime;
    vector< double > meltedwidth;
    meltednodes.push_back({0,0});
    meltedtime.push_back(0);
    meltedwidth.push_back(0);

    double nodelengths = 0;
    int updateInterval = 2;

    vector< vector<double> > MPcoords;
    MPcoords.push_back({0,0,0});

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
            

            if (type != "points"){

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
                        p_pointsm[indexval].push_back(1); //radiation loss history
                    }
                }
            }
        }

        // get current beam position 
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

        if (type == "meltpool" or type == "dynamic") {

            vector <double> xdynam;
            vector <double> ydynam;
            vector<double> zdynam;
            nodes.clear();
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
                nodes.push_back(make_node(nodenum,x,y,z,baseT,0));
                nodenum+=1;
            }

            xdynam.clear();
            ydynam.clear();
            zdynam.clear();

        }
        else{
            rads.push_back(0);
            rads.push_back(fabs(domain[0][domain[0].size()-1]- domain[0][0])/(domain[0].size()-1));
            rads.push_back(fabs(domain[1][domain[1].size()-1]- domain[1][0])/(domain[1].size()-1));

        }

        //melt pool nodes integration (constant properties above liquidus)
        if (type == "meltpool" or type == "dynamic") {
            #pragma omp parallel for num_threads(nthreads)
            for (int iter = 0; iter < nodes.size(); iter++)
            {
                double T = 0;
	            double X, Y, Z;
	            X  = nodes[iter].x;
                Y = nodes[iter].y;
	            Z = nodes[iter].z;

                for (int i = 0; i < weights_tot.size(); i++)
                {
 		            T+=weights_tot[i]*integrand_dim(prefactor,Ppoints[i],p_points[i],int_points[i],X,Y,Z,send,a,sig,vxpoints[i],vypoints[i], posxpoints[i], posypoints[i], poszpoints[i], spoints[i],modelBC,BCvals,delta);
                }

                double Tn = T + baseT;
                nodes[iter].T = Tn;
            
                if (nodes[iter].MP != 1 && Tn>= meltT){
                    nodes[iter].MP = 1;
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

        if (type == "meltpool" or type == "dynamic")  {
            for (int iter = 0; iter < nodes.size(); iter++ ){
                if (nodes[iter].T >= meltT){
                    if (nodes[iter].z == Zpos){
                        tempT.push_back(nodes[iter].T);
                        tempx.push_back(nodes[iter].x);
                        tempy.push_back(nodes[iter].y);
                        tempz.push_back(nodes[iter].z);
                    }
                }
                if (nodes[iter].T >= meltT){
                    currT.push_back(nodes[iter].T);
                    melted.push_back(nodes[iter].MP);
                    currx.push_back(nodes[iter].x);
                    curry.push_back(nodes[iter].y);
                    currz.push_back(nodes[iter].z);
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

        if (type != "points" && send != 0){
            vector <double >tempmelted;
            tempmelted.reserve(2);
            nodelengths = tempT.size();

            if (nodelengths > 3){
                tempmelted.push_back(nodelengths * rads[1] * rads[2]); //calculate melt pool area
                double maxz = *max_element(currz.begin(), currz.end());
                double minz = *min_element(currz.begin(), currz.end());
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
                if (type == "dynamic"){
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


    for (int i = 0 ; i < Tatdt.size() ; i ++){

        /*
        if (type != "points"){
            myfile << meltedtime[i] << ',' << meltednodes[i][0] << ',' << meltednodes[i][1] << ',' << MPcoords[i][0] << ',' << MPcoords[i][1] << ',' << MPcoords[i][2] << endl;
        }
        */

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
