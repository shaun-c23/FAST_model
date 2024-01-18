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



double Pn(double x, int n){
    if(n==0){
        return 1;
    }else if(n==1){
        return x;
    }else{
        return ((2*n-1)*x*Pn(x,n-1)-(n-1)*Pn(x,n-2))/n;
    }
}

double Pn_prime(double x, int n){
    if(n==0){
        return 0;
    }
    if(n==1){
        return 1;
    }
    else{
        return (n/(pow(x,2)-1))*(x*Pn(x,n)-Pn(x,n-1));
    }
}


void findRoots(int n, vector<double> &temp_roots, vector<double> &temp_weights){
    float x0, x1, f0, f1, g0, e, Pi;
	int step = 1, N;
    N = 100;
    e = 1e-5;
    Pi = M_PI;

    for(int i = 1; i<=n; i++)
    {
        x0 = cos(Pi*(i-.25)/(n+.5));
        do
	    {
		    g0 = Pn_prime(x0, n);
		    f0 = Pn(x0, n);
		    if(g0 == 0.0)
		    {
			    cout<<"Mathematical Error.";
			    exit(0);
		    }
		    x1 = x0 - f0/g0;
		    //cout<<"Iteration-"<< step<<":\t x = "<< setw(10)<< x1<<" and f(x) = "<< setw(10)<< Pn(x1,n)<< endl;
		    x0 = x1;

		    step = step+1;

		    if(step > N)
		    {
			    cout<<"Not Convergent.";
			    exit(0);
		    }

		    f1 = Pn(x1, n);

	    }while(fabs(f1)>e);
        temp_roots.push_back(x1);
        temp_weights.push_back(2.0/((1-pow(x1,2))*pow(Pn_prime(x1,n),2)));
    }
}