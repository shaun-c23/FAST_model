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
    float x0, x1, f0, f1, g0, e;
	int step = 1, N;
    N = 100;
    e = 1e-5;
    static double Pi = 3.14159265358979323846;

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
