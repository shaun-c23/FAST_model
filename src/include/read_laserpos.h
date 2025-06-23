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


int laserpath(string name, vector<double> &t, vector<double> &x,\
    vector<double> &vx,vector<double> &y, vector<double> &vy, vector<double> &z,vector<double> &p)
{
    //cout << "before las" << endl;
    fstream fin;
    fin.open(name); 
    string line;
    vector<string> result;
    string substr;

    // Open an existing file 
    while (!fin.eof( ))      //if not at end of file, continue reading numbers
    {
        result.clear(); //clear vector for every new line
        getline(fin, line); //get line from file
        //cout << "line: " << line << endl;
        stringstream s_stream(line);

        while(s_stream.good()) //sperate delim values while not at end of line
        {
            getline(s_stream,substr,','); //get first string delimited by comma
            result.push_back(substr);
        }

        //append values to end of vector
        try 
        {
            t.push_back(stof(result.at(0)));
            x.push_back(stof(result.at(1)));
            vx.push_back(stof(result.at(2)));
            y.push_back(stof(result.at(3)));
            vy.push_back(stof(result.at(4)));
            z.push_back(stof(result.at(5)));
            p.push_back(stof(result.at(6))*1);
        } 
        
        catch (const std::invalid_argument& e) 
        {

        } 

    }

    fin.close();
    return 0;
}
