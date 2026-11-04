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