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
#include <map>
#include <any>
#include <regex>

using namespace std;

struct ConfigPath {
    string BeamPath;
    string Config;
    string Domain;
    string Material;
    string Meltpool;
    bool meltpooldepth;
    string Meshtype;
    bool modelBC;
    string BCPath;
};

struct ConfigVars {
    double AmbientTemp;
    double CoolingTime;
    double Efficiency;
    double InitialTemp;
    double MaxThreads;
    double Sigma;
    double Delta;
    double TimeStep;
};

struct ConfigMPpoints {
    double BoxSizeMult;
    double SurfSpacing;
    double DepthSpacing;
};

struct ConfigDomain {
    double nx,ny,nz;
    double x1,y1,z1;
    double x2,y2,z2;
};

struct ConfigBC {
    double x1,y1,x2,y2,zbot;
};

bool stringToBool(const std::string& str) {
    std::istringstream iss(str);
    bool result;
    iss >> std::boolalpha >> result;
    if (iss.fail()) {
        // Handle conversion failure if needed
        std::cerr << "Conversion failed for string: " << str << std::endl;
    }

    return result;
}


template <typename T>
void printVector(const std::string& name, const std::vector<T>& vec) {
    std::cout << name << ": ";
    for (const auto& element : vec) {
        std::cout << element << " ";
    }
    std::cout << std::endl;
}

int read_input(string fullFilePath, std::map<std::string, double> &configValues)
{
    // Check if the file exists
    if (!std::filesystem::exists(fullFilePath)) {
        std::cerr << "File not found: " << fullFilePath << std::endl;
        return 1;
    }

    // Open the file
    std::ifstream configFile(fullFilePath);

    // Check if the file is open
    if (!configFile.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    // Define a map to store the key-value pairs


    std::string line;
    while (std::getline(configFile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string key;
        char equalSign;
        double value;

        // Read key, equal sign, and value
        if (iss >> key >> equalSign >> value) {
            configValues[key] = value;
        } else {
            std::cerr << "Error parsing line: " << line << std::endl;
        }
    }

    configFile.close();

    return 0;
}

int read_input_files(string fullFilePath, std::map<std::string, std::string> &configValues)
{
// Check if the file exists
    if (!std::filesystem::exists(fullFilePath)) {
        std::cerr << "File not found: " << fullFilePath << std::endl;
        return 1;
    }

    // Open the file
    std::ifstream configFile(fullFilePath);

    // Check if the file is open
    if (!configFile.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    // Define a map to store the key-value pairs
    //std::map<std::string, std::string> configValues;

    std::string line;
    while (std::getline(configFile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        std::istringstream iss(line);
        std::string key, equalSign, value;

        // Read key, equal sign, and value
        if (iss >> key >> equalSign >> value) {
            configValues[key] = value;
        } else {
            std::cerr << "Error parsing line: " << line << std::endl;
        }
    }

    configFile.close();

    return 0;

}


void filenames(string FilePath, string fileName, ConfigPath &configPath){
    std::map<std::string, std::string> configValues;
    read_input_files(FilePath + fileName, configValues);

    vector<string> values;

    for (const auto& pair : configValues) {
        std::string key = pair.first;
        std::string value = pair.second;
        values.push_back(value);
     }
    configPath.BCPath = values[0];
    configPath.BeamPath = values[1];
    configPath.modelBC = stringToBool(values[2]);
    configPath.Config = values[3];
    configPath.Domain = values[4];
    configPath.Material = values[5];
    configPath.Meltpool = values[6];
    configPath.meltpooldepth = stringToBool(values[7]);
    configPath.Meshtype = values[8];

}

void vars(string FilePath, string fileName, ConfigVars &configVars){

    std::map<std::string, double> configValues;
    read_input(FilePath + fileName, configValues);

    vector<double> values;

    for (const auto& pair : configValues) {
        //std::string key = pair.first;
        double value = pair.second;
        values.push_back(value);
     }

    configVars.AmbientTemp = values[0];
    configVars.Delta = values[1];
    configVars.CoolingTime = values[2];
    configVars.Efficiency = values[3];
    configVars.InitialTemp = values[4];
    configVars.MaxThreads = int(values[5]);
    configVars.Sigma = values[6];
    configVars.TimeStep = values[7];
}

void BCs(string FilePath, string fileName, ConfigBC &configBCs, vector<double> &values){

    std::map<std::string, double> configValues;
    read_input(FilePath + fileName, configValues);

    for (const auto& pair : configValues) {
        //std::string key = pair.first;
        double value = pair.second;
        values.push_back(value);
     }
    configBCs.x1 = values[0];
    configBCs.x2 = values[1];
    configBCs.y1 = values[2];
    configBCs.x2 = values[3];
    configBCs.zbot = values[4];

}

void MP_vars(string FilePath, string fileName, ConfigMPpoints &configMP){

    std::map<std::string, double> configValues;
    read_input(FilePath + fileName, configValues);

    vector<double> values;

    for (const auto& pair : configValues) {
        double value = pair.second;
        values.push_back(value);

     }

    configMP.BoxSizeMult = values[0];
    configMP.DepthSpacing = values[1];
    configMP.SurfSpacing= values[2];

}

void get_domain(string FilePath, string fileName, ConfigDomain &configDomain){

    std::map<std::string, double> configValues;
    read_input(FilePath + fileName, configValues);

    vector<double> values;

    for (const auto& pair : configValues) {
        double value = pair.second;
        values.push_back(value);
     }

    configDomain.nx = values[0];
    configDomain.ny = values[1];
    configDomain.nz = values[2];
    configDomain.x1 = values[3];
    configDomain.x2 = values[4];
    configDomain.y1 = values[5];
    configDomain.y2 = values[6];
    configDomain.z1 = values[7];
    configDomain.z2 = values[8];

}

