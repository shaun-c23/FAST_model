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
#include <cmath>
#include <array>
#include <map>
#include <any>
#include <regex>

using namespace std;

struct ConfigPath {
    string BeamPath;
    string BeamPath2;
    string BeamPath3;
    string BeamPath4;
    string BeamPath5;
    string BeamPath6;
    string Config;
    string Domain;
    string DataPath;
    string Material;
    string Meltpool;
    bool meltpooldepth;
    string Meshtype;
    bool modelBC;
    string BCPath;
    string SPPath;
    bool ifSP;
    bool secondBeamON;
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
    double everyNframes;
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

struct ConfigSP {
    vector<double> x,y,z;
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
    configPath.BeamPath2 = values[2];
    configPath.BeamPath3 = values[3];
    configPath.BeamPath4 = values[4];
    configPath.BeamPath5 = values[5];
    configPath.BeamPath6 = values[6];
    configPath.modelBC = stringToBool(values[7]);
    configPath.Config = values[8];
    configPath.DataPath = values[9];
    configPath.Domain = values[10];
    configPath.Material = values[11];
    configPath.Meltpool = values[12];
    configPath.meltpooldepth = stringToBool(values[13]);
    configPath.Meshtype = values[14];
    configPath.secondBeamON = stringToBool(values[15]);
    configPath.SPPath = values[16];
    configPath.ifSP = stringToBool(values[17]);
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
    configVars.everyNframes = values[6];
    configVars.Sigma = values[7];
    configVars.TimeStep = values[8];

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

void readStaticPoints(const std::string &FilePath, const std::string &fileName, ConfigSP& config) {
    std::ifstream inputFile(FilePath + fileName);

    if (!inputFile.is_open()) {
        std::cerr << "Error opening the file: " << FilePath + fileName << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        // Ignore comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        if (line.find("x =") != std::string::npos) {
            std::stringstream ss(line.substr(line.find("{") + 1));
            double value;
            while (ss >> value) {
                config.x.push_back(value);
                if (ss.peek() == ',') ss.ignore();
            }
        } else if (line.find("y =") != std::string::npos) {
            std::stringstream ss(line.substr(line.find("{") + 1));
            double value;
            while (ss >> value) {
                config.y.push_back(value);
                if (ss.peek() == ',') ss.ignore();
            }
        } else if (line.find("z =") != std::string::npos) {
            std::stringstream ss(line.substr(line.find("{") + 1));
            double value;
            while (ss >> value) {
                config.z.push_back(value);
                if (ss.peek() == ',') ss.ignore();
            }
        }
    }

    inputFile.close();
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

