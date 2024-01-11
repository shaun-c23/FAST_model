ABOUT THIS CODE PACKAGE

This code is written in c++ and is best run using g++ compiler with the following dependencies:
- Apple clang version 15.0.0 (clang-1500.1.0.2.5)
- enusre openmp is installed
- running version c++17

There are a number of files located in this package and their descriptions are below:
	- mat_props.h: includes the material properties
	- poly_fit.h: curve fit for weighted average properties
	- domain.h: sets up the domain based on input
	- dynamic_mesh.h: dynamic meshing functions
	- gauss_quad.h: gaussian-quadrature integration
	- mesh.h: contains functions for creating “mesh”
	- radiation.h: iterative radiation losss functions
	- read_laserpos.h: reads the processed gcode file
	- runcode.h: core code which computes the thermal history
	- main.cpp: main code to run the program and define variables and calculation points

There are also a number of Python file for pre- and post-processing:
	- gcode_interface.py: used to process gcode
	- gcode_library.py: libraries to process gcode
	- process.py: post processing for data

HOW TO RUN CODE

It is important that the processed gcode file is created with the attached Python program: gcode_interface.py. See the “about_gcode.txt” file for more information on this.

Beam information, material, initial temperature, gcodefile, etc., can all be changed in main.cpp file before running.

To run the code, assuming the dependencies are installed, the executable can be created with the following on unix systems or equivalently:

g++ -std=c++17 -Xpreprocessor -fopenmp main.cpp -o main -lomp

After the executable is created, the program can be run with following inputs:

./main

OUTPUT DATA

The results are saved in .csv format which can be adjusted in the run code.h file at the end of the code to another format. The format of the .csv files is as follows:

time1
T1, melted1, x1, y1, z1
T2, melted2, x2, y2, z2
.
.
.
.
Tn, meltedn, xn, yn, zn

Where time1 is the tilmestep for the csv file containing the calculation points for the domain, T is the temperature in ˚C, melted is if the node has exceed the solidus point (1 for melted, 0 for unmelted) and, x,y, and z, are the spatial coordinates in mm.


DATA VISUALIZATION

Software such as paraview or a Python program can be used to visualize the data. Please see example Python code (readData.py) to visualize the data or convert to a file readable for paraview.

