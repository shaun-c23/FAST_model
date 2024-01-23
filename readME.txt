---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

ABOUT THIS CODE PACKAGE

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

There are also a number of Python files for pre- and post-processing:
	- gcode_interface.py: used to process gcode
	- gcode_library.py: libraries to process gcode
	- view_data.ipynb: post processing for data

Input files located in the input folder can be chnaged in order to specify the material, beam parameters, trajectory, etc. without having to recompile the code. See the about_inputfiles.txt for more information.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

HOW TO RUN CODE

It is important that the processed gcode file is created with the attached Python program: gcode.py. See the “about_gcode.txt” file for more information on this.

Beam information, material, initial temperature, gcodefile, etc., can all be changed in main.cpp file before running.

**ensure you have version 10 or higher for gcc compiler with openmp installed**

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Running on Windows:

See also for info on using with Vscode:
https://code.visualstudio.com/docs/cpp/config-mingw

Install Msys2

Then run the following to install the latest gcc compiler:

pacman -S mingw-w64-ucrt-x86_64-gcc

pacman -S mingw-w64-x86_64-openmp

Add to your environment variables:

C:\msys64\ucrt64\bin

and

C:\mingw64\bin

compile the code:

g++ -std=c++17 -fopenmp main.cpp -o main

After the executable is created, the program can be run with following inputs:

main.exe

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Running on Mac:

- Apple clang version 15.0.0 (clang-1500.1.0.2.5) (or equivalent with latest gcc compiler) with openmp support

compile the code:

g++ -std=c++17 -Xpreprocessor -fopenmp main.cpp -o main -lomp 

After the executable is created, the program can be run with following inputs:

./main

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Running on Linux:

Similarily, ensure the latest gcc compiler is installed with openmp support.

compile the code:

g++ -std=c++17 -fopenmp main.cpp -o main -lomp 

After the executable is created, the program can be run with following inputs:

./main

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

Where time1 is the tilmestep in s for the csv file containing the calculation points for the domain, T is the temperature in ˚C, melted is if the node has exceed the solidus point (1 for melted, 0 for unmelted) and, x,y, and z, are the spatial coordinates in mm.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

DATA VISUALIZATION

Software such as paraview or a Python program can be used to visualize the data. Please see example Python code (view_data.py) to visualize the data or you can convert the data use for visualization software such as paraview.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
