-------------------------------------------------------------------------------------------------------------------------------

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

-------------------------------------------------------------------------------------------------------------------------------

HOW TO RUN CODE

It is important that the processed gcode file is created with the attached Python program: gcode.py. See the “about_gcode.txt” file for more information on this. There are example files located within the gcode>gcodefiles directory.

Beam information, material, initial temperature, gcodefile, etc., can all be changed in the inputs files before running (see about_inputfiles.txt).

**ensure you have version 10 or higher for gcc compiler using c++17 or higher with openmp installed**

Go into the "src" directory and compule "main.cpp" file after your system is properly setup to run this code.

-------------------------------------------------------------------------------------------------------------------------------

Running on Windows:

See also for info on using with Vscode:
https://code.visualstudio.com/docs/cpp/config-mingw

Install Msys2

Then run the following to install the latest gcc compiler:

pacman -S mingw-w64-ucrt-x86_64-gcc

pacman -S mingw-w64-x86_64-openmp

Add to your environment variables:

C:\msys64\ucrt64\bin

C:\msys64\usr\bin

C:\mingw64\bin

compile the code:

g++ -std=c++17 -fopenmp main.cpp -o main

After the executable is created, the program can be run with following inputs:

main.exe

-------------------------------------------------------------------------------------------------------------------------------

Running on Mac:

Ensure you have Xcode and the command line tools installed on your Mac.

If you do not have openmp installed follow these steps:

Install homebrew:

/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

Follow steps from terminal to add home-brew to your PATH

Run:  brew install llvm libomp

compile the code:

g++ -std=c++17 -Xpreprocessor -fopenmp main.cpp -o main -lomp 

If you have problems compiling the code because of -fopenmp or -lomp then try the following:

Update your shell profile file (~/.bash_profile for Bash or ~/.zshrc for Zsh) to include the Homebrew Clang directory at the beginning of the PATH. Add the following line:

nano ~/.bash_profile # or ~/.zshrc for Zsh users

Add the following lines at the end of the file:

export PATH="/usr/local/opt/llvm/bin:$PATH"
export LDFLAGS="-L/usr/local/opt/libomp/lib"
export CPPFLAGS="-I/usr/local/opt/libomp/include"

Save and exit and then run:
source ~/.bash_profile   # or source ~/.zshrc for Zsh users

Now  specify the library path explicitly using -L to compile:

g++ -std=c++17 -Xpreprocessor -fopenmp main.cpp -o main -L/usr/local/opt/libomp/lib -lomp

After the executable is created, the program can be run with the following:

./main
-------------------------------------------------------------------------------------------------------------------------------

Running on Linux:

Similarily, ensure the latest gcc compiler is installed with openmp support.

compile the code:

g++ -std=c++17 -fopenmp main.cpp -o main

Sometimes you may need -lomp, so if issues try adding that and trying again. Follow similar steps to Mac for troubleshooting if you have issues compiling.

After the executable is created, the program can be run with following inputs:

./main

-------------------------------------------------------------------------------------------------------------------------------

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

The meltpool data is also printed for each beam as MP1.csv, MP2.csv, etc.
MP1 (the main beam from BeamPath1) is updated live while the others are generated at the end. They all have the same format:

time1, MP area, MP depth, x, y, z #x,y,z are the position of the beam center
time2, MP area, MP depth, x, y, z
.
.
.
.
timen, MP area, MP depth, x, y, z
-------------------------------------------------------------------------------------------------------------------------------

DATA VISUALIZATION

Software such as paraview or a Python program can be used to visualize the data. Please see example Python code (view_data.py) to visualize the data or you can convert the data use for visualization software such as paraview.

-------------------------------------------------------------------------------------------------------------------------------
