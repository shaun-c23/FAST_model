# -*- coding: utf-8 -*-
"""
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
"""

'''both this file and 'G1_code_EBAM' need to be in the same folder
   -> gcode_interface_EBAM.py is used to update path, gcode files and get final results
   -> G1_code_EBAM is the library to get this data. Do not modify this 
      file unless you know whatever it is that needs to be fixed
   -> contact Shaun Cooke if any issues that you can't resolve (src@student.ubc.ca)'''

'''add path to files if necessary:
    Spyder->tools-> PYTHONPATH manager -> add path to folderwhere these files are'''

from gcode_library import get_laser_coords
from gcode_library import param_eqn

'''
Gcode Summary

GCODE COMMANDS FOR THIS CODE ARE ABSOLUTE, NOT INCREMENTAL 
(always uses global coords, not relative to previous position)

G0 - rapid traverse (moving while not printing)
G1 - linear feed (printing along straight line from X0 Y0 Z0 to X1 Y1 Z1)
F after G1 - feed rate in mm/min
G4 - dwell time
F after G4 - time to stay at previous position in seconds
L - Beam current
M3 - laser on
M5 - laser off
'''

''' beam voltage '''
beamV = 100 #kV

#path_gcode can be in .txt, .csv, or .gcode format

path_folder = "gcodefiles/"
path_gcode = path_folder+"triangle_gcode.txt"

'''LEAM can move up to 10cm in 250 us --> 400000 mm/s --> 24e6 mm/min'''
travspeed = 24e6 #G0 travel speed in mm/min

gcode_coords = get_laser_coords(path_gcode,travspeed)
print(gcode_coords)
start_line = 0 #line to start reading gcode_coords at

'''LEAM can move up to 10cm in 250 us --> 250us can be assumed to be control speed delay'''
dt = 250e-6

laser_position = param_eqn(gcode_coords,start_line,beamV,dt)
laser_position.to_csv(path_folder+"laser_pos.txt",index = False, header = None) #path folder and name of processed gcode file to be created
