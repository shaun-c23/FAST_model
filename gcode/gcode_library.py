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

import re
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

#---------------------read Gcode----------------------------------------------
#-------laser coordinate function definition to call from anywhere-------------

def bare_numpy_mat(mat_v, mat_u):
   return np.sqrt(np.sum((mat_v - mat_u) ** 2))

def get_laser_coords(path,travspeed):
  
    gcode_data = {
        'G': [],
        'X': [],
        'Y': [],
        'L': [],
        'F': [],
        'M': []
                 }

    ind = 0
    syms = ['(','G', 'X', 'Y', 'L', 'F', 'M']
    with open(path) as gcode:
        for line in gcode:
            ind+=1
            line = line.strip()
            #print(line)
            for command in syms:
                c_ = line.find(command)
                if command == '(' and c_ != -1:
                    #(line)
                    #print("comment")
                    break
                elif  command != '(' and c_ == -1:
                    C = '?'
                elif command != '(' and c_ != -1:
                    C = command
                    
                    if line[c_+ 1] == ' ':
                        count = 2
                    else:
                        count = 1

                    try:
                        while (line[c_+count] != ' '):
                            C+=line[c_+count]
                            count+=1
                    except:
                        count = 0
                    
                    #print(C
                if  command != '(':
                    gcode_data[command].append(C)
                        #else:
                            #gcode_data[sym].append('?')   

    g_ind = 0
    for val in gcode_data['G']:
        if val == 'G0':
            break
        g_ind+=1
    
    g = gcode_data['G'][g_ind:]
    x = gcode_data['X'][g_ind:]
    y = gcode_data['Y'][g_ind:]
    l = gcode_data['L'][g_ind:]
    f = gcode_data['F'][g_ind:]
    m = gcode_data['M'][g_ind:]

    index = 0
    for index in range(len(g)):
        if  g[index] == 'G0' and f[index] == '?':
            f[index] = str(travspeed)

        if g[index] == '?' and m[index] != '?':
            g[index] = m[index]
    
        if x[index] == '?':
            x[index] = x[index-1]
        if y[index] == '?':
            y[index] = y[index-1]
        if l[index] == '?':
            l[index] = l[index-1]
        if f[index] == '?':
            f[index] = f[index-1]
    
        x[index] = x[index].strip('X')
        y[index] = y[index].strip('Y')
        l[index] = l[index].strip('L')
        f[index] = f[index].strip('F')

    index = 0
    while index < len(g):
        if  g[index] != 'G0' and g[index] != 'G1' and  g[index] != 'G4' and  g[index] != 'M3' and g[index] != 'M5':
            g.pop(index) 
            x.pop(index)
            y.pop(index) 
            l.pop(index) 
            f.pop(index) 
            m.pop(index) 
            index-=1

        index+=1

    gcode_ = pd.DataFrame(zip(g,x,y,l,f),columns = ['C','X','Y','L','F'])
    gcode_['X'] = gcode_['X'].astype(float)
    gcode_['Y'] = gcode_['Y'].astype(float)
    gcode_['L'] = gcode_['L'].astype(float)
    gcode_['F'] = gcode_['F'].astype(float)
    
    index = 0
    while index < len(gcode_):
        if gcode_['C'][index] == 'G1':
            val = gcode_['F'].loc[index]
            gcode_['F'].loc[index] = val/60
        index+=1

    index = 0
    while index < len(gcode_):
        if gcode_['C'][index] == 'G0':
            val = gcode_['F'].loc[index]
            gcode_['F'].loc[index] = travspeed/60
        index+=1

    index = 0
    while index < len(gcode_):
        if gcode_['C'][index] == 'M3' or gcode_['C'][index] == 'M5':
            val = gcode_['F'].loc[index]
            gcode_['F'].loc[index] = 0
        index += 1
    
    index = 0
    while index < len(gcode_):
        if gcode_['C'][index] == 'M5':
            val = gcode_['L'].loc[index]
            gcode_['L'].loc[index] = 0
        index += 1
        
    return gcode_

def param_eqn(gcode_coords,start_line,beamV,dt):
    
    x=[]
    y=[]
    z=[]
    t=[]
    p=[]
    vxi = []
    vyi = []

    #start with num
    x.append(0)#gcode_coords['X'][start_line])
    y.append(0)#gcode_coords['Y'][start_line])
    z.append(0)
    t.append(0)
    vxi.append(0)
    vyi.append(0)
    p.append(gcode_coords['L'][start_line]*beamV*1000)


    for row_gcode in np.arange(start_line,len(gcode_coords)): #len(gcode_coords)

        if gcode_coords['C'][row_gcode] == 'G0' or gcode_coords['C'][row_gcode] == 'G1':
            if gcode_coords['C'][row_gcode] == 'G1':
                power = gcode_coords['L'][row_gcode]*beamV
                speed = gcode_coords['F'][row_gcode]
            if gcode_coords['C'][row_gcode] == 'G0':
                power = 0
                speed = gcode_coords['F'][row_gcode]

                
            xy = [x[-1],y[-1]]

            d_i = round(bare_numpy_mat(xy,gcode_coords.iloc[row_gcode,1:3]),3)
            d_ix = round(gcode_coords.iloc[row_gcode,1:2][0]-float(xy[0]),3)
            d_iy = round(gcode_coords.iloc[row_gcode,2:3][0]-float(xy[1]),3)
            
            theta = 0
            if d_ix > 0 and d_iy > 0:
                theta = math.atan(d_iy/d_ix)
                
            elif d_ix > 0 and d_iy < 0:
                theta =360*(np.pi/180) + math.atan(d_iy/d_ix)
                
            elif d_ix < 0 and d_iy > 0:
                theta = 180*(np.pi/180) + math.atan(d_iy/d_ix)
                
            elif d_ix < 0 and d_iy < 0:
                theta = 180*(np.pi/180)  + math.atan(d_iy/d_ix)
                
            elif round(d_ix,2) == 0 and d_iy < 0:
                theta = 270*(np.pi/180) 
                
            elif round(d_ix,2) == 0 and d_iy > 0:
                theta = 90*(np.pi/180) 
                
            elif d_ix > 0 and round(d_iy,2)==0:
                theta = 0
                
            elif d_ix < 0 and round(d_iy,2)==0:
                theta = 180*(np.pi/180) 
                
            vx = round(speed*math.cos(theta),4)
                
            vy = round(speed*math.sin(theta),4)
            
            check = True
            if gcode_coords['X'][row_gcode] == x[-1] and gcode_coords['Y'][row_gcode]== y[-1]:
                check = False
            
            if check:

                x.append(gcode_coords['X'][row_gcode])
                y.append(gcode_coords['Y'][row_gcode])
                z.append(0)
                t.append(d_i/speed + t[-1])
                vxi.append(vx)
                vyi.append(vy)
                p.append(power*1000)
        
        elif gcode_coords['C'][row_gcode] == 'G4':

            x.append(gcode_coords['X'][row_gcode])
            y.append(gcode_coords['Y'][row_gcode])
            z.append(0)
            t.append(gcode_coords['F'][row_gcode] + t[-1])
            vxi.append(0)
            vyi.append(0)
            p.append(p[-1])
            
        elif gcode_coords['C'][row_gcode] == 'M3':
            x.append(x[-1])
            y.append(y[-1])
            z.append(0)
            t.append(t[-1]+dt)
            vxi.append(vxi[-1])
            vyi.append(vyi[-1])
            p.append(gcode_coords['L'][row_gcode]*beamV*1000)

        elif gcode_coords['C'][row_gcode] == 'M5':
            x.append(x[-1])
            y.append(y[-1])
            z.append(0)
            t.append(t[-1]+dt)
            vxi.append(vxi[-1])
            vyi.append(vyi[-1])
            p.append(0)

        laser_pos_array = pd.DataFrame(zip(t,x,vxi,y,vyi,z,p),columns = ['t','x','vx','y','vy','z','p'])
        
    return laser_pos_array

