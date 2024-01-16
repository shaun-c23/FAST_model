#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Shaun

"""

import re
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt

#---------------------read Gcode----------------------------------------------
#-------laser coordinate function definition to call from anywhere-------------
#make input for units so it converts to m only when needed

def bare_numpy_mat(mat_v, mat_u):
   return np.sqrt(np.sum((mat_v - mat_u) ** 2))

'''LEAM can move up to 10cm in 250 us --> 400000 mm/s --> 24e6 mm/min'''
#travspeed = 24e6
def get_laser_coords(path,travspeed):
    #path_folder = "/Users/sc/Desktop/rosen_files/FE_models/qa/layer/gcode/"
    #path = path_folder+"hatch_layer_2.txt"
    gcode_data = {
        'G': [],
        'X': [],
        'Y': [],
        'Z': [],
        'L': [],
        'F': [],
        'M': []
                 }

    ind = 0
    syms = ['(','G', 'X', 'Y','Z', 'L', 'F', 'M']
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
    z = gcode_data['Z'][g_ind:]
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
        if z[index] == '?':
            z[index] = z[index-1]
        if l[index] == '?' and index == 0:
            l[index] = '0'
            
        elif l[index] == '?':
            l[index] = l[index-1]
            
        if f[index] == '?':
            f[index] = f[index-1]
    
        x[index] = x[index].strip('X')
        y[index] = y[index].strip('Y')
        z[index] = z[index].strip('Z')
        l[index] = l[index].strip('L')
        f[index] = f[index].strip('F')

    index = 0
    while index < len(g):
        if  g[index] != 'G0' and g[index] != 'G1' and  g[index] != 'G4' and  g[index] != 'M3' and g[index] != 'M5':
            g.pop(index) 
            x.pop(index)
            y.pop(index) 
            z.pop(index) 
            l.pop(index) 
            f.pop(index) 
            m.pop(index) 
            index-=1

        index+=1

    gcode_ = pd.DataFrame(zip(g,x,y,z,l,f),columns = ['C','X','Y','Z','L','F'])
    gcode_['X'] = gcode_['X'].astype(float)
    gcode_['Y'] = gcode_['Y'].astype(float)
    gcode_['Z'] = gcode_['Z'].astype(float)
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


#------------------------Gcode position discretization array-------------------
#chops up g-code command into x-y-z position of laser depending on timestep
#and scanning speed


x = 0
y = 0
import sys

timestep = 0.005 #[s]
elsize = 0.05 #element size [mm]
#gcode_coords = []
start_line = 0
# sys
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

        #row_gcode = 0
        #if round(gcode_coords['x'][row_gcode],3) != round(gcode_coords['x'][row_gcode-1],3) or \
        #gcode_coords['y'][row_gcode] != gcode_coords['y'][row_gcode-1]:

            
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
                z.append(gcode_coords['Z'][row_gcode])
                t.append(d_i/speed + t[-1])
                vxi.append(vx)
                vyi.append(vy)
                p.append(power*1000)

        
        elif gcode_coords['C'][row_gcode] == 'G4':

            x.append(gcode_coords['X'][row_gcode])
            y.append(gcode_coords['Y'][row_gcode])
            z.append(gcode_coords['Z'][row_gcode])
            t.append(gcode_coords['F'][row_gcode] + t[-1])
            vxi.append(0)
            vyi.append(0)
            p.append(p[-1])

            
        elif gcode_coords['C'][row_gcode] == 'M3':
            x.append(x[-1])
            y.append(y[-1])
            z.append(z[-1])
            t.append(t[-1]+dt)
            vxi.append(vxi[-1])
            vyi.append(vyi[-1])
            p.append(gcode_coords['L'][row_gcode]*beamV*1000)

        elif gcode_coords['C'][row_gcode] == 'M5':
            x.append(x[-1])
            y.append(y[-1])
            z.append(z[-1])
            t.append(t[-1]+dt)
            vxi.append(vxi[-1])
            vyi.append(vyi[-1])
            p.append(0)

        laser_pos_array = pd.DataFrame(zip(t,x,vxi,y,vyi,z,p),columns = ['t','x','vx','y','vy','z','p'])
        
    return laser_pos_array

