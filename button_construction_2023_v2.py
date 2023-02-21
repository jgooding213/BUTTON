'''Create a PMTINFO RATDB table with positions of PMTs arranged in a box.
Must be run from the util folder. Will create a files in the data folder.
'''
from cmath import sqrt
import numpy as np
import math
import os
import array as arr
import re
import csv
import pandas as pd
#from itertools import izip
### Default values to change ratdb geometry files

###rPMT    = 6700.0
###rPMT    = 5065.0
###rPMT    = 4065.0
#rPMT    = 5465.0
#rPMT    = 4465.0
#rPMT    = 6300.0
rPMT     = 1537.5

###zPMT    = 6700.0
###zPMT    = 5065.0
###zPMT    = 4065.0
#zPMT    = 5465.0
#zPMT    = 4465.0
#zPMT    = 6700.0
zPMT     = 1300.0

#dFIDVol = -300#-150.0 ## Arbitrary 1m buffer
dFIDVol = -350#-150.0 ## Arbitrary 1m buffer
tFIDVol = 0.0
dPSUP   = 385    
tPSUP   = 6.
tBSHEET = 2.5
#dTANK   = 935.0
#dTANK  = 535.0
#dTANK   = 1300.0
dTANK   = 250.0
tTANK   = 10.0
oTANK   = 200.
dIBEAM  = 500.
tIBEAM  = 27.0
dAIR    = 1400.0 
dCONC   = 100.0
tCONC   = 25000.0
dROCK   = 2000.0  


## Values to change for PMT arrangement. (PMTINFO)
#photocoverage = 0.10
#photocoverage = 0.1505
#photocoverage = 0.205
photocoverage = 0.15
pmtRad        = 126.5

#rPMT+tBSHEET
def geoFile(rPMT = 3498.125,  zPMT = 3498.125,\
dFIDVol = -1000.0,\
tFIDVol = 0.0,\
dPSUP   = 100.0,\
tPSUP   = 10.0,\
tBSHEET = 5.0,\
#tBSHEET = 25.0,\
dTANK   = 501.875,    
tTANK   = 500.0,\
oTANK   = 200.0,\
dAIR    = 1000.0 ,\
dIBEAM  = 500.,\
tIBEAM  = 27.0,\
dCONC   = 500.0,\
tCONC   = 25000.0,\
dROCK   = 2000.0,\
pmtCnt  =  4000.0   ):                            
## dtank With respect to PMTs

    return  f"""{{
name: "GEO",
index: "world",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "", // world volume has no mother
type: "tube",
r_max: {rPMT+dTANK+dAIR+dROCK+dAIR+rPMT}, // changed to accommodate 0.5m-thick layer of concrete on walls (L. Kneale)
size_z: {zPMT+dTANK+dAIR+dROCK+dAIR+zPMT},
position: [0.0, 0.0, 0.0],
material: "air", //rock?
invisible: 1
}}

///////////////////// Define the rock volumes. Thin slab of rock is assumed ////////////////////////

//Create a 1-m rock layer around a cylindrical cavern
{{
name: "GEO",
index: "rock_1",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "world", // world volume has no mother
type: "tube",
r_max: {rPMT+dTANK+dAIR+dROCK}, // changed to accommodate 0.5m-thick layer of concrete on walls (L. Kneale)
size_z: {zPMT+dTANK+dAIR+dROCK},
position: [0.0, 0.0, 0.0], //this will allow for the concrete layer on the floor and not on the ceiling
material: "rock",
invisible: 1
//color: [1.0,0.6,0.0,1.0],
//drawstyle: "solid"
}}


//Create a 0.5m concrete layer on the walls and base
{{
name: "GEO",
index: "rock_2",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "rock_1",
type: "tube",
r_max: {rPMT+dTANK+dAIR+dCONC}, // changed to accommodate 0.5m-thick layer of concrete on walls (L. Kneale)
size_z: {zPMT+dTANK+dAIR+dCONC},
position: [0.0, 0.0, 0.0], // this will give a concrete layer on the floor and not on the ceiling
material: "rock", // changed from "gunite" (L. Kneale)
invisible: 1
//color: [0.8,0.8,0.8,0.8],
//drawstyle: "solid"
}}
//Create the cavern space between the tank and concrete
{{
name: "GEO",
index: "cavern_1",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "rock_2",
type: "tube",
r_max: {rPMT+dTANK+dAIR}, // changed to accommodate 0.5m-thick layer of concrete on walls (L. Kneale)
size_z: {zPMT+dTANK+dAIR},
position: [0.0, 0.0, 0.0],
material: "air",
invisible: 1
}}
{{name:"GEO",
index: "ibeam",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "cavern_1",
type: "tube",
r_max: {rPMT+dTANK+dAIR-dIBEAM+tIBEAM},
size_z: {zPMT+dTANK+dAIR-dIBEAM+tIBEAM},
position: [0.0, 0.0,0.0],
material: "stainless_steel",
color: [0.96,0.95,0.27,1.0],
drawstyle: "solid",
invisible: 1
}}
{{
name: "GEO",
index: "cavern_2",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "ibeam",
type: "tube",
r_max: {rPMT+dTANK+dAIR-dIBEAM},
size_z: {zPMT+dTANK+dAIR-dIBEAM}, 
position: [0.0, 0.0, 0.0],
material: "air",
color: [0.85, 0.72, 1.0, 0.5],
invisible: 1
}}
////////////////////////////////// Define the rock volumes done.///////////////////////////////////
{{
name: "GEO",
index: "tank",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "cavern_2",
type: "tube",
r_max: {rPMT+dTANK+tTANK}, // changed to accommodate 0.5m-thick layer of concrete on walls (L. Kneale)
size_z: {zPMT+dTANK+tTANK},
position: [0.0, 0.0, 0.0],
material: "stainless_steel",
color: [0.43,0.70,0.90,1.0],
drawstyle: "solid"
invisible: 1 // omitted for visualization
}}
{{
name: "GEO",
index: "detector_veto1",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "tank",
type: "tube",
r_max: {rPMT+dTANK}, // changed to accommodate 0.5m-thick layer of concrete on walls (L. Kneale)
size_z: {zPMT+dTANK},
position: [0.0, 0.0, 0.0],
material: "doped_water",
color: [0.2,0.2,0.9,0.2],
drawstyle: "solid"
invisible: 1 // omitted for visualization
}}
{{
name: "GEO",
index: "black_sheet",
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "detector_veto1",
type: "tube",
r_max: {rPMT+(tBSHEET*50)}, // changed to accommodate 0.5m-thick layer of concrete on walls (L. Kneale)
size_z: {zPMT}, // changed to put tarp behind washer plates // r + z both changed to put tarp completely around encapsulations (Temp?) JG
position: [0.0, 0.0, 0.0],
material: "polypropylene",
//color: [1.,1.,0.,1.0],
color: [0.2,0.2,0.9,0.2],
drawstyle: "solid",
//invisible: 1 // omitted for visualization
}}
{{
name: "GEO",
index: "detector_target_gb",// gb: gamma buffer
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "black_sheet",
type: "tube",
r_max: {rPMT+(tBSHEET*49)}, // changed to put tarp behind washer plates
size_z: {zPMT}, // r + z both changed to put tarp completely around encapsulations (Temp?) JG
position: [0.0, 0.0, 0.0],
material: "doped_water",
color: [0.2,0.2,0.9,0.2],
drawstyle: "solid"
invisible: 1 // omitted for visualization
}}
{{
name: "GEO",
index: "detector_target_fv", // fv : fiducial volume
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "detector_target_gb", // gb : gamma buffer
type: "tube",
r_max: {rPMT+dFIDVol}, // changed to accommodate 0.5m-thick layer of concrete on walls (L. Kneale)
size_z: {zPMT+dFIDVol},
position: [0.0, 0.0, 0.0],
material: "doped_water",
color: [1,0.2,0.9,0.2],
drawstyle: "solid"
invisible: 1 // omitted for visualization
}}
{{
//Bergevin: Set the interface were reflection can occur. Must make sure volume1 and volume2
//are in the correct order
name: "GEO",
index: "midsurface_black_sheet",
valid_begin: [0, 0],
valid_end: [0, 0],
invisible: 1, // omitted for visualization
mother: "black_sheet", //not used but needs to be a valid name, parent of 'a' and 'b' would be best choice
type: "border",
volume1: "detector_target_gb",
volume2: "black_sheet",
reverse: 1, //0 only considers photons from a->b, 1 does both directions
surface: "nonreflective_tarp",
}}
{{
//Bergevin: Set the interface were reflection can occur. Must make sure volume1 and volume2
//are in the correct order
name: "GEO",
index: "midsurface_tank",
valid_begin: [0, 0],
valid_end: [0, 0],
invisible: 1, // omitted for visualization
mother: "tank", //not used but needs to be a valid name, parent of 'a' and 'b' would be best choice
type: "border",
volume1: "detector_veto1",
volume2: "tank",
reverse: 1, //0 only considers photons from a->b, 1 does both directions
surface: "reflective_tarp",
}}
{{
name: "GEO",
index: "inner_pmts",
enable: 1,
valid_begin: [0, 0],
valid_end: [0, 0],
mother: "detector_target_gb",
type: "pmtarray",
end_idx: {int(pmtCnt-1)}, //idx of the last pmt
//end_idx: 0, //idx of the last pmt
start_idx: 0, //idx of the first pmt
pmt_model: "r7081pe",
mu_metal: 0,
mu_metal_material: "aluminum",
mu_metal_surface: "aluminum",
light_cone: 0,
light_cone_material: "aluminum",
light_cone_surface: "aluminum",
light_cone_length: 17.5,
light_cone_innerradius: 12.65,
light_cone_outerradius: 21,
light_cone_thickness: 0.2,
black_sheet_offset: 300.0, //30 cm default black tarp offset
black_sheet_thickness: 10.0, //1 cm default black tarp thickness
pmt_detector_type: "idpmt",
sensitive_detector: "/mydet/pmt/inner",
efficiency_correction: 0.90000,
pos_table: "PMTINFO", //generated by positions.nb
orientation: "manual",
orient_point: [0.,0.,0.],
color: [0.3,0.5, 0.0, 0.2],
//invisible: 1 // omitted for visualization,
encapsulation: 1
}}

{{
name: "GEO",
index: "uprightAssemblies",
enable: 1,
valid_begin: [0, 0],
valid_end: [0, 0],
//mother: "detector_veto1",
mother: "detector_target_gb",
//mother: "tank",
type: "tubearray",
end_idx: 8, //idx of the last pmt
start_idx: 0, //idx of the first pmt
r_max:  50.00,
r_min:  40.40,
size_z: 1250.0,
pos_table: "FRAMEUPRIGHTSINFO", //generated by button-construction (past 19/Aug/22)
orientation: "manual",
orient_point: [0.,0.,0.],
material: "stainless_steel",
drawstyle: "solid",
//color: [1.0,0.0,0.0,0.0],
//invisible: 1 // omitted for visualization
}}


"""





def cylinder(rPMT,zPMT,_type=1,inv = 1.0, inner = 1.0,spacing = 540.,\
delta = 250.,tolerance = 200., photocoverage=0.205, pmtRad = 126.5):
    x,y,z,dx,dy,dz,type = [],[],[],[],[],[],[]
    cnt = 0 
    #setup variables
    pmtArea = 3.14159265359*pmtRad*pmtRad
    length = 540
    print('Length/spacing: ',length,spacing) 
    heights = [-1201.+18., -785., -295., 295., 785., 1201.-18.]
    spacing = 540 
    rangeX = rangeY = int((2.0*rPMT)/spacing)
    height =  zPMT-50
    radius  = rPMT-50
    radius_end = rPMT

    innerSquare = 300.0
    innerRect = 1025.0
    outerSquare = 662.5

    _cntB = 0
    #top (inner square)
    _cntB = 0
    _x = innerSquare
    _y = innerSquare
    x.append(_x)
    y.append(_y)
    x.append(-_x)
    y.append(_y)
    x.append(_x)
    y.append(-_y)
    x.append(-_x)
    y.append(-_y)
    _cntB+=4
    cnt+=4

    #top (horizontal rectangle)
    _x = innerRect
    _y = innerSquare
    x.append(_x)
    y.append(_y)
    x.append(-_x)
    y.append(_y)
    x.append(_x)
    y.append(-_y)
    x.append(-_x)
    y.append(-_y)
    _cntB+=4
    cnt+=4
    
    #top (vertical rectangle)
    _x = innerSquare
    _y = innerRect
    x.append(_x)
    y.append(_y)
    x.append(-_x)
    y.append(_y)
    x.append(_x)
    y.append(-_y)
    x.append(-_x)
    y.append(-_y)
    _cntB+=4
    cnt+=4

    #top (outer square)
    _x = outerSquare
    _y = outerSquare
    x.append(_x)
    y.append(_y)
    x.append(-_x)
    y.append(_y)
    x.append(_x)
    y.append(-_y)
    x.append(-_x)
    y.append(-_y)
    _cntB+=4
    cnt+=4
    
    height = heights[5] #JG
    print("HEIGHT " ,height)
    #z,dx,dy,dz,type = [],[],[],[],[]
    print(_cntB)
    for i in range(_cntB):
        z.append(height)
        dx.append(0.0)
        dy.append(0.0)
        dz.append(-1.0*inv)
        type.append(_type)

    _cntB = 0
    ##bottom (inner square)
    _x = innerSquare
    _y = innerSquare
    x.append(_x)
    y.append(_y)
    x.append(-_x)
    y.append(_y)
    x.append(_x)
    y.append(-_y)
    x.append(-_x)
    y.append(-_y)
    _cntB+=4
    cnt+=4

    #bottom (horizontal rectangle)
    _x = innerRect
    _y = innerSquare
    x.append(_x)
    y.append(_y)
    x.append(-_x)
    y.append(_y)
    x.append(_x)
    y.append(-_y)
    x.append(-_x)
    y.append(-_y)
    _cntB+=4
    cnt+=4
    
    #bottom (vertical rectangle)
    _x = innerSquare
    _y = innerRect
    x.append(_x)
    y.append(_y)
    x.append(-_x)
    y.append(_y)
    x.append(_x)
    y.append(-_y)
    x.append(-_x)
    y.append(-_y)
    _cntB+=4
    cnt+=4

    ##bottom (outer square)
    _x = outerSquare
    _y = outerSquare
    x.append(_x)
    y.append(_y)
    x.append(-_x)
    y.append(_y)
    x.append(_x)
    y.append(-_y)
    x.append(-_x)
    y.append(-_y)
    _cntB+=4
    cnt+=4

    height = heights[0] #JG
    print("HEIGHT " ,height)
    #z,dx,dy,dz,type = [],[],[],[],[]
    print(_cntB)
    for i in range(_cntB):
        z.append(height)
        dx.append(0.0)
        dy.append(0.0)
        dz.append(1.0*inv)
        type.append(_type)

    #side
    nRings  = 2
    colums = 8
    
    _dTheta = 2.0 * math.pi / colums
    print('Columns: ',colums,'\nNrings: ',nRings,'\ndelta Theta: ',_dTheta)
 
    radius-=88.
    for _i in range(nRings):
        #        print((_i*500.)+250.,-((_i*500.)+250.))
        for _j in range(colums):
            _theta,_z = _j*_dTheta,(_i*spacing)+delta
            _z = heights[_i+1] #JG
            print("from the bottom: ",_z)
            print("from the top: ",-_z)            
            x.append(radius * math.cos(_theta) + math.sin(_theta) * spacing/2)
            y.append(radius * math.sin(_theta) - math.cos(_theta) * spacing/2)
            z.append(_z)
            dx.append(-math.cos(_theta)*inv)
            dy.append(-math.sin(_theta)*inv)
            dz.append(0.0)
            x.append(radius * math.cos(_theta) - math.sin(_theta) * spacing/2)
            y.append(radius * math.sin(_theta) + math.cos(_theta) * spacing/2)
            z.append(_z)
            dx.append(-math.cos(_theta)*inv)
            dy.append(-math.sin(_theta)*inv)
            dz.append(0.0)
            type.append(_type)
            type.append(_type)
            #_z = heights[4-_i] #JG
            #print("from the bottom: ",heights[4-_i])
            x.append(radius * math.cos(_theta) + math.sin(_theta) * spacing/2)
            y.append(radius * math.sin(_theta) - math.cos(_theta) * spacing/2)
            z.append(-_z)
            dx.append(-math.cos(_theta)*inv)
            dy.append(-math.sin(_theta)*inv)
            dz.append(0.0)
            x.append(radius * math.cos(_theta) - math.sin(_theta) * spacing/2)
            y.append(radius * math.sin(_theta) + math.cos(_theta) * spacing/2)
            z.append(-_z)
            dx.append(-math.cos(_theta)*inv)
            dy.append(-math.sin(_theta)*inv)
            dz.append(0.0)
            type.append(_type)
            type.append(_type)
            cnt+=4
            
    # print(len(dy))
#    print(cnt)
    return x,y,z,dx,dy,dz,type,cnt

def pmtinfoFile(_x=rPMT, _z=zPMT,\
photocoverage=photocoverage, pmtRad = pmtRad):

    x,y,z,dx,dy,dz,type,cnt = cylinder(rPMT,zPMT,photocoverage=photocoverage, pmtRad = pmtRad)
    
    pmt_info = "{\n"
    pmt_info += f"//// Total number of inner PMTs : {cnt}\n"
    pmt_info += f"//// Total number of veto PMTs : 0\n"
    pmt_info += f"\"name\": \"PMTINFO\",\n"
    pmt_info += f"\"valid_begin\": [0, 0],\n"
    pmt_info += f"\"valid_end\": [0, 0],\n"
    pmt_info += f"\"x\":     {x},\n"
    pmt_info += f"\"y\":     {y},\n"
    pmt_info += f"\"z\":     {z},\n"
    pmt_info += f"\"dir_x\": {dx},\n"
    pmt_info += f"\"dir_y\": {dy},\n"
    pmt_info += f"\"dir_z\": {dz},\n"
    pmt_info += f"\"type\": {type},\n"
    pmt_info += "}"




    xfloat = list(np.float_(x))

    #rowString = 'this is a string  ' + str(x[0]) + ', ' + str(y[0]) + ', ' + str(z[0]) + ', ' + str(dx[0]) + ', ' + str(dy[0]) + ', ' + str(dz[0])



    #print(rowString)
    topNBottomPMTS = 32

    encapsulationCSV = open('data/button_frame/encapsulationCSV.csv','w')
    #writer = csv.writer(encapsulationCSV)
    #encapsulationCSV.write("x[0:31],y[0:31],z[0:31],dx[0:31],dy[0:31],dz[0:31], \n")
    encapsulationCSV.write("x,y,z,dx,dy,dz, \n")   
    #for i in range(len(x)):
    for i in range(topNBottomPMTS):
        rowString =  str(x[i]) + ', ' + str(y[i]) + ', ' + str(z[i]) + ', ' + str(dx[i]) + ', ' + str(dy[i]) + ', ' + str(dz[i]) +'\n'
        #print(rowString)
        encapsulationCSV.write(rowString)

    encapsulationCSV.close()

# lets save this in a nice readable format for later...

    return pmt_info,cnt


######################################################## FRAME INFO






def cylinder2(rPMT,zPMT,_type=1,inv = 1.0, inner = 1.0,spacing = 540.,\
delta = 250.,tolerance = 200., photocoverage=0.205, pmtRad = 126.5, frameRad = 1660):
    #x,y,z,dx,dy,dz,type = [],[],[],[],[],[],[]
    # tubeLocx,tubeLocy,tubeLocz  = [],[],[]
    # tubeDirx,tubeDiry,tubeDirz  = [],[],[]

    # moduleLocx,moduleLocy,moduleLocz = [],[],[]
    # moduleDirx,moduleDiry,moduleDirz = [],[],[]

    Locx,Locy,Locz  = [],[],[]
    Dirx,Diry,Dirz  = [],[],[]
    type = []

    cnt = 0 
    #setup variables
    pmtArea = 3.14159265359*pmtRad*pmtRad
    length = 540
    print('Length/spacing: ',length,spacing) 
    spacing = 540 
    rangeX = rangeY = int((2.0*rPMT)/spacing)
    height =  zPMT-50
    radius  = rPMT-50
    radius_end = rPMT

    innerSquare = 300.0
    innerRect = 1025.0
    outerSquare = 662.5


    #Uprights
    thetaUprights = math.pi/8
    #tubeLocx = frameRad * math.sin(thetaUprights)
    #tubeLocy = frameRad * math.cos(thetaUprights)
    #tubeLocz = 0 




    cnt = 0
    for i in range(8):
          Locx.append(frameRad * math.sin((thetaUprights*2*i - thetaUprights)))
          Locy.append(frameRad * math.cos((thetaUprights*2*i - thetaUprights)))         
          Locz.append(0.)
          #Dirx.append(math.sin((thetaUprights*2*i - thetaUprights - math.pi)))
          #Diry.append(math.cos((thetaUprights*2*i - thetaUprights - math.pi)))
          Dirx.append(0.)
          Diry.append(0.)
          Dirz.append(0.)          
          type.append(0)      
          cnt+=1



    # print(Locx)
    # print(Locy)
    # print(Locz)     
    # print(Dirx)
    # print(Diry)
    # print(Dirz) 
    # print(type)              
    return Locx,Locy,Locz,Dirx,Diry,Dirz,type,cnt
    #return tubeLocx,tubeLocy,tubeLocz,tubeDirx,tubeDiry,tubeDirz,type,cnt




def frameUinfoFile(_x=rPMT, _z=zPMT,\
photocoverage=photocoverage, pmtRad = pmtRad):

    Locx,Locy,Locz,Dirx,Diry,Dirz,type,cnt = cylinder2(rPMT,zPMT,photocoverage=photocoverage, pmtRad = pmtRad)

    frameU_info = "{\n" 
    frameU_info += f"//// Total number of frame elements : {cnt}\n"
    frameU_info += f"//// Total number of veto PMTs : 0\n"
    frameU_info += f"\"name\": \"FRAMEUPRIGHTSINFO\",\n"
    frameU_info += f"\"valid_begin\": [0, 0],\n"
    frameU_info += f"\"valid_end\": [0, 0],\n"
    frameU_info += f"\"x\":     {Locx},\n"
    frameU_info += f"\"y\":     {Locy},\n"
    frameU_info += f"\"z\":     {Locz},\n"
    frameU_info += f"\"dir_x\": {Dirx},\n"
    frameU_info += f"\"dir_y\": {Diry},\n"
    frameU_info += f"\"dir_z\": {Dirz},\n"
    frameU_info += f"\"type\": {type},\n"
    frameU_info += "}"

    return frameU_info,cnt






######### MODULES

def cylinder3(rPMT,zPMT,_type=1,inv = 1.0, inner = 1.0,spacing = 540.,\
delta = 250.,tolerance = 200., photocoverage=0.205, pmtRad = 126.5, frameRad = 1660):


    LocMx,LocMy,LocMz  = [],[],[]
    DirMx,DirMy,DirMz  = [],[],[]
    typeM = []

    cnt = 0 

    #Uprights
    thetaUprights = math.pi/8


    cnt = 0



    moduleHeight = arr.array('i', [625,-625])


    nRings  = 2
    for j in range(2):
        for i in range(8):
            LocMx.append(frameRad * math.sin((i)*2*thetaUprights))
            LocMy.append(frameRad * math.cos((i)*2*thetaUprights))
            LocMz.append(moduleHeight[j])             
            DirMx.append(math.sin((thetaUprights*2*i - math.pi)))
            DirMy.append(math.cos((thetaUprights*2*i - math.pi)))         
            typeM.append(1)  
            cnt+=1


    # print(LocMx)
    # print(LocMy)
    # print(LocMz)     
    # print(DirMx)
    # print(DirMy)
    # print(DirMz) 
    # print(typeM)              
    return LocMx,LocMy,LocMz,DirMx,DirMy,DirMz,typeM,cnt
    #return tubeLocx,tubeLocy,tubeLocz,tubeDirx,tubeDiry,tubeDirz,type,cnt

def frameMinfoFile(_x=rPMT, _z=zPMT,\
photocoverage=photocoverage, pmtRad = pmtRad):

    LocMx,LocMy,LocMz,DirMx,DirMy,DirMz,typeM,cnt = cylinder3(rPMT,zPMT,photocoverage=photocoverage, pmtRad = pmtRad)

    frameM_info = "{\n" 
    frameM_info += f"//// Total number of frame elements : {cnt}\n"
    frameM_info += f"//// Total number of veto PMTs : 0\n"
    frameM_info += f"\"name\": \"FRAMEMODULESINFO\",\n"
    frameM_info += f"\"valid_begin\": [0, 0],\n"
    frameM_info += f"\"valid_end\": [0, 0],\n"
    frameM_info += f"\"x\":     {LocMx},\n"
    frameM_info += f"\"y\":     {LocMy},\n"
    frameM_info += f"\"z\":     {LocMz},\n"
    frameM_info += f"\"dir_x\": {DirMx},\n"
    frameM_info += f"\"dir_y\": {DirMy},\n"
    frameM_info += f"\"dir_z\": {DirMz},\n"
    frameM_info += f"\"type\": {typeM},\n"
    frameM_info += "}"

    return frameM_info,cnt






###########################################################




#####  main code ###########
_frameUinfo,cnt = frameUinfoFile(_x=rPMT, _z=zPMT,\
photocoverage=photocoverage, pmtRad = pmtRad)

_frameMinfo,cnt = frameMinfoFile(_x=rPMT, _z=zPMT,\
photocoverage=photocoverage, pmtRad = pmtRad)

_pmtinfo,cnt = pmtinfoFile(_x=rPMT, _z=zPMT,\
photocoverage=photocoverage, pmtRad = pmtRad)

_geoFile = geoFile(rPMT = rPMT, zPMT = zPMT,dFIDVol = dFIDVol,tFIDVol = tFIDVol,dPSUP   = dPSUP,tPSUP = tPSUP,\
tBSHEET = tBSHEET, dTANK = dTANK,tTANK = tTANK,dAIR = dAIR ,dCONC = dCONC,tCONC = tCONC, dROCK = dROCK, pmtCnt = cnt )





# Module looks like so cross section is 40x40 mm
#
#     |------------------------horL-------------------------|
#     =======================================================
#     #                         |---------horSpace---------| --
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  verL
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  |
#     #                         #                          #  --
#     =======================================================








topnBase = arr.array('i', [1250,-1250])







positions = 4


sides = 8
horL = 1140.
verL = 980.
# xposU = 0
# yposU = 0 
# lenU = 0
crossSec = 40.
nRings = 2


frameRad = 1487.5 #1554. #1660.
#Uprights
thetaUprights = math.pi/sides # 45 degrees
counter = 0
counterSide = 0
outerFrameRad = 1437.1


encapArrX = arr.array('d',[275.,275.,-275.,-275.])
encapArrZ = arr.array('d',[240.,-240.,-240.,240.])

so = 550 * math.sin(math.pi/4)     #389 # Side offset
print("SO ")
print(so)
p45 = frameRad * math.sin(math.pi/4) # 1051.8 # position xy at 45 degree


lposX = arr.array('d', [550, p45-so , p45+so , 1487.5, 1487.5, p45+so, p45-so, 550, -550, -p45+so, -p45-so,  -1487.5, -1487.5,  -p45-so, -p45+so, -550, 550, p45-so , p45+so , 1487.5, 1487.5, p45+so, p45-so, 550, -550, -p45+so, -p45-so,  -1487.5, -1487.5,  -p45-so, -p45+so, -550 ]) # 1051.8 is basic 45 degree
lposY = arr.array('d', [1487.5, p45+so , p45-so ,  550, -550, -p45+so, -p45-so, -1487.5, -1487.5, -p45-so,-p45+so,  -550, + 550, p45-so,  p45+so, 1487.5, 1487.5, p45+so , p45-so ,  550, -550, -p45+so, -p45-so, -1487.5, -1487.5, -p45-so,-p45+so,  -550, + 550, p45-so,  p45+so, 1487.5])
dirX = arr.array('d', [0. ,1. ,1., 0.,0. ,1. ,1., 0.,0. ,1. ,1., 0.,0. ,1. ,1., 0., 0. ,1. ,1., 0.,0. ,1. ,1., 0.,0. ,1. ,1., 0.,0.,1. ,1. ,0. ])
dirY = arr.array('d', [0. ,1. ,1., 0.,0. ,1. ,1., 0.,0. ,1. ,1., 0.,0. ,1. ,1., 0., 0. ,1. ,1., 0.,0. ,1. ,1. ,0.,0. ,1. ,1., 0.,0.,1. ,1. ,0. ])


dir = arr.array('d', [0. , 45.,45. ,90.,90., 135., 135., 180., 180. , 225., 225. , 270., 270. , 315., 315., 0., 0. , 45.,45. ,90.,90., 135., 135., 180., 180. , 225., 225. , 270., 270. , 315., 315., 0.  ])


counter2 = 0



moduleHeight = arr.array('i', [540,-540])
for j in range(nRings):    
    for i in range(sides):
        Uz = (moduleHeight[j])    
        #DirUx = math.sin((thetaUprights*2*i - math.pi/4))
        #DirUy = math.cos((thetaUprights*2*i - math.pi/4))   
        #DirUz = 0# i*math.pi
        DirUx = 70.
        DirUy = 0.
        DirUz = 0.#i*45#(180.0/math.pi)*math.pi/4#(math.sin((math.pi/8)))
        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"sideUpright{counter2}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"
        _geoFile +=  f"size: [{crossSec/2},{crossSec/2} , {verL/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{lposX[counter2]}, {lposY[counter2]}, {Uz}],\n"  
        #_geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz}],\n"   
        #_geoFile +=  f"rotation: [{dirX[counter2]}, {dirY[counter2]}, 0.],\n"          
        _geoFile +=  f"rotation: [0., 0., {dir[counter2]}],\n"         
        #_geoFile +=  "position: [0.0, 0.0, 0.0],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        #_geoFile +=  "drawstyle: \"wireframe\",\n"
        _geoFile +=  "drawstyle: \"solid\",\n"        
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"


        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"sideUpright{counter2+1}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"
        _geoFile +=  f"size: [{crossSec/2},{crossSec/2} , {verL/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{lposX[counter2+1]}, {lposY[counter2+1]}, {Uz}],\n"  
        #_geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz}],\n"  
        #_geoFile +=  "rotation: [0.,0.,0.],\n"          
        #_geoFile +=  f"rotation: [{dirX[counter2]}, {dirY[counter2]}, 0.],\n"    
        _geoFile +=  f"rotation: [0., 0., {dir[counter2+1]}],\n"     
        #_geoFile +=  "position: [0.0, 0.0, 0.0],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        #_geoFile +=  "drawstyle: \"wireframe\",\n"
        _geoFile +=  "drawstyle: \"solid\",\n"        
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"        
        counter2+=2





# #HERE
for j in range(nRings):    
    for i in range(sides):

        Ux = (frameRad * math.sin((i)*2*thetaUprights))
        Uy = (frameRad * math.cos((i)*2*thetaUprights))
        OUx = (outerFrameRad * math.sin((i)*2*thetaUprights))
        OUy = (outerFrameRad * math.cos((i)*2*thetaUprights))        
        Uz = (moduleHeight[j])             
        DirUx = 0.0#180*(math.sin((thetaUprights*2*i - math.pi)))/math.pi
        DirUy = 0.0#180*(math.cos((thetaUprights*2*i - math.pi)))/math.pi
        DirUz = i*(180.0/math.pi)*math.pi/4#(math.sin((math.pi/8)))
        print(DirUz)
        LRXFRAME = 550 * math.sin((math.pi/4)-((i)*2*thetaUprights))
        LRYFRAME = 550 * math.cos((math.pi/4)-((i)*2*thetaUprights))

        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"frameBar{counter}_{j}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"      
        #_geoFile +=  f"size: [{horL/2}, {crossSec/2}, {((2*crossSec)+(verL))/2}], // mm, half-length\n"
        _geoFile +=  f"size: [{1167./2}, {crossSec/2},  {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{OUx}, {OUy}, {topnBase[j]}],\n"  
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz}],\n"
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"wireframe\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"

        
        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"centerUpright{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"
        _geoFile +=  f"size: [{crossSec/2},{crossSec/2} , {verL/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"  
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz}],\n"      
        #_geoFile +=  "position: [0.0, 0.0, 0.0],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"wireframe\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"



        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"topHorizontal{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"
        _geoFile +=  f"size: [{horL/2}, {crossSec/2},{crossSec/2}], // mm, half-length\n"
        # _geoFile +=  f"position: [0.0, 0.0, {(verL /2) + 20.}],\n"      
        _geoFile +=  f"position: [{Ux}, {Uy}, {Uz+(verL /2) + 20.}],\n" 
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz}],\n"                      
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"wireframe\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"

        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"bottomHorizontal{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"
        _geoFile +=  f"size: [{horL/2}, {crossSec/2},{crossSec/2}], // mm, half-length\n"
        # _geoFile +=  f"position: [0.0, 0.0, {(verL /2) + 20.}],\n"      
        _geoFile +=  f"position: [{Ux}, {Uy}, {Uz-(verL /2) - 20.}],\n" 
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz}],\n"    
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"wireframe\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"






        # _geoFile += "{\n"
        # _geoFile += "name: \"GEO\",\n"
        # _geoFile +=  f"index: \"rightUpright{counter}\",\n"       
        # _geoFile += "valid_begin: [0, 0],\n"
        # _geoFile += "valid_end: [0, 0],\n"                
        # _geoFile +=  f"mother: \"frameModule{counter}\",\n"
        # _geoFile +=  "type: \"box\",\n"
        # _geoFile +=  f"size: [{crossSec/2},{crossSec/2} , {verL/2}], // mm, half-length\n"
        # _geoFile +=  f"position: [{(horL - crossSec) /2}, 0.0, 0.0],\n"        
        # _geoFile +=  "material: \"stainless_steel\"\n"
        # _geoFile +=  "drawstyle: \"wireframe\",\n"
        # #_geoFile +=  "invisible: 1,\n" 
        # _geoFile += "}\n"




        # _geoFile += "{\n"
        # _geoFile += "name: \"GEO\",\n"
        # _geoFile +=  f"index: \"topHorizontal{counter}\",\n"       
        # _geoFile += "valid_begin: [0, 0],\n"
        # _geoFile += "valid_end: [0, 0],\n"                
        # _geoFile +=  f"mother: \"frameModule{counter}\",\n"
        # _geoFile +=  "type: \"box\",\n"
        # _geoFile +=  f"size: [{horL/2}, {crossSec/2},{crossSec/2}], // mm, half-length\n"
        # _geoFile +=  f"position: [0.0, 0.0, {(verL /2) + 20.}],\n"        
        # _geoFile +=  "material: \"stainless_steel\"\n"
        # _geoFile +=  "drawstyle: \"wireframe\",\n"
        # #_geoFile +=  "invisible: 1,\n" 
        # _geoFile += "}\n"


        # _geoFile += "{\n"
        # _geoFile += "name: \"GEO\",\n"
        # _geoFile +=  f"index: \"bottomHorizontal{counter}\",\n"       
        # _geoFile += "valid_begin: [0, 0],\n"
        # _geoFile += "valid_end: [0, 0],\n"                
        # _geoFile +=  f"mother: \"frameModule{counter}\",\n"
        # _geoFile +=  "type: \"box\",\n"
        # _geoFile +=  f"size: [{horL/2}, {crossSec/2},{crossSec/2}], // mm, half-length\n"
        # _geoFile +=  f"position: [0.0, 0.0, {-(verL + crossSec) /2}],\n"        
        # _geoFile +=  "material: \"stainless_steel\"\n"
        # _geoFile +=  "drawstyle: \"wireframe\",\n"
        # #_geoFile +=  "invisible: 1,\n"       
        # _geoFile += "}\n"

        counter+=1





centerSqu = 250
innerSqu = 1400
outerSqu = 2426
connectSqu = 470
connectSpacing = 420.

print("5 + is WEIRD LITTLE FIX OF OVERLAP, DOUBLE CHECK")
#
conOffset = 5+(connectSpacing/2)#*math.cos(math.pi/4)
#conOffset = (connectSpacing/2)#*math.cos(math.pi/4)



sides = 4

counter = 0


for j in range(nRings):    
    for i in range(sides):

        CUx = ((centerSqu-crossSec) * math.sin((i)*0.5*math.pi))/2
        CUy = ((centerSqu-crossSec)  * math.cos((i)*0.5*math.pi))/2

        #HI = (innerSqu-crossSec)/2
        IUx = (innerSqu-crossSec)/2 * math.cos((math.pi*i/2) + math.pi/4) # 
        IUy = (innerSqu-crossSec)/2 * math.sin((math.pi*i/2) + math.pi/4) #  
        OUx = (outerSqu-crossSec)/2 * math.cos((math.pi*i/2) + math.pi/4) # 
        OUy = (outerSqu-crossSec)/2 * math.sin((math.pi*i/2) + math.pi/4) #        

        CoUx = (IUx + OUx )/2
        CoUy = (IUy + OUy )/2

####### bottom and top base

        Uz = (topnBase[j])    

        DirUx = 0.0
        DirUy = 0.0
        DirUz = i*90
        print(i)
        print("//// xyz : ",IUx, IUy, Uz)
        print("//// Dir : ", DirUz+45)    


# #back to structure

        print("THIS IS THE INNER SQU NUMBER ", counter )
        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"innerSquSide{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{(1179.)/2}, {crossSec/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{IUx}, {IUy}, {Uz}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz+45}],\n"
#        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"solid\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"





        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"outerSquSide{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"   
        #_geoFile +=  f"size: [{(outerSqu -(2*crossSec))/2}, {crossSec/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"size: [{1561./2}, {crossSec/2}, {crossSec/2}], // mm, half-length\n"        
        _geoFile +=  f"position: [{OUx}, {OUy}, {Uz}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz+45}],\n"
#        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"solid\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"




# check the location dimension of this
        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"connectSqu{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{connectSqu/2}, {(connectSpacing+ crossSec)/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{CoUx}, {CoUy}, {Uz}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz-45}],\n"
       # _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"wireframe\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"



# check the location dimension of this
        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"connectSquL{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  f"mother: \"connectSqu{counter}\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{connectSqu/2}, {crossSec/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [0., {conOffset}, {0}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
#        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz-45}],\n"
        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"solid\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"


# check the location dimension of this
        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"connectSquR{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  f"mother: \"connectSqu{counter}\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{connectSqu/2}, {crossSec/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [0., {-conOffset},  {0}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
#        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz-45}],\n"
        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"solid\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"


        counter+=1






# now to do differences between top and bottom structures
    sqSides = 4


    topBarLen = 1284.
    topBarMid =  (250 + topBarLen)/2

    bottomBarLen = 1360.
    bottomBarMid = 50 + (bottomBarLen/2)
    bottom_width = 100.

    counter=0
    for i in range(sqSides):
        Uz = (topnBase[0])   
        CUx = ((centerSqu-crossSec) * math.sin((i)*0.5*math.pi))/2
        CUy = ((centerSqu-crossSec)  * math.cos((i)*0.5*math.pi))/2

        TUx = ((topBarMid) * math.sin((i)*math.pi/2))
        TUy = ((topBarMid)  * math.cos((i)*math.pi/2))       
        BUx = ((bottomBarMid) * math.sin((i)*math.pi/2))
        BUy = ((bottomBarMid)  * math.cos((i)*math.pi/2))   

        DirUx = 0.0
        DirUy = 0.0
        DirUz = i*90
        print("this is the number", i , "  pos of bar " , TUx, " , " , TUy)


        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"centerSquSide{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{(centerSqu -(2*crossSec))/2}, {crossSec/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{CUx}, {CUy}, {Uz}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz}],\n"
#        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"solid\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"

        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"topCross{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{(topBarLen)/2}, {crossSec/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{TUx}, {TUy}, {topnBase[0]}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz+90}],\n"
#        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"solid\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"

        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"bottomCross{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  "mother: \"tank\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{(bottomBarLen)/2}, {bottom_width/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{BUx}, {BUy}, {topnBase[1]}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
        _geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz+90}],\n"
#        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"doped_water\"\n"
        _geoFile +=  "drawstyle: \"wireframe\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"

        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"bottomCrossL{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  f"mother: \"bottomCross{counter}\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{bottomBarLen/2}, {(crossSec)/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{0.}, {30.}, {0.}],\n"        
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
        #_geoFile +=  f"rotation: [{DirUx}, {DirUy}, {}],\n"
#        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"solid\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"

        _geoFile += "{\n"
        _geoFile += "name: \"GEO\",\n"
        _geoFile +=  f"index: \"bottomCrossR{counter}\",\n"       
        _geoFile += "valid_begin: [0, 0],\n"
        _geoFile += "valid_end: [0, 0],\n"                
        _geoFile +=  f"mother: \"bottomCross{counter}\",\n"
        _geoFile +=  "type: \"box\",\n"   
        _geoFile +=  f"size: [{bottomBarLen/2}, {(crossSec)/2}, {crossSec/2}], // mm, half-length\n"
        _geoFile +=  f"position: [{0.}, {-30.}, {0.}],\n"         
        #_geoFile +=  f"position: [{Ux}, {Uy}, {Uz}],\n"     
        #_geoFile +=  f"rotation: [{DirUx}, {DirUy}, {DirUz+90}],\n"
#        _geoFile +=  f"rotation: [0., 0., 0.],\n"        
        _geoFile +=  "material: \"stainless_steel\"\n"
        _geoFile +=  "drawstyle: \"solid\",\n"
        #_geoFile +=  "invisible: 1,\n" 
        _geoFile += "}\n"

        #counter+=1 




        counter+=1 


                                                                                        
print("real values need to be given")                                                                                                                                                    
                  
    
outerDia = 200.
innerDia = 192.    
flangeDia = 243.
flangeZ = 20.                                                                                                                              
                                                                                                                                              



data = pd.read_csv("data/button_frame/encapsulationCSV.csv")


# THE BIG OL LOOP
print("LENGTH: " , len(data))


#print(_geoFile)
#print()
#print(_pmtinfo)


# try:
#     os.mkdir(f"../../data/button_{int((rPMT+dTANK)*2.0/1000)}m_{int((zPMT+dTANK)*2.0/1000)}m")
#     print('Created', f"../../data/button_{int((rPMT+dTANK)*2.0/1000)}m_{int((zPMT+dTANK)*2.0/1000)}m")
# except OSError as error:  
#     print(error)   

# geofile = open(f"../../data/button_{int((rPMT+dTANK)*2.0/1000)}m_{int((zPMT+dTANK)*2.0/1000)}m/button_{int((rPMT+dTANK)*2.0/1000)}m_{int((zPMT+dTANK)*2.0/1000)}m_{int(photocoverage*100)}pct.geo","w+")
# geofile.writelines(_geoFile)
# geofile.close

# pmtfile = open(f"../../data/button_{int((rPMT+dTANK)*2.0/1000)}m_{int((zPMT+dTANK)*2.0/1000)}m/PMTINFO.ratdb","w+")



pmtfile = open(f"data/button_frame/PMTINFO.ratdb","w+")
pmtfile.writelines(_pmtinfo)
pmtfile.close


# values = re.findall(r'\d+', _pmtinfo)
# print(values)



dataStr = _pmtinfo.partition('\n')[2]
data = dataStr.partition(',')[2]
print(data[30])
print(cnt)





geofile = open(f"data/button_frame/button_frame.geo","w+")
geofile.writelines(_geoFile)
geofile.close


frameUfile = open(f"data/button_frame/FRAMEUPRIGHTSINFO.ratdb","w+")
frameUfile.writelines(_frameUinfo)
frameUfile.close


frameMfile = open(f"data/button_frame/FRAMEMODULESSINFO.ratdb","w+")
frameMfile.writelines(_frameMinfo)
frameMfile.close



surfaceArea =  (2.*3.14159265359*rPMT)*(2*zPMT) + 2.0*3.14159265359*rPMT*rPMT
pmtArea     = float(cnt)*3.14159265359*pmtRad*pmtRad

print("//// Total number of inner PMTs : ",cnt)
print("//// Target Photocoverage (%) : ",photocoverage*100.)
print("//// Actual Photocoverage (%) : ",pmtArea/surfaceArea*100.)
print("//// Detector height (m) : ",(zPMT+dTANK)*2.0/1000.)
print("//// Detector diameter (m) : ",(rPMT+dTANK)*2.0/1000.)
print("////")
