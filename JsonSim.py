#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 13:38:43 2023

@author: ygaillard
"""
import os
import json
import re

import numpy as np

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
  

class JsonSim:
  
    siminformation={}
    fileExisted=False
    
    overAllFileName='siminformation.json'
    def __init__(self,fileName='siminformation.json',simName='',access='r'):
        self.overAllFileName=fileName
        self.access=access
        if 'w' in access:
            print('Open writable')
        else:
            print('Read only mode')
        try:
            f = open(self.overAllFileName)     
            # returns JSON object as 
            # a dictionary
            self.siminformation = json.load(f)
            # Closing file
            f.close()
            self.fileExisted=True
        except:
            print("No file found")
            self.siminformation={}
            path=os.getcwd()
            
            if simName=='':
                simName=re.sub('/.*/', '', path,flags=re.DOTALL)
            self.siminformation['name']=simName
            
    def __del__(self):
        if 'w' in self.access:
            jsonFile=json.dumps(self.siminformation,cls=NpEncoder)
            f = open(self.overAllFileName,'w')     
            f.write(jsonFile)
            f.close()
        else:
            print('Read only!')
        
            
    def write(self,dictName,value):
        self.siminformation[dictName]=value
        
    def delete(self,dictName):
        self.siminformation.pop(dictName)
        
    def exist(self,dictName):
        if dictName in self.siminformation.keys():
            return True
        else:
            return False
        
        

    

    
