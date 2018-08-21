#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 20:10:47 2018

@author: nisheng
"""
import numpy as np
import math

def readVectorFromFile(filePath):
    
    featureList=[]
    for line in open(filePath):
        s = line.split()
        
        l = []
        for tmp in s:
            l.append((float)(tmp))
        featureList.append(l)
    return featureList

def readFeatures():
    drugFeature = []
    proteinFeature = []
    
    fileDir = 'feature/'
    
    drugFeature = readVectorFromFile(fileDir+"drug_vector_d100.txt")
    
    proteinFeature = readVectorFromFile(fileDir+"protein_vector_d400.txt")
    
    return np.array(drugFeature),np.array(proteinFeature)



def readDrugProteinInteractionMat():
    fileDir = 'data/'
    drugProteinMat = np.array(readVectorFromFile(fileDir+"mat_drug_protein.txt"),dtype=int)
    return drugProteinMat
    