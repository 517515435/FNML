#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 20:41:28 2018

@author: nisheng
"""
import Input
import copy
import math
import numpy as np
import scipy
import scipy.special
import time
import random
from sklearn.model_selection import KFold  
from sklearn import metrics  
import jpype #from jpype import *  
import os.path  

def getVisitedPart(drugProteinMat):
    s1 = []
    s2 = []
    for i in range(drugProteinMat.shape[0]):
        for j in range(drugProteinMat.shape[1]):
            if drugProteinMat[i,j]==1:
                s1.append((i,j))
            else:
                s2.append((i,j))
    return s1,s2
    
def randomSelectNegSample(s,n):
    return random.sample(s,n)

def digamma(x):
        return scipy.special.psi(x)
    
def logistic(x):   
    tmp = 1.0/(1+math.exp(-x))
    return tmp
        
def f(t):
     return (logistic(t)-0.5)/(2*t)

jarpath = os.path.join(os.path.abspath('.'), '/root/')  

jpype.startJVM(jpype.getDefaultJVMPath(),"-ea", "-Djava.class.path=%s" % (jarpath + 'java.jar'))  
   
class FNML:

    sigma = 1
    eta_00 = 1
    eta_01 = 1
    eta_10 = 1
    eta_11 = 1
    
    


    mu = []

    v = []

    rho_00 = 0
    rho_01 = 0
    rho_10 = 0
    rho_11 = 0

    lambda_ = []

    xi = []
    
    def __init__(self,X,Y,drugProteinMat,R,l1,l2):
        self.X = X
        self.Y = Y
        self.N_d = X.shape[0]
        self.f_d = X.shape[1]
        
        self.N_t = Y.shape[0]
        self.f_t = Y.shape[1]
        
        self.P = drugProteinMat
        self.R = R
        self.l1 = l1
        self.l2 = l2
    
        print self.N_d,self.f_d,self.N_t,self.f_t

    
    
    def initParameter(self):
       
    
        self.lambda_ = np.random.random( (self.N_d,self.N_t) )
        self.xi = np.random.random( (self.N_d,self.N_t) )
       
       
        self.mu = np.zeros((self.f_d,self.f_t))
        self.v = np.zeros((self.f_d,self.f_t))
        
                
    def updateZ(self,l1_trainIndex,l2_trainIndex):
        startTime = time.time()
        
        allLambdaXiValue = np.ones( (self.N_d,self.N_t) )
        for i in range(self.N_d):
            for j in range(self.N_t):
                allLambdaXiValue[i,j] = f(self.xi[i,j])
        s=(1.0/(self.sigma^2)) 
      
        JDClass = jpype.JClass("com.Main")  

        jd = JDClass()  
    
        jd.calc(self.X,self.Y,self.P,self.l1,self.l2,l1_trainIndex,l2_trainIndex,self.lambda_,allLambdaXiValue,self.N_d,self.N_t,self.f_d,self.f_t,s)
        
        tmp_mu = jd.getMu()
        tmp_v = jd.getV()
        for i in range(self.f_d):
            for j in range(self.f_t):
                self.mu[i,j] = tmp_mu[i][j]
                self.v[i,j] = tmp_v[i][j]
               
        endTime = time.time()
        print "time",endTime,startTime,endTime-startTime
        
    def updateRho(self,l1_trainIndex,l2_trainIndex):
        
        self.rho_00 = self.eta_00
        self.rho_01 = self.eta_01
        self.rho_10 = self.eta_10
        self.rho_11 = self.eta_11
     
        for index in l2_trainIndex:
                i = self.l2[index][0]
                j = self.l2[index][1]
             
                self.rho_01 += (1-self.lambda_[i,j])
                self.rho_11 += self.lambda_[i,j]
        for index in l1_trainIndex:
                self.rho_10 += 1
    
    def updateLambda(self,l1_trainIndex,l2_trainIndex):
    
        counter = 0
        for index in l2_trainIndex:
                i = self.l2[index][0]
                j = self.l2[index][1]
                s1 = 0.0
                s2 = 0.0
                
                tmp1 = digamma(self.rho_11)
                tmp2 = -1*digamma(self.rho_10+self.rho_11)
                tmp3 = self.X[i].dot(self.mu).dot(self.Y[j].T)
                s1 = tmp1+tmp2+tmp3
              
                
                tmp1 = digamma(self.rho_01)
                tmp2 = -1*digamma(self.rho_00+self.rho_01)
                s2 = tmp1+tmp2
                
                if s1>10 or s2>10:
                    s1 -= max(s1,s2)
                    s2 -= max(s1,s2)
                    counter+=1
                tmp1 = math.exp(s1)
                tmp2 = math.exp(s2)
                self.lambda_[i,j] = tmp1/(tmp1+tmp2)
                
        print "counter:",counter
    
    def updateXi(self,l1_trainIndex,l2_trainIndex):
        
        for index in l1_trainIndex:
                i = self.l1[index][0]
                j = self.l1[index][1]
               
                tmp=self.X[i].dot(self.mu)
                
                self.xi[i,j] = abs(tmp.dot(self.Y[j].T))
                
        for index in l2_trainIndex:
                i = self.l2[index][0]
                j = self.l2[index][1]
               
                tmp=self.X[i].dot(self.mu)
                
                self.xi[i,j] = abs(tmp.dot(self.Y[j].T))
    def trainOnce(self,l1_trainIndex,l2_trainIndex):
        self.updateZ(l1_trainIndex,l2_trainIndex)
        self.updateRho(l1_trainIndex,l2_trainIndex)
        self.updateLambda(l1_trainIndex,l2_trainIndex)
        self.updateXi(l1_trainIndex,l2_trainIndex)
        
    def predictScore(self,l1_testIndex,l2_testIndex):
        predictScoreList = []
        trueScoreList = []
      
        for index in l1_testIndex:
            i = self.l1[index][0]
            j = self.l1[index][1]
            
            score = logistic(self.X[i].dot(self.mu).dot(self.Y[j].T))
            predictScoreList.append(score)
            trueScoreList.append(1)
          
        for index in l2_testIndex:
            i = self.l2[index][0]
            j = self.l2[index][1]
            
            score = logistic(self.X[i].dot(self.mu).dot(self.Y[j].T))
            predictScoreList.append(score)
            trueScoreList.append(0)
          
        return predictScoreList,trueScoreList
 
    
            
def readFeatures(filepath1,filepath2):
    drugFeature = []
    proteinFeature = []
    
    
    drugFeature = Input.readVectorFromFile(filepath1)
    
    proteinFeature = Input.readVectorFromFile(filepath2)
    
    return np.array(drugFeature),np.array(proteinFeature)



    
def model():
    X1,Y1 = readFeatures("feature/drug_vector_d300.txt","feature/protein_vector_d300.txt")
    
    drugProteinMat = Input.readDrugProteinInteractionMat()
    R = copy.deepcopy(drugProteinMat)
    l1,l2 = getVisitedPart(drugProteinMat)
  
    aveAucList = []
    aveAuprList = []
    for iteration in range(5):
        
        
        F1 = FNML(X1,Y1,drugProteinMat,R,l1,l2)
        
        
        l1_trainIndexList = []
        l2_trainIndexList = []
        l1_testIndexList = []
        l2_testIndexList = []
      
        kf = KFold(n_splits=10,shuffle=True)  
        for train_index , test_index in kf.split(l1): 
            l1_trainIndexList.append(train_index)
            l1_testIndexList.append(test_index)
            
        for train_index , test_index in kf.split(l2):  
            l2_trainIndexList.append(train_index)
            l2_testIndexList.append(test_index)
        
        

        aucList = []
        auprList = []
        
      
        for k in range(len(l1_trainIndexList)):
            
            F1.initParameter()
           
    
            
            l1_trainIndex = l1_trainIndexList[k]
            l2_trainIndex = l2_trainIndexList[k]
            
            l1_testIndex = l1_testIndexList[k]
            l2_testIndex = l2_testIndexList[k]
            
           
            print len(l1),len(l2),len(l1_trainIndex),len(l1_testIndex),len(l2_trainIndex),len(l2_testIndex)
            
            out = open("parameter"+(str)(k)+".txt","w")
            

            lastAUPR = 0
            for i in range(30):
                print k,"iteration:",i
                out.write("iteration "+(str)(i)+":")
                
                F1.trainOnce(l1_trainIndex,l2_trainIndex)
               
               
                predictScoreList1,trueScoreList1 = F1.predictScore(l1_testIndex,l2_testIndex)
                
                auc = metrics.roc_auc_score(trueScoreList1, predictScoreList1)
                
            
                aupr = metrics.average_precision_score(trueScoreList1, predictScoreList1)
                
                
                print "auc:",auc
                print "aupr:",aupr
                out.write("auc:"+(str)(auc)+"\n")
                out.write("aupr:"+(str)(aupr)+"\n")
                
                out.write("\n---------------------------------\n")
                out.flush()
                if lastAUPR>aupr:
                    break
                lastAUPR = aupr
                
            aucList.append(auc)
            auprList.append(aupr)

            
            for auc in aucList:
                print "auc:",auc
            for aupr in auprList:
                print "aupr:",aupr
                
            aveAuc = sum(aucList)/len(aucList)
            aveAupr = sum(auprList)/len(auprList)
           
            print "aveAuc:",aveAuc,"aveAupr:",aveAupr 
            out.write("aveAuc:"+(str)(aveAuc)+" aveAupr:"+(str)(aveAupr)+"\n")
         
            out.flush()
        aveAucList.append(aveAuc)
        aveAuprList.append(aveAupr)
    print "aveAUC:",sum(aveAucList)/len(aveAucList)
    print "aveAUPR:",sum(aveAuprList)/len(aveAuprList)
    

    


if __name__ == '__main__':
    model()       
    