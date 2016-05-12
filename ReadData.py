#!/usr/bin/env python

# Originally written by Bryan Clark
# modified by ChangMo Yang to work with Python 3

from numpy import *
import re

def removeBlanksFromList(theList):
    newList=[]
    for counter in range(len(theList)):
        if theList[counter]!='':
            newList.append(theList[counter])
    return newList


def load(fileName):
    loadAscii(fileName)


#def saveAscii(theData,fileName):
#    fileOut=open(fileName,'w')
#    for counter in range(len(theData)):
#        for counter2 in range(len(theData[counter])):
#            print(theData[counter,counter2],end='',file=fileOut)
#        print(file=fileOut)
#    fileOut.close()


#Reads in an ascii file of data much like octave
def loadAscii(fileName):
    infile=open(fileName)
    theLines=infile.readlines()
    infile.close()
    splitter=re.compile('\s|,')
    lineSize=len(removeBlanksFromList(splitter.split(theLines[0])))
    numLines=len(theLines)
    if lineSize > 1:
        newData=zeros((numLines,lineSize))
    else:
        newData=zeros(numLines)
    for counter in range(numLines):
        theLine=removeBlanksFromList(splitter.split(theLines[counter]))
        if len(theLine)!=lineSize:
            print("All lines are not the same size in",fileName)
            #abort()
        else:
            if lineSize>1:
                for counter2 in range(lineSize):
                    newData[counter,counter2]=float(theLine[counter2])
            else:
                for counter2 in range(lineSize):
                    newData[counter]=float(theLine[counter2])
    return newData


#Prints a 2d real array of data
#def Print2d(theData):
#    for counter in range(len(theData)):
#        for counter2 in range(len(theData[counter])):
#            print(theData[counter,counter2],end='')
#        print()


#take a set of data and if there are multiple
#x values that are the same, average them together
def HistogramAverageData(theData):
    newData=[]
    theData=theData.tolist()
    theData.sort()
    theData=array(theData)
    counter=0
    (lenX,lenY)=shape(theData)
    while (counter<lenX):
      totalY=theData[counter,1]
      totalError=theData[counter,2]
      numY=1
      counter=counter+1
      while (counter<lenX and theData[counter,0]-theData[counter-1,0]<1e-5):
           totalY=totalY+theData[counter,1]
           totalError=totalError+theData[counter,2]
           numY=numY+1
           counter=counter+1
      newData.append([theData[counter-1,0],totalY/(numY+0.0),totalError/(math.sqrt(numY+0.0))])
#      if len(theData)==3:
 #         newData.append([theData[counter-1,0],totalY/(numY+0.0),totalError/(numY+0.0)])
 #     else:
 #         newData.append([theData[counter-1,0],totalY/(numY+0.0)])
#Removed HACK!      newData.append([theData[counter-1,0],totalY/(1+0.0)])
    return array(newData)


#myData=load('borisData')
#histData=HistogramAverageData(myData)
#Print2d(myData)
