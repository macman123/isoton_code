# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:44:02 2016

@author: iryna
"""

import pandas as pd
import os
import glob  # to list files with specific names
import numpy as np
from subprocess import call
import math
import paramiko #to execute ssh

#### this script requires transcripts/features file with lines=features columns=patients. A phenotype file with 2 phenotypes coded 1 and 0. 

usrname='macman'
pswd='passSYSBIO'

# Test data set provided?
test_set_provided = True

dirShared ="/home/macman/isotonic_regression/GSE19491_tuberculosis_train_10/" #directory for your analysis


dirIsotonInputs = dirShared + "data/"
patternTranscripts=dirIsotonInputs+"*EXPRESSION*.txt"
patternPhenotypes=dirIsotonInputs+"*PHENOTYPE*"



#fannot = "/home/iryna/isotonic_regression/Dengue/HTA-2_0_probeset_annotations.txt"

isoton="/home/macman/isotonic_regression/201607_isoton/201607_isoton"
mathematica='math' #'/Applications/Mathematica.app/Contents/MacOS/MathKernel'
binarize='writeBinaryData_iryna_17jun16.m' #dirShared+'code/writeBinaryData_iryna_17jun16.m'


dirIsotonOutputs = dirShared + "isotonicResults/"
relativePathIntSetCv="cvInternalSets/"
relativePathTestSetCv="cvTestSet/"
relativePathReal="real/"
relativePathCvIntReal=relativePathReal+relativePathIntSetCv
SuffixCatAndSort="_concatenatedAndSorted/"
dirCvTestSet=dirIsotonInputs+relativePathTestSetCv
PrefixFCvSet="Forisoton_"
PrefixFOutCvSet="Fromisoton_"
SuffixFCvTestSet="_LOTest"
PrefixFPhenoCvSet="Phenotype_"

dirCvIntSet=dirIsotonInputs+relativePathIntSetCv
SuffixFCvIntSet="_LOInt"

dirCvIntReal=dirIsotonInputs+relativePathReal+relativePathIntSetCv


fRealScript=dirIsotonInputs+'real.sh'
fTestSetScript=dirIsotonInputs+'testSet.sh' 
fIntSetScript=dirIsotonInputs+'intSet.sh'    
fRealIntSetScript=dirIsotonInputs+'realIntSet.sh'

nApproxJobsRealData=1000
nApproxJobsCvTestSetData=2000
nApproxJobsCvIntSetData=3000
nApproxJobsCvIntRealData=1000


"""Def functions"""
def matchNamePattern(pattern):
    lFiles=glob.glob(pattern)
    assert (len(lFiles)==1), "The expression pattern " + pattern+ " doesn't match exactely 1 element. It matches "+ str(len(lFiles))+" elements."

    nameFile=lFiles[0]
    return nameFile

def catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fNameForMathema):
    dirOutCat=dirToCat[:(len(dirToCat)-1)]+SuffixCatAndSort
    if not os.path.exists(dirOutCat):
        os.makedirs(dirOutCat)
    
    resFile=dirOutCat+fNameForMathema
    
    bashCommand="cat "+dirToCat+filePatternToCat+" > "+resFile+"_temp1.txt" 
    bashCommand2="sort -nk 1 -nk 2 "+resFile+"_temp1.txt"+" > " +resFile+"_temp2.txt"
    bashCommand3="awk 'NR==1{a[$3]=0;a[$4]=0;print $0;next}{if(!($3 in a || $4 in a)) {print $0; a[$3]=0;a[$4]=0}}' "+resFile+"_temp2.txt > "+resFile+"_temp3.txt"
    bashCommand4="head -1000 "+resFile+"_temp3.txt > "+resFile 
    
    os.system(bashCommand)    
    os.system(bashCommand2)    
    os.system(bashCommand3)
    os.system(bashCommand4)    
        
    os.remove(resFile+"_temp1.txt")
    os.remove(resFile+"_temp2.txt")
    os.remove(resFile+"_temp3.txt")
    return resFile


def createScriptFile(fScript, nFeatures, dirIsotonOutputs, relativePathSetCv, PrefixFCvSet, SuffixFCvSet, dirCvSet, PrefixFPhenoCvSet,nApproxJobsCvSetData):
    dirOutCv=dirIsotonOutputs+relativePathSetCv
    
    if not os.path.exists(dirIsotonOutputs):
            os.makedirs(dirIsotonOutputs)
    if not os.path.exists(dirOutCv):
            os.makedirs(dirOutCv)
    SetFiles=glob.glob(dirCvSet+PrefixFCvSet+"*"+SuffixFCvSet+".txt.bin")
    SetFiles.sort()
    
    
    SetPheno=glob.glob(dirCvSet+PrefixFPhenoCvSet+"*"+SuffixFCvSet+".txt")
    SetPheno.sort()
    
    
    nSetFiles=int(len(SetFiles))
    nPhenoSetFiles=int(len(SetPheno))
    assert (nSetFiles==nPhenoSetFiles), "Not the same number of phenotype and gene expression files in set directory. nPhenotypes="+str(nPhenoSetFiles)+"nSetFiles="+str(nSetFiles)+". Files are defined as: "+dirCvSet+PrefixFCvSet+"*"+SuffixFCvSet+".bin"
    
    nApproxJobsPerCvSetFile=math.ceil(nApproxJobsCvSetData/nSetFiles)
    
    
     # determining the parameter b defined such as 2^b is nFeaturesPerJob  in isoton script:
    bIsoton = math.ceil(math.log(nFeatures * (1 + math.sqrt(1 + 8 * nApproxJobsPerCvSetFile)) / (4 * nApproxJobsPerCvSetFile), 2))
        
    mIsoton = math.ceil(nFeatures / math.pow(2, bIsoton))
        
        #  nJobsRealData corresponds to the parameter we called "s" in isoton. It is equal to half of a square of length mIsoton with its diagonal)
    nJobsPerSetFile=int(((mIsoton*mIsoton +mIsoton)/2))    
    nJobsSetData=int(nJobsPerSetFile*nSetFiles)
        
    
    if os.path.isfile(fScript):
        os.remove(fScript)
           
    with open(fScript, 'a') as the_file:
        the_file.write("#!/bin/bash\n#$ -cwd \n#$ -S /bin/bash\n#$ -o /dev/null\n#$ -t 1-"+str(math.floor(nJobsSetData))+"\n\n")
        the_file.write("files=('"+"'\t'".join(SetFiles)+"')\n\n")
        the_file.write("phenotypes=('"+"'\t'".join(SetPheno)+"')\n\n")
        the_file.write("nJobsPerSetFile="+str(nJobsPerSetFile)+"\n")
        
    #file that we're deling with    
        the_file.write("nFile=$[($SGE_TASK_ID-1)/$nJobsPerSetFile]"+"\n")
    #subjob number for that file
        the_file.write("nSubJob=$[$SGE_TASK_ID-$nFile*$nJobsPerSetFile]"+"\n")
        the_file.write("file=${files[$nFile]}"+"\n")    
        the_file.write("'"+isoton+"' --exhaustiveCV -k 1  -q 2000 -b "+str(bIsoton)+" -t $nSubJob " + "$file ${phenotypes[$nFile]} > '"+dirOutCv+"'$SGE_TASK_ID."+"${file##*/}"+"'.'"+"$nSubJob"+"'.txt'"+" \n")    
    return
    
    
    




def runScriptFileOnSpice(usrname,pswd,fScript):

###ssh spice and submit to queue the fScript file analysis

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect('spice', username = usrname , password = pswd)
    
    
    
    #With -sync -y , I'm waiting for all jobs to end before the command finishes.
    
    stdin, stdout, stderr= ssh.exec_command('source /etc/profile; qsub -sync y '+fScript)
    exit_status = stdout.channel.recv_exit_status()
    stdout_output = stdout.read().decode('utf8').rstrip('\n')
    print ("Stdout: ", stdout_output)
    stderr_output = stderr.read().decode('utf8').rstrip('\n')
    print ("Stderr: ", stderr_output)
    assert (exit_status == 0), "While executing "+fScript+", qsub generated the non-zero exit status: "+str(exit_status)
    return exit_status



def makeCvFilesLeaveOneOut(dftranscripts, dfphenotypes,dirOut, suffixOut,PrefixFExprCv,PrefixFPhenoCv):
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)
    nColsdfTranscripts=len(dftranscripts.columns)
    for i in range (0, nColsdfTranscripts):
        keptCols=list(range(0,i))+list(range(i+1,nColsdfTranscripts))
        dfKeptCols = dftranscripts.ix[:,keptCols]
        dfPhenoKeptCols = dfphenotypes.iloc[keptCols]
        fOut=dirOut+PrefixFExprCv+str(i)+suffixOut+".txt"
        fPhenoOut=dirOut+PrefixFPhenoCv+str(i)+suffixOut+".txt"
        dfKeptCols.to_csv( fOut, sep="\t", float_format='%.6f')
        dfPhenoKeptCols.columns=[str(int(len(dfPhenoKeptCols)))]
        dfPhenoKeptCols.to_csv( fPhenoOut, float_format='%.0f', index=False, sep="\t")
        
        
def makeCvFilesLeaveThreeOut(dftranscripts, dfphenotypes,dirOut, suffixOut,PrefixFExprCv,PrefixFPhenoCv,lim):
    if not os.path.exists(dirOut):
        os.makedirs(dirOut)
        
    nCasesThisSet = int(dfphenotypes.sum())
    nControlsThisSet = int(len(dfphenotypes)-nCasesThisSet)
    
    nColsdfTranscripts=len(dftranscripts.columns)
    for i in range (0, lim):
        #if I need to change leave out patients, change here:
        keptCols=list(range(0,i))+list(range(i+1,i+nCasesThisSet))+list(range(i+nCasesThisSet+1,i+nCasesThisSet+lim))+list(range(i+nCasesThisSet+lim+1,nCasesThisSet+nControlsThisSet))
        dfKeptCols = dftranscripts.ix[:,keptCols]
        dfPhenoKeptCols = dfphenotypes.iloc[keptCols]
        fOut=dirOut+PrefixFExprCv+str(i)+"_"+str(i+nCasesThisSet)+"_"+str(i+nCasesThisSet+lim)+suffixOut+".txt"
        fPhenoOut=dirOut+PrefixFPhenoCv+str(i)+"_"+str(i+nCasesThisSet)+"_"+str(i+nCasesThisSet+lim)+suffixOut+".txt"
        dfKeptCols.to_csv( fOut, sep="\t", float_format='%.6f')
        dfPhenoKeptCols.columns=[str(int(len(dfPhenoKeptCols)))]
        dfPhenoKeptCols.to_csv( fPhenoOut, float_format='%.0f', index=False, sep="\t")
    
##def vars
ftranscripts = matchNamePattern(patternTranscripts)
nameftranscripts=ftranscripts.rsplit('/', 1)[1]
fphenotypes=matchNamePattern(patternPhenotypes)

dfphenotypes = pd.read_csv(fphenotypes, skiprows=1, names=["Phenotype"])
nCases = int(dfphenotypes.sum())
nControls = int(len(dfphenotypes)-nCases)

dftranscripts = pd.read_csv(ftranscripts, sep="\t", index_col=0)
nFeatures = len(dftranscripts)
nColsdfTranscripts = len(dftranscripts.columns)
realCols = ['s1','s2','t1','t2','dir1','dir2']
lim =int(min((nControls-1) / 2, nCases-1))# len(dftranscripts.columns)-1 
limRealIntCv=int(min((nControls-1) / 2, nCases-1))#len(dftranscripts.columns)
assert (len(dftranscripts.columns)==len(dfphenotypes)), "We do not have"+str(len(dftranscripts.columns))==str(len(dfphenotypes))
fAnnotOutRealOnce=dirIsotonOutputs+nameftranscripts+"_isotonResult_annot.txt"

###Making files for cv test set: leaving 1 element out in turn

assert (len(dftranscripts.columns) == nCases+nControls), "nCases + nControls doesn't correspont to the length of columns of dftranscripts. "
if test_set_provided == False:
    makeCvFilesLeaveOneOut (dftranscripts, dfphenotypes,dirIsotonInputs+relativePathTestSetCv, SuffixFCvTestSet,PrefixFCvSet,PrefixFPhenoCvSet)

print("files for cv test set generated")

####Making files for internal cv (to determine k): leaving 1 element out in turn for each of the previous cvs
if test_set_provided == False:
    patternCvTestSetInput=dirIsotonInputs+relativePathTestSetCv+PrefixFCvSet+"*"+SuffixFCvTestSet+".txt"

    for j in range(0,len(glob.glob(patternCvTestSetInput))):
        fCvTestSetInput=dirIsotonInputs+relativePathTestSetCv+PrefixFCvSet+str(j)+SuffixFCvTestSet+".txt"
    
        fCvTestSetPhenoInput=dirIsotonInputs+relativePathTestSetCv+PrefixFPhenoCvSet+str(j)+SuffixFCvTestSet+".txt"
        dfCvTestSetInput=pd.read_csv(fCvTestSetInput, sep="\t", index_col=0)
        dfCvTestSetPhenoInput=pd.read_csv(fCvTestSetPhenoInput,skiprows=1, names=["Phenotype"])
        prefixExprIntCv=PrefixFCvSet+str(j)+SuffixFCvTestSet+"_"
        prefixPhenoIntCv=PrefixFPhenoCvSet+str(j)+SuffixFCvTestSet+"_"
    
        makeCvFilesLeaveThreeOut (dfCvTestSetInput, dfCvTestSetPhenoInput,dirIsotonInputs+relativePathIntSetCv, SuffixFCvIntSet,prefixExprIntCv,prefixPhenoIntCv, lim)
        print (j)
        
else:
    makeCvFilesLeaveOneOut (dftranscripts, dfphenotypes,dirIsotonInputs+relativePathIntSetCv, SuffixFCvIntSet,PrefixFCvSet,PrefixFPhenoCvSet)


   




####Making files for int  cv real data (no left out patient for testing)

if not os.path.exists(dirIsotonInputs+relativePathReal):
        os.makedirs(dirIsotonInputs+relativePathReal)

if test_set_provided == False:
    makeCvFilesLeaveThreeOut (dftranscripts, dfphenotypes, dirIsotonInputs+relativePathCvIntReal, SuffixFCvIntSet,PrefixFCvSet,PrefixFPhenoCvSet,lim)
else:
    if not os.path.exists(dirIsotonInputs+relativePathCvIntReal):
        os.makedirs(dirIsotonInputs+relativePathCvIntReal)
    dftranscripts.to_csv(dirIsotonInputs+relativePathCvIntReal+PrefixFCvSet+"complete"+SuffixFCvIntSet+".txt", float_format='%.0f', index=False, sep="\t")
    dfphenotypes.to_csv(dirIsotonInputs+relativePathCvIntReal+PrefixFPhenoCvSet+"complete"+SuffixFCvIntSet+".txt", float_format='%.0f', index=False, sep="\t")
"""
call([mathematica, '-script', binarize, dirIsotonInputs+relativePathCvIntReal,PrefixFCvSet+"*"+".txt" ])

call([mathematica, '-script', binarize, dirIsotonInputs,nameftranscripts ])


#Making files for misclassification errors' error bars
#
#
#
#End of Making files for misclassification errors' error bars  
#
#
#Call mathematica script: binarize files for cross validation
# 
#
call([mathematica, '-script', binarize, dirCvIntSet,PrefixFCvSet+"*"+".txt" ])
if test_set_provided == False:
    call([mathematica, '-script', binarize, dirCvTestSet,PrefixFCvSet+"*"+".txt" ])
else:
    call([mathematica, '-script', binarize, dirCvTestSet,"*.txt" ])
#
#
#

####End of binarize



#####create real.sh file for cluster:
call([mathematica , '-script', binarize, dirIsotonInputs, nameftranscripts ]) # '/Users/benno/dev/mmatest/test.ma'])


if not os.path.exists(dirIsotonOutputs):
    os.makedirs(dirIsotonOutputs)
    
if not os.path.exists(dirIsotonOutputs+relativePathReal ):
    os.makedirs(dirIsotonOutputs+relativePathReal)
    
# determining the parameter b defined such as 2^b is nFeaturesPerJob  in isoton script:
bIsoton = math.ceil(math.log(nFeatures * (1 + math.sqrt(1 + 8 * nApproxJobsRealData)) / (4 * nApproxJobsRealData), 2))
    
mIsoton = math.ceil(nFeatures / math.pow(2, bIsoton))
    
#  nJobsRealData corresponds to the parameter we called "s" in isoton. It is equal to half of a square of length mIsoton with its diagonal)
nJobsRealData=(mIsoton*mIsoton +mIsoton)/2
        

if os.path.isfile(fRealScript):
    os.remove(fRealScript)

with open(fRealScript, 'a') as the_file:
    the_file.write("#!/bin/bash\n#$ -cwd \n#$ -S /bin/bash\n#$ -o /dev/null\n#$ -t 1-"+str(math.floor(nJobsRealData))+"\n\n")
    the_file.write(isoton+" --exhaustiveCV -k 1  -q 2000 -b "+str(bIsoton)+" -t $SGE_TASK_ID " + ftranscripts+".bin "+ fphenotypes+" > "+dirIsotonOutputs+relativePathReal+"$SGE_TASK_ID.realResults.txt \n")    
    


    
### run  real.sh 


####ssh spice and submit to queue the real file analysis
runScriptFileOnSpice(usrname,pswd,fRealScript)






#### create and run cvTestSet.sh file for cluster
 
createScriptFile(fTestSetScript, nFeatures, dirIsotonOutputs, relativePathTestSetCv, PrefixFCvSet, SuffixFCvTestSet, dirCvTestSet, PrefixFPhenoCvSet, nApproxJobsCvTestSetData)     
runScriptFileOnSpice(usrname,pswd,fTestSetScript)




####create and run fIntSetScript file for the cluster

createScriptFile(fIntSetScript, nFeatures, dirIsotonOutputs, relativePathIntSetCv, PrefixFCvSet, SuffixFCvIntSet, dirCvIntSet, PrefixFPhenoCvSet, nApproxJobsCvIntSetData)     
runScriptFileOnSpice(usrname,pswd,fIntSetScript)

####create and run fRealIntSetScript file for the cluster
createScriptFile(fRealIntSetScript, nFeatures, dirIsotonOutputs, relativePathReal+relativePathIntSetCv, PrefixFCvSet, SuffixFCvIntSet, dirCvIntReal, PrefixFPhenoCvSet, nApproxJobsCvIntRealData)     

runScriptFileOnSpice(usrname,pswd,fRealIntSetScript)



##########Preparing real file for mathematica



dirToCat=dirIsotonOutputs+relativePathReal
filePatternToCat="*.txt"

fForMathema=PrefixFOutCvSet+"all_once.txt"

fForMathemaReal=catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fForMathema)

print("real file for mathematica: ready") 


    
##########Preparing Test Set files for mathematica 

for i in range (0, nColsdfTranscripts):
    subpattern=str(i)+SuffixFCvTestSet
    filePatternToCat="*"+PrefixFCvSet+subpattern+"*.txt"  
    dirToCat=dirIsotonOutputs+relativePathTestSetCv
    
    fForMathema=PrefixFOutCvSet+subpattern+"_once.txt"
    
    catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fForMathema)

print("Test Set files for mathematica: ready") 



#####preparing Internal Set files for mathematica:

dirToCat=dirIsotonOutputs+relativePathIntSetCv


for j in range (0,nColsdfTranscripts):
    
    for i in range(0, lim):
        subPattern=str(j)+SuffixFCvTestSet+"_"+str(i)+SuffixFCvIntSet
        filePatternToCat="*"+PrefixFCvSet+subPattern+"*.txt"
        fForMathema=PrefixFOutCvSet+subPattern+"_once.txt"
        catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fForMathema)
        
        
print("Internal Set files for mathematica: ready")        
        
        
        
        
        
#        
#        with open(catFile, 'w') as outfile:
#            for fname in filenames:
#                with open(fname) as infile:
#                    for line in infile:
#                        outfile.write(line)
#        
#        bashCommand="cat "+catFile+"  | sort -nk 1 -nk 2 | awk 'NR==1{a[$3]=0;a[$4]=0;print $0;next}{if(!($3 in a || $4 in a)) {print $0; a[$3]=0;a[$4]=0}}' | head -1000 > "+ fForMathema 
#        os.system(bashCommand)
#        os.remove(catFile)
#

 


###Preparing real intSet files for mathematica


dirToCat=dirIsotonOutputs+relativePathReal+relativePathIntSetCv

      
for i in range(0, limRealIntCv):
    subPattern="_"+str(i)+SuffixFCvIntSet
    filePatternToCat="*"+PrefixFCvSet+subPattern+"*.txt"

    patternfOut=PrefixFOutCvSet+subPattern
    
    
    fForMathema=patternfOut+"_once.txt"
    catSortOnce(dirToCat,SuffixCatAndSort,filePatternToCat,fForMathema)
        
print("Real Internal Set files for mathematica: ready")        
    
    
    
    
#    
#    with open(catFile, 'w') as outfile:
#        for fname in filenames:
#            with open(fname) as infile:
#                for line in infile:
#                    outfile.write(line)
#    
#    bashCommand="cat "+catFile+"  | sort -nk 1 -nk 2 | awk 'NR==1{a[$3]=0;a[$4]=0;print $0;next}{if(!($3 in a || $4 in a)) {print $0; a[$3]=0;a[$4]=0}}' | head -1000 > "+ fForMathema 
#    os.system(bashCommand)
#    os.remove(catFile)
#
#
#

#####Annotating real file: not finished


dfOutRealOnce=pd.read_csv("/home/iryna/isotonic_regression/Dengue/20160817_2009Classif/isotonicResults/real_concatenatedAndSorted/Fromisoton_all_once.txt", sep=" ", header=None, names=realCols) #fForMathemaReal
#dfOutRealOnce=pd.read_csv(, sep=" ", header=None, names=realCols)

def annotate(dfreal, ftranscripts, fannot, outannotreal):#annotate real result:
    
    dftranscripts= pd.read_csv(ftranscripts, sep="\t", header=0, index_col=0)
    


    dfannot=pd.read_csv(fannot, sep="\t", header=0)#, index_col='Probeset_ID')
    
    dfreal['t1ProbeID']=dftranscripts.index.values[dfreal['t1']]
    dfreal['t2ProbeID']=dftranscripts.index.values[dfreal['t2']]
    
    
    colDfAnnotSaved=dfannot.columns
    dfannot.columns =['t1ProbeID','t1EntrezID','t1GeneSymbol']
    dfreal=pd.merge(dfreal,dfannot, on=['t1ProbeID'])
    dfannot.columns =['t2ProbeID','t2EntrezID','t2GeneSymbol']
    dfreal=pd.merge(dfreal,dfannot, on=['t2ProbeID'])
    dfannot.columns=colDfAnnotSaved         
    
    
    dfreal.to_csv(outannotreal, sep="\t", columns=['s1','s2','t1GeneSymbol', 't2GeneSymbol',
    't1ProbeID', 't2ProbeID','t1','t2','t1EntrezID', 't2EntrezID','dir1','dir2'])      
    return       
 
annotate( dfOutRealOnce, ftranscripts, fannot, fAnnotOutRealOnce)




### deleting the useless partners
dfOutRealOnce=dfOutRealOnce=pd.read_csv("/home/iryna/isotonic_regression/Dengue/20160817_2009Classif/isotonicResults/real_concatenatedAndSorted/Fromisoton_all_once.txt", sep=" ", header=None, names=realCols)
k=10

dfToBeFixed=pd.concat([dfOutRealOnce.ix[0:(k-1),2], dfOutRealOnce.ix[0:(k-1),3]], axis=0)

dfToBeFixed.to_csv(dirIsotonOutputs+"ToBeFixed.txt", index=False)



"""   
