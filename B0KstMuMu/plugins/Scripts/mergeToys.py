###################################################################
# Program to merge root files generated with ExtractYield.cc toys #
# To be used after running the program runBatchToys.py            #
###################################################################

from os import system, chdir
import sys
import glob
import string

if len(sys.argv) < 4:
    print "Synopsis: python mergeToys.py [q2bin nJobs dirName(\"name_nJobs_q2bin\" without \"_nJobs_q2bin\")]"
    sys.exit("ERROR !")

q2bin   = sys.argv[1]
nJobs   = int(sys.argv[2])
dirName = sys.argv[3]

print "\n@@@ Input parameters : ", q2bin, nJobs, dirName, " @@@"

#################################
# Read the number of root files #
#################################
newDir = dirName + "_1_" + q2bin
chdir(newDir)
nFiles = len(glob.glob("*.root"))
chdir("..")

fileList = [["" for x in xrange(nFiles)] for x in xrange(nJobs)]


################################
# Load all the root file names #
################################
for j in range(0, nJobs):
    newDir = dirName + "_" + str(j+1) + "_" + q2bin
    chdir(newDir)
    fileList[j] = glob.glob("*.root")
    print fileList[j]
    chdir("..")


#########################
# Add all the root file #
#########################
for i in range(0, nFiles):
    outputFile = fileList[0][i]
    outputFile = outputFile.replace("_1_","_",1)
    cmd = "hadd " + outputFile + " "
    
    for j in range(0, nJobs):
        if len(fileList[j]) != 0:
            cmd += dirName + "_" + str(j+1) + "_" + q2bin + "/" + fileList[j][i] + " "
    print cmd
    system(cmd)
