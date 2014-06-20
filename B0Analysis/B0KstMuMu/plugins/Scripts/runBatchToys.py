#################################
# Program to run toy-MC studies #
#################################

from os import system, chdir
import sys

if len(sys.argv) < 4:
    print "Synopsis: python runBatchToy.py type[BF,FLAFB] nBins[0-8,-1] nJobs[> 0]"
    sys.exit()

par = sys.argv[1]
if par == 'BF':
    fitType = 21
elif par == 'FLAFB':
    fitType = 26
else:
    print "Wrong parameter: ", par
    print "Synopsis: python runBatchToy.py type[BF,FLAFB] nBins[0-8,-1] nJobs[> 0]"

nBins = int(sys.argv[2])
if nBins != -1:
    binList = [nBins,]
else:
    binList = [0,1,2,3,5,7,8]

nJobs    = int(sys.argv[3])
totToys  = 1000
nToysJob = totToys / nJobs

print "Total number of toys: ", totToys
system("unset DISPLAY")

for i in binList:
    nJobs = int(sys.argv[3])
    while nJobs > 0:
        print "\nInstantiating job: ", nJobs, " for bin: ", i

        dir = "./" + par + "_" + str(nJobs) + "_" + str(i) + "/"

        mkdir = "mkdir " + dir
        print mkdir
        system(mkdir)

        cp = "cp ../ExtractYield.cc " + dir
        print cp
        system(cp)

        cp = "cp ../ExtractYield " + dir
        print cp
        system(cp)

        chdir(dir)

        toRun = "Qsub -l lnxfarm -e -o toyMC_" + par + "_" + str(nJobs) + "_" + str(i) + ".log -N T" + str(nJobs) + par + str(i) + " ./ExtractYield " + str(fitType)
        toRun = toRun + " toyMC_" + par + "_" + str(nJobs) + "_" + str(i) + ".root yesEffCorr " + str(i) + " " + str(nToysJob) + " " + str(nJobs) + " ../../python/ParameterFile.txt"
        print toRun
        system(toRun)

        chdir("..")

        nJobs = nJobs - 1
