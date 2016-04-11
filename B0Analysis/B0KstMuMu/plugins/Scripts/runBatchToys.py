#################################
# Program to run toy-MC studies #
#################################

from os import system, chdir
import sys

if len(sys.argv) < 4:
    print "Synopsis: python runBatchToy.py type[BF,FLAFB] nBins[0-8,-1] nJobs(per bin)[> 0]"
    sys.exit("ERROR !")

par = sys.argv[1]
if 'BF' in par:
    fitType = 21
elif 'FLAFB' in par:
    fitType = 26
else:
    print "Wrong parameter: ", par
    print "Synopsis: python runBatchToy.py type[BF,FLAFB] nBins[0-8,-1] nJobs(per bin)[> 0]"
    sys.exit("ERROR !")

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

        cmd = "mkdir " + dir
        print cmd
        system(cmd)

        cmd = "cp ../ExtractYield.cc " + dir
        print cmd
        system(cmd)

        cmd = "cp ../ExtractYield " + dir
        print cmd
        system(cmd)

        cmd = "cp ../../python/ParameterFile.txt " + dir
        print cmd
        system(cmd)

        chdir(dir)

        system("sleep 2")
        toRun = "Qsub -l lnxfarm -e -o toyMC_" + par + "_" + str(nJobs) + "_" + str(i) + ".log -N T" + str(nJobs) + par + str(i) + " ./ExtractYield " + str(fitType)
        toRun = toRun + " toyMC_" + par + "_" + str(nJobs) + "_" + str(i) + ".root yesEffCorr " + str(i) + " " + str(nToysJob) + " " + str(nJobs) + " ParameterFile.txt"

        print toRun
        system(toRun)

        chdir("..")

        nJobs = nJobs - 1

chdir("..")
