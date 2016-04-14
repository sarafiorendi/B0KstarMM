"""
##################################################
MVA implementation with Perceptron Neural Networks
                                  by Mauro Dinardo
##################################################
"""
####################
# @TMP@            #
# - Inject noise   #
# - Batch learning #
####################
from argparse  import ArgumentParser
from random    import seed, random

from ROOT      import gROOT, gStyle, TCanvas, TGraph, TGaxis

from NeuralNet import NeuralNet


def ArgParser():
    ###################
    # Argument parser #
    ###################
    parser = ArgumentParser()
    parser.add_argument("-nv", "--Nvars",        dest = "Nvars",        type = int,  help = "Number of variables",              required=True)
    parser.add_argument("-np", "--Nperceptrons", dest = "Nperceptrons", type = int,  help = "Number of perceptrons",            required=True)
    parser.add_argument("-nn", "--Nneurons",     dest = "Nneurons",     type = int,  help = "Number of neurons per perceptron", required=True,  nargs='*')

    options = parser.parse_args()
    if options.Nvars:
        print "I'm reading the variable number: ", options.Nvars

    if options.Nperceptrons:
        print "I'm reading the perceptron number: ", options.Nperceptrons

    if options.Nneurons:
        print "I'm reading the neuron number per perceptron: ", options.Nneurons

    return options


def SetStyle():
    #######################
    # ROOT style settings #
    #######################
    gROOT.SetStyle("Plain")
    gROOT.ForceStyle()
    gStyle.SetTextFont(42)
    
    gStyle.SetOptTitle(0)
    gStyle.SetOptFit(1112)
    gStyle.SetOptStat(1110)
    
    gStyle.SetPadRightMargin(0.08)
    gStyle.SetPadTopMargin(0.11)
    gStyle.SetPadBottomMargin(0.12)
    
    gStyle.SetTitleFont(42,"x")
    gStyle.SetTitleFont(42,"y")
    gStyle.SetTitleFont(42,"z")
    
    gStyle.SetTitleOffset(1.05,"x")
    gStyle.SetTitleOffset(0.95,"y")
    
    gStyle.SetTitleSize(0.05,"x")
    gStyle.SetTitleSize(0.05,"y")
    gStyle.SetTitleSize(0.05,"z")
    
    gStyle.SetLabelFont(42,"x")
    gStyle.SetLabelFont(42,"y")
    gStyle.SetLabelFont(42,"z")
    
    gStyle.SetLabelSize(0.05,"x")
    gStyle.SetLabelSize(0.05,"y")
    gStyle.SetLabelSize(0.05,"z")
    
    TGaxis.SetMaxDigits(3)
    gStyle.SetStatY(0.9)


def findLevelCurveNN(NN,graphNN,xstart,xend,ystart,yend,nbins):
    for i in xrange(nbins):
        myMin = 1
        myX   = i * (xend-xstart) / nbins + xstart
        found = False
        for j in xrange(nbins):
            y = j * (yend-ystart) / nbins + ystart
            NNoutput = NN.eval([myX,y])
            if abs(NNoutput[0]) < myMin:
                myMin = abs(NNoutput[0])
                myY   = y
                found = True
        if found is True:
            graphNN.SetPoint(i,myX,myY)


################
# Main program #
################
cmd = ArgParser()


#############################
# Graphics layout and plots #
#############################
gROOT.Reset()
SetStyle()

cCost  = TCanvas('cCost', 'The Cost Function', 800, 0, 700, 500)
cSpace = TCanvas('cSpace', 'The Parameter Space', 0, 0, 700, 500)

graphNNCost = TGraph()
graphNNCost.SetTitle('NN cost function;Epoch [#];Cost Function')
graphNNCost.SetMarkerStyle(20)
graphNNCost.SetMarkerSize(0.5)
graphNNCost.SetMarkerColor(1)

graphS = TGraph()
graphS.SetTitle('NN test in 2D;x;y')
graphS.SetMarkerStyle(20)
graphS.SetMarkerSize(0.5)
graphS.SetMarkerColor(4)

graphB = TGraph()
graphB.SetTitle('NN test in 2D;x;y')
graphB.SetMarkerStyle(20)
graphB.SetMarkerSize(0.5)
graphB.SetMarkerColor(2)

graphNN = TGraph()
graphNN.SetTitle('NN test in 2D;x;y')
graphNN.SetMarkerStyle(20)
graphNN.SetMarkerSize(0.7)
graphNN.SetMarkerColor(1)


###################################
# Use always the same random seed #
###################################
seed(0)
nruns    = 200000
scrstart =  10000
scrlen   =  10000


##############################
# Neural net: initialization #
##############################
NN = NeuralNet(cmd.Nvars,cmd.Nperceptrons,cmd.Nneurons)
NN.printParams()




xstart = -4.
xend   =  4.
ystart = -4.
yend   =  4.
allneurons = sum(cmd.Nneurons)
########################
# Neural net: training #
########################
for i in xrange(nruns):
    x = random() * (xend-xstart) + xstart
    y = random() * (yend-ystart) + ystart

    if 6*x*x*x > y:
        target = +1
        graphS.SetPoint(i,x,y)
    else:
        target = -1
        graphB.SetPoint(i,x,y)

    # @TMP@
    if i > scrstart and (i-scrstart) % scrlen == 0 and (i-scrstart)/scrlen-1 < allneurons:
        if (i-scrstart)/scrlen-1 == 0:
            NN.scramble({0:[0]})
            print "AAA : 0"
        elif (i-scrstart)/scrlen-1 == 1:
            NN.scramble({0:[1]})
            print "AAA : 1"

        elif (i-scrstart)/scrlen-1 == 2:
            NN.scramble({1:[0]})
            print "AAA : 2"
        elif (i-scrstart)/scrlen-1 == 3:
            NN.scramble({1:[1]})
            print "AAA : 3"
        elif (i-scrstart)/scrlen-1 == 4:
            NN.scramble({1:[2]})
            print "AAA : 4"
        elif (i-scrstart)/scrlen-1 == 5:
            NN.scramble({1:[3]})
            print "AAA : 5"

        elif (i-scrstart)/scrlen-1 == 6:
            NN.scramble({2:[0]})
            print "AAA : 6"
        elif (i-scrstart)/scrlen-1 == 7:
            NN.scramble({2:[1]})
            print "AAA : 7"
        elif (i-scrstart)/scrlen-1 == 8:
            NN.scramble({2:[2]})
            print "AAA : 8"
        elif (i-scrstart)/scrlen-1 == 9:
            NN.scramble({2:[3]})
            print "AAA : 9"

        elif (i-scrstart)/scrlen-1 == 10:
            NN.scramble({3:[0]})
            print "AAA : 10"
        elif (i-scrstart)/scrlen-1 == 11:
            NN.scramble({3:[1]})
            print "AAA : 11"

        elif (i-scrstart)/scrlen-1 == 12:
            NN.scramble({4:[0]})
            print "AAA : 12"

    cost = NN.learn([x,y],[target])
    graphNNCost.SetPoint(i,i,cost)
NN.printParams()
NN.save("NeuralNet.txt")


####################
# Neural net: test #
####################
count = 0.
for i in xrange(nruns):
    x = random() * (xend-xstart) + xstart
    y = random() * (yend-ystart) + ystart
    NNoutput = NN.eval([x,y])
    if ((NNoutput[0] > 0 and 6*x*x*x > y) or (NNoutput[0] < 0 and 6*x*x*x < y)):
        count += 1
print "\n--> I got it right:", count / nruns * 100., "% <--"




####################
# Neural net: show #
####################
findLevelCurveNN(NN,graphNN,xstart,xend,ystart,yend,100)


##############
# Show plots #
##############
cCost.cd()
graphNNCost.Draw('AL')
cCost.Modified()
cCost.Update()

cSpace.cd()
cSpace.DrawFrame(xstart - 0.5,ystart - 0.5,xend + 0.5,yend + 0.5);
graphS.Draw('P')
graphB.Draw('P same')
graphNN.Draw('P same')
cSpace.Modified()
cSpace.Update()


############################
# Wait for keyborad stroke #
############################
raw_input("\n---press enter---")
