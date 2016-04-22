"""
###############################################
MVA implementation with Neural Networks
                       by Mauro Dinardo
###############################################
Before running double-check:
  - number of perceptrons & neurons
  - learnRate
  - scramble
To-do:
  - batch learning
  - porting in pyCUDA
###############################################
e.g.: python MVA.py -nv 2 -np 6 -nn 2 3 4 3 2 1
###############################################
"""
from argparse  import ArgumentParser
from random    import seed, random, gauss

from ROOT      import gROOT, gStyle, TCanvas, TGraph, TH1D, TGaxis, TLegend

from NeuralNet import NeuralNet


def ArgParser():
    """
    ###############
    Argument parser
    ###############
    """
    parser = ArgumentParser()
    parser.add_argument("-nv", "--Nvars",        dest = "Nvars",        type = int,  help = "Number of variables",              required=True)
    parser.add_argument("-np", "--Nperceptrons", dest = "Nperceptrons", type = int,  help = "Number of perceptrons",            required=True)
    parser.add_argument("-nn", "--Nneurons",     dest = "Nneurons",     type = int,  help = "Number of neurons per perceptron", required=True,  nargs='*')
    parser.add_argument("-sc", "--Scramble",     dest = "doScramble",   type = bool, help = "Perform NN scrambling",            required=False, default=False)

    options = parser.parse_args()

    print ""
    if options.Nvars:
        print "--> I'm reading the variable number: ", options.Nvars

    if options.Nperceptrons:
        print "--> I'm reading the perceptron number: ", options.Nperceptrons

    if options.Nneurons:
        print "--> I'm reading the neuron number per perceptron: ", options.Nneurons

    if options.doScramble:
        print "--> I'm reading the option scramble: ", options.doScramble

    return options


def SetStyle():
    """
    ###################
    ROOT style settings
    ###################
    """
    gROOT.SetStyle("Plain")
    gROOT.ForceStyle()
    gStyle.SetTextFont(42)
    
    gStyle.SetOptTitle(0)
    gStyle.SetOptFit(0)
    gStyle.SetOptStat(0)
    
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


    
    
"""
############
Main program
############
"""
cmd = ArgParser()


"""
#########################
Graphics layout and plots
#########################
"""
gROOT.Reset()
SetStyle()

cCost  = TCanvas('cCost',  'NN Cost Function', 0, 0, 700, 500)
cNNin  = TCanvas('cNNin',  'NN Input',         0, 0, 700, 500)
cNNout = TCanvas('cNNout', 'NN Output',        0, 0, 700, 500)
cNNval = TCanvas('cNNval', 'NN Values',        0, 0, 700, 500)
cDeriv = TCanvas('cDeriv', 'NN Derivatives',   0, 0, 700, 500)
cDelta = TCanvas('cDelta', 'NN Deltas',        0, 0, 700, 500)

graphNNCost = TGraph()
graphNNCost.SetTitle('NN cost function;Epoch [#];Cost Function')

graphSin = TGraph()
graphSin.SetTitle('NN input;x;y')
graphSin.SetMarkerStyle(20)
graphSin.SetMarkerSize(0.5)
graphSin.SetMarkerColor(4)

graphBin = TGraph()
graphBin.SetTitle('NN input;x;y')
graphBin.SetMarkerStyle(20)
graphBin.SetMarkerSize(0.5)
graphBin.SetMarkerColor(2)

graphSout = TGraph()
graphSout.SetTitle('NN output;x;y')
graphSout.SetMarkerStyle(20)
graphSout.SetMarkerSize(0.5)
graphSout.SetMarkerColor(4)

graphBout = TGraph()
graphBout.SetTitle('NN output;x;y')
graphBout.SetMarkerStyle(20)
graphBout.SetMarkerSize(0.5)
graphBout.SetMarkerColor(2)

graphNNerr = TGraph()
graphNNerr.SetTitle('NN output;x;y')
graphNNerr.SetMarkerStyle(20)
graphNNerr.SetMarkerSize(0.5)
graphNNerr.SetMarkerColor(3)

histoNNS = TH1D('histoNNS','histoNNS',100,-1,1)
histoNNS.SetTitle('NN signal output;NN output;Entries [#]')
histoNNS.SetLineColor(4)

histoNNB = TH1D('histoNNB','histoNNB',100,-1,1)
histoNNB.SetTitle('NN background output;NN output;Entries [#]')
histoNNB.SetLineColor(2)

legDeriv = TLegend(0.8, 0.12, 0.92, 0.89, "")
legDeriv.SetTextSize(0.03)
legDeriv.SetFillStyle(0)

legDelta = TLegend(0.8, 0.12, 0.92, 0.89, "")
legDelta.SetTextSize(0.03)
legDelta.SetFillStyle(0)

graphDeriv = []
graphDelta = []


"""
###################
Internal parameters
###################
"""
seed(0)
nRuns     = 10000000
scrStart  =   100000
scrLen    =   100000
saveEvery =     1000
xRange    = 3.
xOffset   = 3.
yRange    = 3.
yOffset   = 0.
noiseBand = 0.1
loR       = 0.5
hiR       = 1.
xyCorr    = lambda x,y: ((x-xOffset)*(x-xOffset)+(y-yOffset)*(y-yOffset))


"""
##########################
Neural net: initialization
##########################
"""
seed(0)
NN = NeuralNet(cmd.Nvars,cmd.Nperceptrons,cmd.Nneurons)
NN.printParams()


"""
####################
Neural net: training
####################
"""
print "\n=== Training neural netowrk ==="
for n in xrange(nRuns):
    x = random() * xRange + xOffset - xRange/2
    y = random() * yRange + yOffset - yRange/2

    if gauss(loR,noiseBand) <= xyCorr(x,y) < gauss(hiR,noiseBand):
        target = +1
        if n % saveEvery == 0:
            graphSin.SetPoint(n,x,y)
    else:
        target = -1
        if n % saveEvery == 0:
            graphBin.SetPoint(n,x,y)


    """
    ######################
    Neural net: scrambling
    ######################
    """
    indx = (n-scrStart)/scrLen-1
    if cmd.doScramble and n > scrStart and (n-scrStart) % scrLen == 0 and indx < cmd.Nperceptrons:
        indx = cmd.Nperceptrons - 1 - indx
        print "=== Scrambling perceptron [", indx, "] ==="
        NN.scramble({indx:[-1]})


    cost = NN.learn([x,y],[target])
    
    if n % saveEvery == 0:
        graphNNCost.SetPoint(n/saveEvery,n,cost)

        
        """
        #################################################
        Neural net: saving activation function derivative
        #################################################
        """
        k = 0
        for j,P in enumerate(NN.FFperceptrons):
            for i,N in enumerate(P.neurons):
                if n == 0:
                    graphDeriv.append(TGraph())
                    leg = "P:" + str(j) +  "; N:" + str(i)
                    legDeriv.AddEntry(graphDeriv[k],leg,"L");
                graphDeriv[k].SetPoint(n/saveEvery,n,N.dafundz)
                k += 1

        """
        #################################
        Neural net: saving delta function
        #################################
        """
        k = 0
        for j,P in enumerate(NN.BPperceptrons):
            for i,N in enumerate(P.neurons):
                if n == 0:
                    graphDelta.append(TGraph())
                    leg = "P:" + str(j) +  "; N:" + str(i)
                    legDelta.AddEntry(graphDelta[k],leg,"L");
                graphDelta[k].SetPoint(n/saveEvery,n,N.afun)
                k += 1

                
NN.printParams()
NN.save("NeuralNet.txt")


"""
################
Neural net: test
################
"""
print "\n=== Testing neural netowrk ==="
count = 0.
for n in xrange(nRuns):
    x = random() * xRange + xOffset - xRange/2
    y = random() * yRange + yOffset - yRange/2

    NNoutput = NN.eval([x,y])
    
    if ((NNoutput[0] > 0 and loR <= xyCorr(x,y) < hiR) or (NNoutput[0] <= 0 and (xyCorr(x,y) < loR or hiR <= xyCorr(x,y)))):
        count += 1

    if n % saveEvery == 0:
        if (NNoutput[0] > 0 and loR <= xyCorr(x,y) < hiR):
            graphSout.SetPoint(n/saveEvery,x,y)
        elif (NNoutput[0] <= 0 and (xyCorr(x,y) < loR or hiR <= xyCorr(x,y))):
            graphBout.SetPoint(n/saveEvery,x,y)
        else:
            graphNNerr.SetPoint(n/saveEvery,x,y)

    if NNoutput[0] > 0:
        histoNNS.Fill(NNoutput[0])
    elif NNoutput[0] <= 0:
        histoNNB.Fill(NNoutput[0])

print "\n=== Success rate:", count / nRuns * 100., "% ==="


"""
#########################
Neural net: control plots
#########################
"""
cCost.cd()
graphNNCost.Draw('AL')
cCost.Modified()
cCost.Update()

cNNin.cd()
cNNin.DrawFrame(xOffset - xRange/2 - xRange/20,yOffset - yRange/2 - xRange/20,xOffset + xRange/2 + xRange/20,yOffset + yRange/2 + xRange/20)
graphSin.Draw('P')
graphBin.Draw('P same')
cNNin.Modified()
cNNin.Update()

cNNout.cd()
cNNout.DrawFrame(xOffset - xRange/2 - xRange/20,yOffset - yRange/2 - xRange/20,xOffset + xRange/2 + xRange/20,yOffset + yRange/2 + xRange/20)
graphSout.Draw('P')
graphBout.Draw('P same')
graphNNerr.Draw('P same')
cNNout.Modified()
cNNout.Update()

cNNval.cd()
histoNNS.Draw()
histoNNB.Draw('same')
cNNval.Modified()
cNNval.Update()

cDeriv.cd()
graphDeriv[0].Draw('AL')
graphDeriv[0].SetTitle('NN activation function derivative;Epoch [#];Activation Function Derivative')
graphDeriv[0].SetLineColor(1)
legDeriv.Draw("same")
for k in xrange(1,len(graphDeriv[:])):
    graphDeriv[k].SetLineColor(k+1)
    graphDeriv[k].Draw('L same')
cDeriv.Modified()
cDeriv.Update()

cDelta.cd()
graphDelta[0].Draw('AL')
graphDelta[0].SetTitle('NN delta function;Epoch [#];Delta Function')
graphDelta[0].SetLineColor(1)
legDelta.Draw("same")
for k in xrange(1,len(graphDelta[:])):
    graphDelta[k].SetLineColor(k+1)
    graphDelta[k].Draw('L same')
cDelta.Modified()
cDelta.Update()


"""
########################
Wait for keyborad stroke
########################
"""
raw_input("\n---press enter---")
