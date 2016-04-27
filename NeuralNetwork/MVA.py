"""
###############################################
MVA implementation with Neural Networks
                            by Mauro E. Dinardo
###############################################
Before running check on hyper-parameter space:
  - number of perceptrons & neurons
  - learn rate
  - regulator
  - scramble
  - Cost functions:
    quadratic / cross-entropy / softmax&logLikelihood

To-do:
  - mini-batch learning
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
    parser.add_argument("-nv", "--Nvars",        dest = "Nvars",        type = int, help = "Number of variables",              required=True)
    parser.add_argument("-np", "--Nperceptrons", dest = "Nperceptrons", type = int, help = "Number of perceptrons",            required=True)
    parser.add_argument("-nn", "--Nneurons",     dest = "Nneurons",     type = int, help = "Number of neurons per perceptron", required=True,  nargs='*')
    parser.add_argument("-sc", "--Nscramble",    dest = "Nscramble",    type = int, help = "Perceptrons to scramble",          required=False, nargs='*', default=[])

    options = parser.parse_args()

    print ""
    if options.Nvars:
        print "--> I'm reading the variable number: ", options.Nvars

    if options.Nperceptrons:
        print "--> I'm reading the perceptron number: ", options.Nperceptrons

    if options.Nneurons:
        print "--> I'm reading the neuron number per perceptron: ", options.Nneurons

    if options.Nscramble:
        print "--> I'm reading the perceptron number to scramble: ", options.Nscramble

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
    gStyle.SetOptStat(1111)
    
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
##########################
Neural net: initialization
##########################
"""
seed(0)
NN = NeuralNet(cmd.Nvars,cmd.Nperceptrons,cmd.Nneurons)
NN.read("NeuralNet_dropout.txt")
NN.printParams()
NN.add({3:[1,3]})


"""
###################
Internal parameters
###################
"""
seed(0)
nRuns     = 50000
saveEvery =   100
scrStart  =     0
scrLen    = 10000
xRange    = 3.
xOffset   = 3.
yRange    = 3.
yOffset   = 0.
noiseBand = 0.1
loR       = 0.5
hiR       = 1.
NNoutMin  = NN.outputMin(NN.Nperceptrons-1)
NNoutMax  = NN.outputMax(NN.Nperceptrons-1)
NNthr     = (NNoutMin+NNoutMax) / 2.
xyCorr    = lambda x,y: ((x-xOffset)*(x-xOffset)+(y-yOffset)*(y-yOffset))
#xyCorr    = lambda x,y: (6*(x-xOffset)*(x-xOffset)*(x-xOffset) - (y-yOffset)) # @TMP@


"""
#########################
Graphics layout and plots
#########################
"""
gROOT.Reset()
SetStyle()

cCost  = TCanvas('cCost',  'NN Cost Function', 0, 0, 700, 500)
cAccu  = TCanvas('cAccu',  'NN Accuracy',      0, 0, 700, 500)
cSpeed = TCanvas('cSpeed', 'NN Speed',         0, 0, 700, 500)
cNNin  = TCanvas('cNNin',  'NN Input',         0, 0, 700, 500)
cNNout = TCanvas('cNNout', 'NN Output',        0, 0, 700, 500)
cNNval = TCanvas('cNNval', 'NN Values',        0, 0, 700, 500)

graphNNcost = TGraph()
graphNNcost.SetTitle('NN cost function;Epoch [#];Cost Function')

graphNNaccuracy = TGraph()
graphNNaccuracy.SetTitle('NN accuracy;Epoch [#];Accuracy [%]')

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

histoNNS = TH1D('histoNNS','histoNNS',100,NNoutMin,NNoutMax)
histoNNS.SetTitle('NN signal output;NN output;Entries [#]')
histoNNS.SetLineColor(4)

histoNNB = TH1D('histoNNB','histoNNB',100,NNoutMin,NNoutMax)
histoNNB.SetTitle('NN background output;NN output;Entries [#]')
histoNNB.SetLineColor(2)

histoNNE = TH1D('histoNNE','histoNNE',100,NNoutMin,NNoutMax)
histoNNE.SetTitle('NN error output;NN output;Entries [#]')
histoNNE.SetLineColor(3)

legNNspeed = TLegend(0.88, 0.17, 1.0, 1.0, "")
legNNspeed.SetTextSize(0.03)
legNNspeed.SetFillStyle(1001)

graphNNspeed = []


"""
####################
Neural net: training
####################
"""
print "\n\n=== Training neural network ==="
NNspeed = [ 0. for j in xrange(NN.Nperceptrons) ]
NNcost  = 0.
count   = 0.
for n in xrange(nRuns):
    """
    ####################
    Neural net: training
    ####################
    """
    x = random() * xRange + xOffset - xRange/2
    y = random() * yRange + yOffset - yRange/2

    if gauss(loR,noiseBand) <= xyCorr(x,y) < gauss(hiR,noiseBand):
#    if xyCorr(x,y) > hiR: # @TMP@
        target = NNoutMax
        if n % saveEvery == 0:
            graphSin.SetPoint(n/saveEvery,x,y)
    else:
        target = NNoutMin
        if n % saveEvery == 0:
            graphBin.SetPoint(n/saveEvery,x,y)


    """
    ######################
    Neural net: scrambling
    ######################
    """
    if cmd.Nscramble and n == scrStart:
        print ""
    indx = NN.Nperceptrons - 1 - (n-scrStart)/scrLen
    if indx in cmd.Nscramble and n >= scrStart and (n-scrStart) % scrLen == 0:
        print "  Scrambling perceptron [", indx, "]"
        NN.scramble({indx:[-1]})

        
    NNcost += NN.learn([x,y],[target])
    NNspeed = [ a + NN.speed(j) for j,a in enumerate(NNspeed) ]

    
    if n % saveEvery == 0:
        graphNNcost.SetPoint(n/saveEvery,n,NNcost / saveEvery)
        NNcost = 0.

        
        """
        ############################################
        Neural net: saving activation function speed
        ############################################
        """
        k = 0
        for j,a in enumerate(NNspeed):
            if n == 0:
                graphNNspeed.append(TGraph())
                leg = "P:" + str(j)
                legNNspeed.AddEntry(graphNNspeed[k],leg,"L");
            graphNNspeed[k].SetPoint(n/saveEvery,n,a / saveEvery)
            k += 1
        NNspeed = [ 0. for j in xrange(NN.Nperceptrons) ]


    """
    ################
    Neural net: test
    ################
    """
    x = random() * xRange + xOffset - xRange/2
    y = random() * yRange + yOffset - yRange/2

    
    NNout = NN.eval([x,y])
     

    if ((NNout[0] > NNthr and loR <= xyCorr(x,y) < hiR) or (NNout[0] <= NNthr and (xyCorr(x,y) < loR or hiR <= xyCorr(x,y)))):
#    if ((NNout[0] > NNthr and xyCorr(x,y) > hiR) or (NNout[0] <= NNthr and xyCorr(x,y) <= hiR)): # @TMP@
        count += 1

    if n % saveEvery == 0:
        graphNNaccuracy.SetPoint(n/saveEvery,n,count / saveEvery * 100)
        count = 0.
NN.printParams()
NN.save("NeuralNet.txt")


"""
################
Neural net: test
################
"""
print "\n\n=== Testing neural network ==="
for n in xrange(nRuns):
    x = random() * xRange + xOffset - xRange/2
    y = random() * yRange + yOffset - yRange/2

    
    NNout = NN.eval([x,y])
     

    if (NNout[0] > NNthr and loR <= xyCorr(x,y) < hiR):
#    if (NNout[0] > NNthr and xyCorr(x,y) > hiR): # @TMP@
        graphSout.SetPoint(n,x,y)
        histoNNS.Fill(NNout[0])
    elif (NNout[0] <= NNthr and (xyCorr(x,y) < loR or hiR <= xyCorr(x,y))):
#    elif (NNout[0] <= NNthr and xyCorr(x,y) <= hiR): # @TMP@
        graphBout.SetPoint(n,x,y)
        histoNNB.Fill(NNout[0])
    else:
        graphNNerr.SetPoint(n,x,y)
        histoNNE.Fill(NNout[0])


"""
#########################
Neural net: control plots
#########################
"""
cCost.cd()
graphNNcost.Draw('AL')
cCost.Modified()
cCost.Update()

cAccu.cd()
cAccu.SetGrid()
graphNNaccuracy.Draw('AL')
cAccu.Modified()
cAccu.Update()

cSpeed.cd()
graphNNspeed[0].Draw('AL')
graphNNspeed[0].SetTitle('NN activation function speed;Epoch [#];Activation Function Speed')
graphNNspeed[0].SetLineColor(1)
for k in xrange(1,len(graphNNspeed[:])):
    graphNNspeed[k].SetLineColor(k+1)
    graphNNspeed[k].Draw('L same')
legNNspeed.Draw("same")
cSpeed.Modified()
cSpeed.Update()

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
histoNNB.Draw('sames')
histoNNE.Draw('sames')
cNNval.Modified()
cNNval.Update()


"""
########################
Wait for keyborad stroke
########################
"""
raw_input("\n\n---press enter---")
