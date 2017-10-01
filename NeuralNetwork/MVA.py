"""""
#######################################################
MVA implementation with Neural Network
MVA.py                              by Mauro E. Dinardo
#######################################################
Before running check hyper-parameter space:
  .number of perceptrons & neurons
  .number of mini-batches
  .RMSprop and regularization
  .learn rate
  .scramble
  .dropout
  .cost function: quadratic / cross-entropy

To-do:
  - add ROC curve
  - test mini-batch
  - test reguarization
  - test RMSprop
  - activation function: ReLU and softmax&logLikelihood
#######################################################
e.g.: python MVA.py -nv 2 -np 5 -nn 2 3 3 2 1 -sc
      Neural Network with two inputs and one output
#######################################################
"""
from argparse  import ArgumentParser
from random    import seed, random, gauss
from math      import sqrt, exp
from os        import system
from timeit    import default_timer

from ROOT      import gROOT, gStyle, TCanvas, TGraph, TH1D, TGaxis, TLegend

from NeuralNet import NeuralNet


def ArgParser():
    """
    ###############
    Argument parser
    ###############
    """
    parser = ArgumentParser()
    parser.add_argument('-in', '--inFile',       dest = 'inFile',       type = str,          help = 'Input neural network',             required=False, default='')
    parser.add_argument('-nv', '--Nvars',        dest = 'Nvars',        type = int,          help = 'Number of variables',              required=False)
    parser.add_argument('-np', '--Nperceptrons', dest = 'Nperceptrons', type = int,          help = 'Number of perceptrons',            required=False)
    parser.add_argument('-nn', '--Nneurons',     dest = 'Nneurons',     type = int,          help = 'Number of neurons per perceptron', required=False, nargs='*')
    parser.add_argument('-sc', '--Scramble',     dest = 'scramble',     action='store_true', help = 'Do scramble',                      required=False)

    options = parser.parse_args()

    print ''
    if options.inFile:
        print '--> I\'m reading the input file:', options.inFile

    if options.Nvars:
        print '--> I\'m reading the variable number:', options.Nvars

    if options.Nperceptrons:
        print '--> I\'m reading the perceptron number:', options.Nperceptrons

    if options.Nneurons:
        print '--> I\'m reading the neuron number per perceptron:', options.Nneurons

    if options.scramble:
        print '--> I\'m reading the scramble flag:', options.scramble

    return options


def SetStyle():
    """
    ###################
    ROOT style settings
    ###################
    """
    gROOT.SetStyle('Plain')
    gROOT.ForceStyle()
    gStyle.SetTextFont(42)

    gStyle.SetOptTitle(0)
    gStyle.SetOptFit(0)
    gStyle.SetOptStat(1111)

    gStyle.SetPadRightMargin(0.08)
    gStyle.SetPadTopMargin(0.11)
    gStyle.SetPadBottomMargin(0.12)

    gStyle.SetTitleFont(42,'x')
    gStyle.SetTitleFont(42,'y')
    gStyle.SetTitleFont(42,'z')
    
    gStyle.SetTitleOffset(1.05,'x')
    gStyle.SetTitleOffset(0.95,'y')

    gStyle.SetTitleSize(0.05,'x')
    gStyle.SetTitleSize(0.05,'y')
    gStyle.SetTitleSize(0.05,'z')
    
    gStyle.SetLabelFont(42,'x')
    gStyle.SetLabelFont(42,'y')
    gStyle.SetLabelFont(42,'z')

    gStyle.SetLabelSize(0.05,'x')
    gStyle.SetLabelSize(0.05,'y')
    gStyle.SetLabelSize(0.05,'z')
    
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

if cmd.inFile:
    NN = NeuralNet()
    NN.read(cmd.inFile)
else:
    NN = NeuralNet(cmd.Nvars,cmd.Nperceptrons,cmd.Nneurons)
NN.printParams()


"""
###############
Hyperparameters
###############
"""
nRunTrainingLn = 100000
nRunTrainingSt = 0
nRunTest       = 10000
miniBatch      = 1

toScramble     = {3:[2]}

learnRateStart = 0.01
learnRateEnd   = 0.01
learnRateTau   = 1


"""
###########################################
Read additional hyper-parameter information
###########################################
"""
if cmd.inFile:
    nRunTrainingSt,miniBatch,learnRateStart,learnRateEnd,learnRateTau,toScramble = NN.readHypPar(cmd.inFile)


"""
#####################################
Internal parameters: problem specific
#####################################
"""
saveEvery = 100

NNoutMin  = NN.aFunMin(NN.Nperceptrons-1)
NNoutMax  = NN.aFunMax(NN.Nperceptrons-1)
NNthr     = (NNoutMin + NNoutMax) / 2.

xRange    = 2.
xOffset   = 0.2
yRange    = 2.
yOffset   = 0.2

noiseBand = 0.05
loR       = 0.2
hiR       = 0.5

xyDist    = lambda x,y: sqrt((x-xOffset)*(x-xOffset) + (y-yOffset)*(y-yOffset))


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
graphBin.SetMarkerColorAlpha(2,0.5)

graphSout = TGraph()
graphSout.SetTitle('NN output;x;y')
graphSout.SetMarkerStyle(20)
graphSout.SetMarkerSize(0.5)
graphSout.SetMarkerColor(4)

graphBout = TGraph()
graphBout.SetTitle('NN output;x;y')
graphBout.SetMarkerStyle(20)
graphBout.SetMarkerSize(0.5)
graphBout.SetMarkerColorAlpha(2,0.5)

graphNNerr = TGraph()
graphNNerr.SetTitle('NN output;x;y')
graphNNerr.SetMarkerStyle(20)
graphNNerr.SetMarkerSize(0.5)
graphNNerr.SetMarkerColorAlpha(3,0.25)

histoNNS = TH1D('histoNNS','histoNNS',100,NNoutMin,NNoutMax)
histoNNS.SetTitle('NN signal output;NN output;Entries [#]')
histoNNS.SetLineColor(4)

histoNNB = TH1D('histoNNB','histoNNB',100,NNoutMin,NNoutMax)
histoNNB.SetTitle('NN background output;NN output;Entries [#]')
histoNNB.SetLineColor(2)

histoNNE = TH1D('histoNNE','histoNNE',100,NNoutMin,NNoutMax)
histoNNE.SetTitle('NN error output;NN output;Entries [#]')
histoNNE.SetLineColor(3)

legNNspeed = TLegend(0.92, 0.17, 1.0, 1.0, '')
legNNspeed.SetTextSize(0.03)
legNNspeed.SetFillStyle(1001)

graphNNspeed = []


"""
####################
Neural net: training
####################
"""
print '\n\n=== Training neural network ==='
NNspeed = [ 0. for j in xrange(NN.Nperceptrons) ]
NNcost  = 0.
count   = 0.
startClock = default_timer()
for n in xrange(nRunTrainingSt+1,nRunTrainingSt+1 + nRunTrainingLn):
    """
    ####################p
    Neural net: training
    ####################
    """
    x = random() * xRange - xRange/2
    y = random() * yRange - yRange/2

    if gauss(loR,noiseBand) <= xyDist(x,y) < gauss(hiR,noiseBand):
        target = NNoutMax
        if n % saveEvery == 0:
            graphSin.SetPoint((n-nRunTrainingSt) / saveEvery - 1,x,y)
    else:
        target = NNoutMin
        if n % saveEvery == 0:
            graphBin.SetPoint((n-nRunTrainingSt) / saveEvery - 1,x,y)


    """
    ######################
    Neural net: scrambling
    ######################
    """
    if cmd.scramble == True and toScramble != None and n == nRunTrainingSt+1:
        NN.release({-1:[]})
        print '  [', n, '] Scrambling', toScramble
        NN.scramble(toScramble)
        NN.fixAllBut(toScramble)


    """
    ##############################
    Neural net: setting learn rate
    ##############################
    """
    NN.setLearnRate(learnRateEnd*(1 - exp(-(n-1.) / learnRateTau)) + learnRateStart*exp(-(n-1.) / learnRateTau))


    """
    ####################
    Neural net: learning
    ####################
    """
    NNcost += NN.learn([x,y],[target],miniBatch) if n % miniBatch == 0 else NN.learn([x,y],[target])


    """
    #####################################################
    Neural net: saving activation function speed and cost
    #####################################################
    """
    NNspeed = [ a + NN.speed(j) for j,a in enumerate(NNspeed) ]
    if n % saveEvery == 0:
        graphNNcost.SetPoint((n-nRunTrainingSt) / saveEvery - 1,n,NNcost / saveEvery)
        NNcost = 0.


        for j,a in enumerate(NNspeed):
            if (n-nRunTrainingSt) / saveEvery == 1:
                graphNNspeed.append(TGraph())
                leg = 'P:' + str(j)
                legNNspeed.AddEntry(graphNNspeed[j],leg,'L');
            graphNNspeed[j].SetPoint((n-nRunTrainingSt) / saveEvery - 1,n,a / saveEvery)
        NNspeed = [ 0. for j in xrange(NN.Nperceptrons) ]


    """
    ################
    Neural net: test
    ################
    """
    x = random() * xRange - xRange/2
    y = random() * yRange - yRange/2

    
    NNout = NN.eval([x,y])


    if ((NNout[0] > NNthr and loR <= xyDist(x,y) < hiR) or (NNout[0] <= NNthr and (xyDist(x,y) < loR or hiR <= xyDist(x,y)))):
        count += 1.

    if n % saveEvery == 0:
        graphNNaccuracy.SetPoint((n-nRunTrainingSt) / saveEvery - 1,n,count / saveEvery * 100)
        count = 0.

endClock = default_timer()
print '--> Training time:', endClock - startClock, '[s]'

NN.printParams()
NN.save('NeuralNet.txt')


"""
###########################################
Save additional hyper-parameter information
###########################################
"""
NN.saveHypPar('NeuralNet.txt',nRunTrainingSt+nRunTrainingLn,miniBatch,learnRateStart,learnRateEnd,learnRateTau,cmd.scramble,toScramble)


"""
################
Neural net: test
################
"""
print '\n\n=== Testing neural network ==='
startClock = default_timer()
for n in xrange(0,nRunTest):
    x = random() * xRange - xRange/2
    y = random() * yRange - yRange/2


    NNout = NN.eval([x,y])


    if NNout[0] > NNthr and loR <= xyDist(x,y) < hiR:
        graphSout.SetPoint(n,x,y)
        histoNNS.Fill(NNout[0])
    elif NNout[0] <= NNthr and (xyDist(x,y) < loR or hiR <= xyDist(x,y)):
        graphBout.SetPoint(n,x,y)
        histoNNB.Fill(NNout[0])
    else:
        graphNNerr.SetPoint(n,x,y)
        histoNNE.Fill(NNout[0])

endClock = default_timer()
print '--> Testing time:', endClock - startClock, '[s]'


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
if len(graphNNspeed[:]) > 0:
    graphNNspeed[0].Draw('AL')
    graphNNspeed[0].SetTitle('NN activation function speed;Epoch [#];Activation Function Speed')
    graphNNspeed[0].SetLineColor(1)
for k in xrange(1,len(graphNNspeed[:])):
    graphNNspeed[k].SetLineColor(k+1)
    graphNNspeed[k].Draw('L same')
legNNspeed.Draw('same')
cSpeed.Modified()
cSpeed.Update()

cNNin.cd()
cNNin.DrawFrame(-xRange/2 - xRange/20,-yRange/2 - xRange/20,+xRange/2 + xRange/20,+yRange/2 + xRange/20)
graphSin.Draw('P')
graphBin.Draw('P same')
cNNin.Modified()
cNNin.Update()

cNNout.cd()
cNNout.DrawFrame(-xRange/2 - xRange/20,-yRange/2 - xRange/20,+xRange/2 + xRange/20,+yRange/2 + xRange/20)
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
system('say \'Neural netowrk optimized\'')
raw_input('\n\n---press enter---')
