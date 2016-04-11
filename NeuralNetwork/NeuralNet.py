####################################################
# Documentation: neuralnetworksanddeeplearning.com #
####################################################
# Cost function: 1/2 (target-result)^2             #
# Learning algorithm: online/incremental learning  #
#                     with gradient descent        #
####################################################
from random     import gauss

from Perceptron import Perceptron

class NeuralNet(object):

    ####################################################################
    # Nvars        = number of input variables                         #
    # Nperceptrons = number of perceptrons of neural net (must be > 1) #
    # Nneurons     = list of number of neurons for each perceptron     #
    ####################################################################
    # e.g.:                                                            #
    #       a network with 2 inputs, a first layer of 2 neurons, an    #
    #       intermediate layer of 2 neurons, and an output of 1 neuron #
    #       2, 3, [2,2,1]                                              #
    ####################################################################
    def __init__(self,Nvars,Nperceptrons,Nneurons):
        self.Nperceptrons = Nperceptrons
        self.FFperceptron = []
        self.BPperceptron = []

        #######################
        # Feedforward network #
        #######################
        ### First layer ###
        self.FFperceptron.append(Perceptron(Nneurons[0],Nvars,False))

        ### Intermediate layers + last layer ###
        self.FFperceptron.extend(Perceptron(Nneurons[j],self.FFperceptron[j-1].Nneurons,False) for j in xrange(1,self.Nperceptrons))

        ###########################
        # Backpropagation network #
        ###########################
        ### First layer ###
        self.BPperceptron.append(Perceptron(self.FFperceptron[self.Nperceptrons-2].Nneurons,self.FFperceptron[self.Nperceptrons-1].Nneurons,True))

        ### Intermediate layers + last layer ###
        self.BPperceptron.extend(Perceptron(self.FFperceptron[self.Nperceptrons-j-1].Nneurons,self.FFperceptron[self.Nperceptrons-j].Nneurons,True) for j in xrange(2,self.Nperceptrons))

    def eval(self,invec):
        ### First layer ###
        self.FFperceptron[0].eval(invec)

        ### Intermediate layers + last layer ###
        for j in xrange(1,self.Nperceptrons):
            self.FFperceptron[j].eval([ self.FFperceptron[j-1].neuron[i].afun for i in xrange(self.FFperceptron[j-1].Nneurons) ])

        return [ self.FFperceptron[self.Nperceptrons-1].neuron[i].afun for i in xrange(self.FFperceptron[self.Nperceptrons-1].Nneurons) ]

    def evalBP(self,invec):
        ### First layer ###
        self.BPperceptron[0].eval(invec)

        ### Intermediate layers + last layer ###
        for j in xrange(1,self.Nperceptrons-1):
            self.BPperceptron[j].eval([ self.BPperceptron[j-1].neuron[i].afun for i in xrange(self.BPperceptron[j-1].Nneurons) ])

    def learn(self,invec,target,stoGrad):
        ###############
        # Feedforward #
        ###############
        result = self.eval(invec)
        error  = [ a - b for a,b in zip(target,result) ]

        # @TMP@
        if stoGrad is True:
            rnd   = [ gauss(0,1./10) for i in xrange(self.FFperceptron[self.Nperceptrons-1].Nneurons) ]
            error = [ a + b for a,b in zip(error,rnd) ]

        ###################
        # Backpropagation #
        ###################
        self.copyFFtoBPweights()
        dcdz = [ error[i] * self.FFperceptron[self.Nperceptrons-1].neuron[i].dafundz for i in xrange(self.FFperceptron[self.Nperceptrons-1].Nneurons) ]
        self.evalBP(dcdz)

        ##########
        # Update #
        ##########
        ### Last layer ###
        self.FFperceptron[self.Nperceptrons-1].adapt([ self.FFperceptron[self.Nperceptrons-2].neuron[i].afun for i in xrange(self.FFperceptron[self.Nperceptrons-2].Nneurons) ],dcdz)

        ### Intermediate layers ###
        if self.Nperceptrons > 2:
            for j in xrange(self.Nperceptrons-2):
                self.FFperceptron[self.Nperceptrons-j-2].adapt([ self.FFperceptron[self.Nperceptrons-j-3].neuron[i].afun for i in xrange(self.FFperceptron[self.Nperceptrons-j-3].Nneurons) ],
                                                               [ self.BPperceptron[j].neuron[i].afun for i in xrange(self.BPperceptron[j].Nneurons) ])

        ### First layer ###
        self.FFperceptron[0].adapt(invec,[ self.BPperceptron[self.Nperceptrons-2].neuron[i].afun for i in xrange(self.BPperceptron[self.Nperceptrons-2].Nneurons) ])

        return sum(1./2 * a*a for a in error)

    def printParams(self):
        print "\n\n===>>> Feedforward Neural Net ===>>>"
        for j in xrange(self.Nperceptrons):
            print "Perceptron[", j, "]"
            self.FFperceptron[j].printParams()

        print "\n<<<=== Backpropagation Neural Net <<<==="
        for j in xrange(self.Nperceptrons-1):
            print "Perceptron[", j, "]"
            self.BPperceptron[j].printParams()

    def reset(self,what):
        for j in xrange(self.Nperceptrons):
            self.FFperceptron[j].reset(what)

        for j in xrange(self.Nperceptrons-1):
            self.BPperceptron[j].reset(what)

    def copyFFtoBPweights(self):
        for j in xrange(self.Nperceptrons-1):
            for i in xrange(self.BPperceptron[j].Nneurons):
                self.BPperceptron[j].neuron[i].dafundz = self.FFperceptron[self.Nperceptrons-j-2].neuron[i].dafundz
                for k in xrange(self.BPperceptron[j].neuron[i].Nvars):
                    self.BPperceptron[j].neuron[i].weights[k] = self.FFperceptron[self.Nperceptrons-j-1].neuron[k].weights[i]
                    
    def save(self):
        f = open("NeuralNet.txt","w")
        f.write("### M.E.D. Neural Network ###\n")
        f.write("# Nvars, Nperceptrons, NNeurons\n")
        f.write(str(self.FFperceptron[0].neuron[0].Nvars) + " " + str(self.Nperceptrons) + " " + str([ self.FFperceptron[j].Nneurons for j in xrange(self.Nperceptrons) ]) + "\n\n")
        for j in xrange(self.Nperceptrons):
            f.write("Perceptron[ " + str(j) + " ]\n")
            self.FFperceptron[j].save(f)
        f.close()

    def read(self):
        f = open("NeuralNet.txt","r")
        fin = f.readline()
        f.close()
