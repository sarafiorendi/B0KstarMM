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

    def learn(self,invec,target):
        ###############
        # Feedforward #
        ###############
        result = self.eval(invec)
        error  = [ a - b for a,b in zip(target,result) ]

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

    def scramble(self,who):
        if who and who[0] == -1:
            who = [ [-1] for j in xrange(self.Nperceptrons) ]
            
        for j in who:
            self.FFperceptron[j].scramble(who[j])

    def copyFFtoBPweights(self):
        for j in xrange(self.Nperceptrons-1):
            for i in xrange(self.BPperceptron[j].Nneurons):
                self.BPperceptron[j].neuron[i].dafundz = self.FFperceptron[self.Nperceptrons-j-2].neuron[i].dafundz
                for k in xrange(self.BPperceptron[j].neuron[i].Nvars):
                    self.BPperceptron[j].neuron[i].weights[k] = self.FFperceptron[self.Nperceptrons-j-1].neuron[k].weights[i]
                    
    def save(self,name):
        f = open(name,"w")
        
        f.write("### M.E.D. Neural Network ###\n")
        f.write("# Nvars, Nperceptrons, Nneurons\n")
        out = str(self.FFperceptron[0].neuron[0].Nvars) + " " + str(self.Nperceptrons) + " ["
        for j in xrange(self.Nperceptrons):
            out += " " + str(self.FFperceptron[j].Nneurons)
        out += " ]\n\n"
        f.write(out)

        for j in xrange(self.Nperceptrons):
            f.write("Perceptron[ " + str(j) + " ]\n")
            self.FFperceptron[j].save(f)
            
        f.close()

    def read(self,name):
        f    = open(name,"r")
        line = f.readline()
        lele = line.split()
        
        while len(lele) == 0 or (len(lele) > 0 and "#" in lele[0]):
            line = f.readline()
            lele = line.split()
        
        self.__init__(int(lele[0]),int(lele[1]),[ int(lele[i]) for i in xrange(2,len(lele)) if lele[i].isdigit() ])

        for j in xrange(self.Nperceptrons):
            self.FFperceptron[j].read(f)

        f.close()
