"""
################################################
Feedforward Neural Network
Cost function: 1/2 (target - result)^2
Learning algorithm: online/incremental learning
                    with gradient descent
################################################
Documentation: neuralnetworksanddeeplearning.com
################################################
"""
from Perceptron import Perceptron

class NeuralNet(object):
    """
    ################################################################
    Nvars        = number of input variables
    Nperceptrons = number of perceptrons of neural net (must be > 1)
    Nneurons     = list of number of neurons for each perceptron    
    ################################################################
    e.g.:
          a network with 2 inputs, a first layer of 2 neurons, an
          intermediate layer of 2 neurons, and an output of 1 neuron
          2, 3, [2,2,1]
    ################################################################
    """
    def __init__(self,Nvars,Nperceptrons,Nneurons):
        self.Nperceptrons  = Nperceptrons
        self.FFperceptrons = []

        """
        ###################
        Feedforward network
        ###################
        """
        ### First layer ###
        self.FFperceptrons.append(Perceptron(Nneurons[0],Nvars,False))

        ### Intermediate layers + last layer ###
        self.FFperceptrons.extend(Perceptron(Nneurons[j],self.FFperceptrons[j-1].Nneurons,False) for j in xrange(1,self.Nperceptrons))

        """
        #######################
        Backpropagation network
        #######################
        """
        self.initBPnetwork()

    def initBPnetwork(self):
        self.BPperceptrons = []

        ### First layer ###
        self.BPperceptrons.append(Perceptron(self.FFperceptrons[self.Nperceptrons-2].Nneurons,self.FFperceptrons[self.Nperceptrons-1].Nneurons,True))

        ### Intermediate layers + last layer ###
        self.BPperceptrons.extend(Perceptron(self.FFperceptrons[self.Nperceptrons-j-1].Nneurons,self.FFperceptrons[self.Nperceptrons-j].Nneurons,True) for j in xrange(2,self.Nperceptrons))
        
    def eval(self,invec):
        ### First layer ###
        self.FFperceptrons[0].eval(invec)

        ### Intermediate layers + last layer ###
        for j in xrange(1,self.Nperceptrons):
            self.FFperceptrons[j].eval([ N.afun for N in self.FFperceptrons[j-1].neurons ])

        return [ N.afun for N in self.FFperceptrons[self.Nperceptrons-1].neurons ]

    def evalBP(self,invec):
        ### First layer ###
        self.BPperceptrons[0].eval(invec)

        ### Intermediate layers + last layer ###
        for j in xrange(1,self.Nperceptrons-1):
            self.BPperceptrons[j].eval([ N.afun for N in self.BPperceptrons[j-1].neurons ])

    def learn(self,invec,target):
        """
        ###########
        Feedforward
        ###########
        """
        result = self.eval(invec)
        error  = [ a - b for a,b in zip(target,result) ]

        """
        ###############
        Backpropagation
        ###############
        """
        self.copyFFtoBPweights()
        dcdz = [ error[i] * N.dafundz for i,N in enumerate(self.FFperceptrons[self.Nperceptrons-1].neurons) ]
        self.evalBP(dcdz)

        """
        ######
        Update
        ######
        """
        ### Last layer ###
        self.FFperceptrons[self.Nperceptrons-1].adapt([ N.afun for N in self.FFperceptrons[self.Nperceptrons-2].neurons ],dcdz)

        ### Intermediate layers ###
        if self.Nperceptrons > 2:
            for j in xrange(self.Nperceptrons-2):
                self.FFperceptrons[self.Nperceptrons-j-2].adapt([ N.afun for N in self.FFperceptrons[self.Nperceptrons-j-3].neurons ],
                                                                [ N.afun for N in self.BPperceptrons[j].neurons ])

        ### First layer ###
        self.FFperceptrons[0].adapt(invec,[ N.afun for N in self.BPperceptrons[self.Nperceptrons-2].neurons ])

        return sum(1./2 * a*a for a in error)

    def printParams(self):
        print "\n\n===>>> Feedforward Neural Net ===>>>"
        for j,P in enumerate(self.FFperceptrons):
            print "Perceptron[", j, "]"
            P.printParams()

        print "\n<<<=== Backpropagation Neural Net <<<==="
        for j,P in enumerate(self.BPperceptrons):
            print "Perceptron[", j, "]"
            P.printParams()

    def reset(self,what):
        for P in self.FFperceptrons:
            P.reset(what)

        for P in self.BPperceptrons:
            P.reset(what)

    def scramble(self,who):
        if -1 in who.keys():
            who = { j:[-1] for j in xrange(self.Nperceptrons) }

        for j in who:
            self.FFperceptrons[j].scramble(who[j])

    def copyFFtoBPweights(self):
        for j in xrange(self.Nperceptrons-1):
            for i in xrange(self.BPperceptrons[j].Nneurons):
                self.BPperceptrons[j].neurons[i].dafundz = self.FFperceptrons[self.Nperceptrons-j-2].neurons[i].dafundz
                for k in xrange(self.BPperceptrons[j].neurons[i].Nvars):
                    self.BPperceptrons[j].neurons[i].weights[k] = self.FFperceptrons[self.Nperceptrons-j-1].neurons[k].weights[i]
                    
    def copy(self):
        out = NeuralNet(self.FFperceptrons[0].neurons[0].Nvars,self.Nperceptrons,[ P.Nneurons for P in self.FFperceptrons ])
        
        for j in xrange(self.Nperceptrons):
            for i in xrange(self.FFperceptrons[j].Nneurons):
                for k in xrange(len(self.FFperceptrons[j].neurons[i].weights[:])):
                    out.FFperceptrons[j].neurons[i].weights[k] = self.FFperceptrons[j].neurons[i].weights[k]

        return out

    def remove(self,who):
        if -1 in who.keys():
            who = { j:[-1] for j in xrange(self.Nperceptrons) }

        for j in who:
            if j < self.Nperceptrons-1:
                if who[j][0] == -1:
                    self.FFperceptrons[j+1].removeW(self.FFperceptrons[j].neurons[0].Nvars)
                    self.Nperceptrons -= 1
                else:
                    self.FFperceptrons[j+1].removeW(who[j])
            elif who[j][0] == -1:
                    self.Nperceptrons -= 1

            self.FFperceptrons[j].removeN(who[j])

        self.FFperceptrons = [ P for j,P in enumerate(self.FFperceptrons) if (j not in who.keys()) or (j in who.keys() and who[j][0] != -1) ]
        
        self.initBPnetwork()

    def add(self,who):
        print "[NeuralNet::add]\tNot implemented yet"
        quit()

    def save(self,name):
        f = open(name,"w")
        
        f.write("### M.E.D. Neural Network ###\n")
        f.write("# Nvars, Nperceptrons, Nneurons\n")
        out = str(self.FFperceptrons[0].neurons[0].Nvars) + " " + str(self.Nperceptrons) + " ["
        for P in self.FFperceptrons:
            out += " " + str(P.Nneurons)
        out += " ]\n\n"
        f.write(out)

        for j,P in enumerate(self.FFperceptrons):
            f.write("Perceptron[ " + str(j) + " ]\n")
            P.save(f)
            
        f.close()

    def read(self,name):
        f    = open(name,"r")
        line = f.readline()
        lele = line.split()
        
        while len(lele) == 0 or (len(lele) > 0 and "#" in lele[0]):
            line = f.readline()
            lele = line.split()
        
        self.__init__(int(lele[0]),int(lele[1]),[ int(lele[i]) for i in xrange(2,len(lele)) if lele[i].isdigit() ])

        for P in self.FFperceptrons:
            P.read(f)

        f.close()
