from Neuron import Neuron

class Perceptron(object):
    """
    ####################################################
    Nneurons = number of neurons of the perceptron
    Nvars    = number of input variables for each neuron
    isBPN    = is a backpropagation network
    ####################################################
    """
    def __init__(self,Nneurons,Nvars,isBPN=False):
        self.Nneurons = Nneurons
        self.neurons  = [ Neuron(Nvars,isBPN) for i in xrange(self.Nneurons) ]

    def eval(self,invec):
        return [ N.eval(invec) for N in self.neurons ]

    def adapt(self,invec,dcdz):
        for i,N in enumerate(self.neurons):
            N.adapt(invec,dcdz[i])

    def printParams(self):
        for i,N in enumerate(self.neurons):
            print "  Neuron[", i, "aFun =", round(N.afun,2), "d(aFun)/dz =", round(N.dafundz,2), "]"
            N.printParams()

    def reset(self,what):
        for N in self.neurons:
            N.reset(what)

    def scramble(self,who):
        if who[0] == -1:
            who = [ i for i in xrange(self.Nneurons) ]

        for i in who:
            self.neurons[i].scramble()

    def removeN(self,who):
        if who[0] == -1:
            who = [ i for i in xrange(self.Nneurons) ]

        self.neurons  = [ N for i,N in enumerate(self.neurons) if i not in who ]
        self.Nneurons = len(self.neurons[:])

    def removeW(self,who):
        for N in self.neurons:
            if type(who) is list:
                N.removeW(who)
            else:
                N.__init__(who)

    def addN(self,who):
        self.Nneurons += len(who[:])
        
        for i,pos in enumerate(who):
            self.neurons.insert(pos+i,Neuron(self.neurons[0].Nvars))

    def addW(self,who):
        for N in self.neurons:
            if type(who) is list:
                N.addW(who)
            else:
                N.__init__(who)

    def save(self,f):
        for i,N in enumerate(self.neurons):
            f.write("  Neuron[ {0:d} ] aFun = {1:20f} d(aFun)/dz = {2:20f}".format(i,N.afun,N.dafundz))
            N.save(f)

    def read(self,f):
        line = f.readline()
        lele = line.split()
        
        while len(lele) == 0 or (len(lele) > 0 and ("#" in lele[0] or "Perceptron[" not in line)):
            line = f.readline()
            lele = line.split()

        for N in self.neurons:
            N.read(f)
