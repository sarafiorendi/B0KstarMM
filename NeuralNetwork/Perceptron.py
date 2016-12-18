from math   import sqrt
from Neuron import Neuron

class Perceptron(object):
    """
    ####################################################
    Nneurons = number of neurons of the perceptron
    Nvars    = number of input variables for each neuron
    afunType = type of activation function "tanh", "BPN"
    ####################################################
    """
    def __init__(self,Nneurons,Nvars,afunType="tanh"):
        self.Nneurons = Nneurons
        self.neurons  = [ Neuron(Nvars,afunType) for i in xrange(self.Nneurons) ]

    def eval(self,invec):
        return [ N.eval(invec) for N in self.neurons ]

    def adapt(self,invec,dCdZ,miniBatch=0):
        for a,N in zip(dCdZ,self.neurons):
            N.adapt(invec,a,miniBatch)

    def cFun(self,target):
        return sum(N.cFun(a) for a,N in zip(target,self.neurons))

    def dcFunDz(self,target):
        return [ N.dcFunDz(a) for a,N in zip(target,self.neurons) ]

    def speed(self):
        return sqrt(sum(N.afun * N.afun for N in self.neurons))

    def printParams(self):
        for i,N in enumerate(self.neurons):
            print "  Neuron[", i, "type =", N.afunType, "aFun =", round(N.afun,2), "d(aFun)/dz =", round(N.dafundz,2), "]"
            N.printParams()

    def reset(self):
        for N in self.neurons:
            N.reset()

    def sum2W(self):
        return sum(N.sum2W() for N in self.neurons)
            
    def scramble(self,who):
        if who[0] == -1:
            who = [ i for i in xrange(self.Nneurons) ]

        for i in who:
            self.neurons[i].scramble()

    def fixAllBut(self,who):
        genExp = (N for (i,N) in enumerate(self.neurons) if i not in who)
        for N in genExp:
            N.amIfixed = True

    def release(self,who):
        if who[0] == -1:
            who = [ i for i in xrange(self.Nneurons) ]

        genExp = (N for (i,N) in enumerate(self.neurons) if i in who)
        for N in genExp:
            N.amIfixed = False

    def aFunMin(self,indx):
        return self.neurons[indx].aFunMin

    def aFunMax(self,indx):
        return self.neurons[indx].aFunMax

    def removeN(self,who):
        if who[0] == -1:
            who = [ i for i in xrange(self.Nneurons) ]

        self.neurons  = [ N for i,N in enumerate(self.neurons) if i not in who ]
        self.Nneurons = len(self.neurons[:])

    def removeW(self,who):
        for N in self.neurons:
            N.removeW(who) if type(who) is list else N.__init__(who,N.afunType)

    def addN(self,who):
        self.Nneurons += len(who[:])
        
        for i,pos in enumerate(who):
            self.neurons.insert(pos+i,Neuron(self.neurons[0].Nvars,self.neurons[0].afunType))

    def addW(self,who):
        for N in self.neurons:
            N.addW(who) if type(who) is list else N.__init__(who,N.afunType)

    def setLearnRate(self,val):
        for N in self.neurons:
            N.learnRate = val

    def save(self,f):
        for i,N in enumerate(self.neurons):
            f.write("  Neuron[ {0:d} ] type = {1:10s} aFun = {2:20f} d(aFun)/dz = {3:20f} Weights:".format(i,N.afunType,N.afun,N.dafundz))
            N.save(f)

    def read(self,f):
        line = f.readline()
        lele = line.split()
        
        while len(lele) == 0 or (len(lele) > 0 and ("#" in lele[0] or "Perceptron[" not in line)):
            line = f.readline()
            lele = line.split()

        for N in self.neurons:
            N.read(f)
