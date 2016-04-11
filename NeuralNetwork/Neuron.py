from random import gauss
from math   import tanh, sqrt, exp
#############################
# Activation function: tanh #
#############################

class Neuron(object):

    ############################################
    # Nvars     = number of input variables    #
    # isBackPrN = is a backpropagation network #
    ############################################
    # Return the value of the activation       #
    # function and its derivative              #
    ############################################
    def __init__(self,Nvars,isBackPrN):
        self.Nvars     = Nvars
        self.isBackPrN = isBackPrN        
        if self.isBackPrN is False:
            self.weights = [ gauss(0,1./sqrt(Nvars)) for k in xrange(self.Nvars+1) ]
        else:
            self.weights = [ 0 for k in xrange(self.Nvars) ]
        self.afun    = 0
        self.dafundz = 0
        self.epoch   = 0
        
    def eval(self,invec):
        wsum = 0
        for k in xrange(self.Nvars):
            wsum += self.weights[k] * invec[k]
        if self.isBackPrN is False:
            wsum += self.weights[self.Nvars]

        self.afun    = self.aFun(wsum)
        self.dafundz = self.daFunDz(wsum)

        return [self.afun, self.dafundz]

    def adapt(self,invec,dcdz):
        for k in xrange(self.Nvars):
                self.weights[k] = self.weights[k] + self.learnRate(self.epoch) * dcdz * invec[k]
        self.weights[self.Nvars] = self.weights[self.Nvars] + self.learnRate(self.epoch) * dcdz
        self.epoch += 1
        
    def aFun(self,val):
        if self.isBackPrN is False:
            return tanh(val)
        else:
            return val * self.dafundz

    def daFunDz(self,val):
        if self.isBackPrN is False:
            return 1 - self.afun * self.afun
        else:
            return self.dafundz

    def printParams(self):
        for k in xrange(len(self.weights)):
            print "    Weight[", k, "] ", round(self.weights[k],2)

    def reset(self,what):
        if what == "all":
            self.__init__(self.Nvars,self.isBackPrN)
        elif what == "epoch":
            self.epoch = 0
        else:
            print "[Neuron::reset]\toption not valid:", what
            quit()

    def save(self,f):
        for k in xrange(len(self.weights)):
            f.write("    Weight[ " + str(k) + " ] " + str(self.weights[k]) + "\n")

    def read(self,f):
        f.readline()

    def learnRate(self,epoch):
        ltstart = 100
        ltend   = 100
        return 1. / (ltstart + ltend * (1 - exp(-epoch / ltend)))
