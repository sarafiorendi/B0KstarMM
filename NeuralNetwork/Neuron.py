from random import gauss
from math   import tanh, sqrt, exp
"""
#########################
Activation function: tanh
#########################
"""
class Neuron(object):
    lrStart = 0.01
    lrEnd   = 0.0001
    tau     = 10000

    """
    ####################################
    Nvars = number of input variables
    isBPN = is a backpropagation network
    ####################################
    """
    def __init__(self,Nvars,isBPN=False):
        self.Nvars = Nvars
        self.isBPN = isBPN
        
        self.weights = [ gauss(0,1/sqrt(self.Nvars+1)) for k in xrange(self.Nvars+1) ]
        if self.isBPN is True:
            self.weights = [ 0 for k in xrange(self.Nvars+1) ]

        self.afun    = 0
        self.dafundz = 0
        self.epoch   = 0

    def eval(self,invec):
        """
        ##################################
        Return the value of the activation
        function and its derivative
        ##################################
        """
        wsum = 0
        for k in xrange(self.Nvars):
            wsum += self.weights[k] * invec[k]
        wsum += self.weights[self.Nvars]

        self.afun    = self.aFun(wsum)
        self.dafundz = self.daFunDz(wsum)

        return [self.afun, self.dafundz]

    def adapt(self,invec,dcdz):
        for k in xrange(self.Nvars):
                self.weights[k] = self.weights[k] + self.learnRate(self.epoch) * dcdz * invec[k]
        self.weights[self.Nvars] = self.weights[self.Nvars] + self.learnRate(self.epoch) * dcdz
        
        self.epoch += 1
        
    ### Activation function ###
    def aFun(self,val):
        if self.isBPN is False:
            return tanh(val)
        else:
            return val * self.dafundz

    ### d(Activation function)/dz ###
    def daFunDz(self,val):
        if self.isBPN is False:
            return 1 - self.afun * self.afun
        else:
            return self.dafundz

    def printParams(self):
        for k,W in enumerate(self.weights):
            print "    Weight[", k, "] ", round(W,2)

    def reset(self,what):
        """
        ##############
        what = "all"
        what = "epoch"
        ##############
        """
        if what == "all":
            self.__init__(self.Nvars,self.isBPN)
        elif what == "epoch":
            self.epoch = 0
        else:
            print "[Neuron::reset]\toption not valid:", what
            quit()

    def scramble(self):
        self.weights = [ gauss(W,(1 - self.dafundz) / sqrt(self.Nvars+1)) for W in self.weights ]

    def removeW(self,who):
            self.weights = [ W for k,W in enumerate(self.weights) if k not in who ]
            self.Nvars   = len(self.weights[:])

    def addW(self,who):
        self.Nvars += len(who[:])        
        for k,pos in enumerate(who):
            self.weights.insert(pos+k,gauss(0,1/sqrt(self.Nvars+1)))
            
    def save(self,f):
        out = ""
        for W in self.weights:
            out += "{0:20f}".format(W)
        out += "\n"
        f.write(out)

    def read(self,f):
        line = f.readline()
        lele = line.split()

        while len(lele) == 0 or (len(lele) > 0 and ("#" in lele[0] or "Neuron[" not in line)):
            line = f.readline()
            lele = line.split()

        w = [ float(lele[i]) for i in xrange(len(lele)) if lele[i].replace(".","").replace("-","").isdigit() ]

        w.pop(0)
        self.afun    = w.pop(0)
        self.dafundz = w.pop(0)
        self.weights = w

    def learnRate(self,epoch):
        return self.lrStart * exp(-epoch / self.tau) + self.lrEnd * (1 - exp(-epoch / self.tau))
