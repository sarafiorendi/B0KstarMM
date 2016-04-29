from random import gauss
from math   import sqrt, log, tanh, atanh
"""
#########################################
.Activation function: tanh
.Cost functions:
  1/2 (result - target)^2
  -log(sqrt(result))-target*atanh(result)
.Regularization: L2
#########################################
"""
class Neuron(object):
    lRate      =  0.001
    regular    =  0.
    aFunMin    = -1.
    aFunMax    = +1.
    daFunDzMax =  1.
    """
    ####################################
    Nvars = number of input variables
    isBPN = is a backpropagation network
    ####################################
    """
    def __init__(self,Nvars,isBPN=False):
        self.Nvars = Nvars
        self.isBPN = isBPN

        self.weights = [ gauss(0,1 / sqrt(self.Nvars)) for k in xrange(self.Nvars) ]
        self.weights.append(gauss(0,1))
        if self.isBPN is True:
            self.weights = [ 0 for k in xrange(self.Nvars+1) ]

        self.afun    = 0
        self.dafundz = 0

    def eval(self,invec):
        """
        ##################################
        Return the value of the activation
        function and its derivative
        ##################################
        """
        wsum  = sum(W * i for W,i in zip(self.weights,invec))
        wsum += self.weights[self.Nvars]

        self.afun    = self.aFun(wsum)
        self.dafundz = self.daFunDz(wsum)

        return [self.afun, self.dafundz]

    def adapt(self,invec,dCdZ):
        for k in xrange(self.Nvars):
            self.weights[k] = self.weights[k] * (1 - self.learnRate() * self.regular) - self.learnRate() * dCdZ * invec[k]
        self.weights[self.Nvars] -= self.learnRate() * dCdZ

    ### Activation function ###
    def aFun(self,val):
        if self.isBPN is False:
            return tanh(val)
        else:
            return val * self.dafundz

    ### d(Activation function) / dz ###
    def daFunDz(self,val):
        if self.isBPN is False:
            return 1 - self.afun * self.afun
        else:
            return self.dafundz

    ### Cost function ###
    def cFun(self,target):
        # @TMP@
        return 1./2 * (self.afun - target) * (self.afun - target)
#        return -log(sqrt(self.dafundz)) - target * atanh(self.afun)

    ### d(Cost function) / dz ###
    def dcFunDz(self,target):
        # @TMP@
        return (self.afun - target) * self.dafundz
#        return self.afun - target

    def printParams(self):
        for k,W in enumerate(self.weights):
            print "    Weight[", k, "] ", round(W,2)

    def reset(self):
        self.__init__(self.Nvars,self.isBPN)

    def sum2W(self):
        return sum(W*W for W in self.weights[:-1])

    def scramble(self):
        for i in xrange(self.Nvars):
            self.weights[i] = self.weights[i] - cmp(self.weights[i],1) * gauss(0,self.afun / self.dafundz / sqrt(self.Nvars))
        self.weights[self.Nvars] = self.weights[self.Nvars] - cmp(self.weights[self.Nvars],1) * gauss(0,self.afun / self.dafundz)

    def removeW(self,who):
        self.weights = [ W for k,W in enumerate(self.weights) if k not in who ]
        self.Nvars   = len(self.weights[:]) - 1

    def addW(self,who):
        self.Nvars += len(who[:])        
        for k,pos in enumerate(who):
            self.weights.insert(pos+k,gauss(0,1 / sqrt(self.Nvars)))
            
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

        w = [ float(a) for a in lele if a.replace(".","").replace("-","").isdigit() ]

        w.pop(0)
        self.afun    = w.pop(0)
        self.dafundz = w.pop(0)
        self.weights = w

    def learnRate(self):
        return self.lRate
