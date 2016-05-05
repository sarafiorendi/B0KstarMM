from random import gauss
from math   import sqrt, log, tanh, atanh
"""
#########################################
.Activation function: tanh, BPN
.Cost functions:
  1/2 (result - target)^2
  -log(sqrt(result))-target*atanh(result)
.Regularization: L2
#########################################
"""
class Neuron(object):
    learnRate  =  0.01
    rmsPrDecay =  0.99 # If = 1 then no RMSprop
    regular    =  0.   # If = 0 then no reguarization

    aFunMin    = -1.
    aFunMax    = +1.
    daFunDzMax =  1.
    """
    ######################################
    Nvars    = number of input variables
    aFunType = type of activation function
               "tanh", "BPN"
    ######################################
    """
    def __init__(self,Nvars,afunType="tanh"):
        self.Nvars = Nvars
        
        if afunType == "tanh" or afunType == "BPN":
            self.afunType = afunType
        else:
            print "[Neuron::__init__]\tWrong option:", afunType
            quit()

        self.weights = [ gauss(0,1 / sqrt(self.Nvars)) for k in xrange(self.Nvars) ]
        self.weights.append(gauss(0,1))
        if self.afunType is "BPN":
            self.weights = [ 0 for k in xrange(self.Nvars+1) ]

        self.afun    = 0
        self.dafundz = 0

        if self.rmsPrDecay == 1:
            self.rmsProp = 1
        else:
            self.rmsProp = 0

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
        self.dafundz = self.daFunDz()

        return [self.afun, self.dafundz]

    def adapt(self,invec,dCdZ):
        self.rmsProp = self.rmsPrDecay * self.rmsProp + (1 - self.rmsPrDecay) * dCdZ * dCdZ

        for k in xrange(self.Nvars):
            self.weights[k] = self.weights[k] * (1 - self.learnRate * self.regular) - self.learnRate * dCdZ * invec[k] / sqrt(self.rmsProp)
        self.weights[self.Nvars] -= self.learnRate * dCdZ

    ### Activation function ###
    def aFun(self,val):
        if self.afunType == "BPN":
            return val * self.dafundz

        if self.afunType == "tanh":
            return tanh(val)

    ### d(Activation function) / dz ###
    def daFunDz(self):
        if self.afunType == "BPN":
            return self.dafundz

        if self.afunType == "tanh":
            return 1 - self.afun * self.afun

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
        self.__init__(self.Nvars,self.afunType)

    def sum2W(self):
        return sum(W*W for W in self.weights[:-1])

    def scramble(self):
        for k in xrange(self.Nvars):
            self.weights[k] = self.weights[k] - cmp(self.weights[k],1) * gauss(0,(1 - self.afun) / self.dafundz / sqrt(self.Nvars))
        self.weights[self.Nvars] = self.weights[self.Nvars] - cmp(self.weights[self.Nvars],1) * gauss(0,(1 - self.afun) / self.dafundz)

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

        self.afunType = next(lele[i+2] for i,a in enumerate(lele) if a == "type")

        w.pop(0)
        self.afun     = w.pop(0)
        self.dafundz  = w.pop(0)
        self.weights  = w
