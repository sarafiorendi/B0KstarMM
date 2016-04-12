from Neuron import Neuron

class Perceptron(object):

    #########################################################
    # Nneurons  = number of neurons of the perceptron       #
    # Nvars     = number of input variables for each neuron #
    # isBackPrN = is a backpropagation network              #
    #########################################################
    def __init__(self,Nneurons,Nvars,isBackPrN):
        self.Nneurons = Nneurons
        self.neuron   = [ Neuron(Nvars,isBackPrN) for i in xrange(self.Nneurons) ]

    def eval(self,invec):
        return [ self.neuron[i].eval(invec) for i in xrange(self.Nneurons) ]

    def adapt(self,invec,dcdz):
        for i in xrange(self.Nneurons):
            self.neuron[i].adapt(invec,dcdz[i])

    def printParams(self):
        for i in xrange(self.Nneurons):
            print "  Neuron[", i, "aFun =", round(self.neuron[i].afun,2), "d(aFun)/dz =", round(self.neuron[i].dafundz,2), "]"
            self.neuron[i].printParams()

    def reset(self,what):
        for i in xrange(self.Nneurons):
            self.neuron[i].reset(what)

    def save(self,f):
        for i in xrange(self.Nneurons):
            f.write("  Neuron[ {0:d} ] aFun = {1:20f} d(aFun)/dz = {2:20f}".format(i,self.neuron[i].afun,self.neuron[i].dafundz))
            self.neuron[i].save(f)

    def read(self,f):
        fin = f.readline()
