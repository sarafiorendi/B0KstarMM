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

    def printWeights(self):
        for i in xrange(self.Nneurons):
            print "------ Neuron[", i, "aFun =", round(self.neuron[i].afun,2), "d(aFun)/dz =", round(self.neuron[i].dafundz,2), "] ------"
            self.neuron[i].printWeights()

    def reset(self):
        for i in xrange(self.Nneurons):
            self.neuron[i].reset()
