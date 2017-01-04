"""
###########################################
.Neural Network type: feedforward
.Learning algorithm:  online/incremental
                      with gradient descent
                        by Mauro E. Dinardo
###########################################
"""
import os
import ast

from math       import sqrt
from Perceptron import Perceptron

class NeuralNet(object):
    """
    ################################################################
    Nvars        = number of input variables
    Nperceptrons = number of perceptrons of neural net (must be > 1)
    Nneurons     = list of number of neurons for each perceptron    
    ################################################################
    e.g.: a network with 2 inputs, a first layer of 2 neurons, an
          intermediate layer of 2 neurons, and an output layer of 1
          neuron 2, 3, [2, 2, 1]
    ################################################################
    """
    def __init__(self,Nvars=2,Nperceptrons=2,Nneurons=[2,1]):
        self.Nperceptrons  = Nperceptrons
        self.FFperceptrons = []
        self.inputDcDz     = []
        
        """
        ###################
        Feedforward network
        ###################
        """
        ### First layer ###
        self.FFperceptrons.append(Perceptron(Nneurons[0],Nvars))

        ### Intermediate layers ###
        self.FFperceptrons.extend(Perceptron(Nneurons[j],self.FFperceptrons[j-1].Nneurons) for j in xrange(1,self.Nperceptrons-1))

        ### Last layer ###
        self.FFperceptrons.append(Perceptron(Nneurons[self.Nperceptrons-1],self.FFperceptrons[self.Nperceptrons-2].Nneurons))

        """
        #######################
        Backpropagation network
        #######################
        """
        self.initBPnetwork()

    def initBPnetwork(self):
        self.BPperceptrons = []

        ### First layer ###
        self.BPperceptrons.append(Perceptron(self.FFperceptrons[self.Nperceptrons-2].Nneurons,self.FFperceptrons[self.Nperceptrons-1].Nneurons,'BPN'))

        ### Intermediate layers + last layer ###
        self.BPperceptrons.extend(Perceptron(self.FFperceptrons[self.Nperceptrons-j-1].Nneurons,self.FFperceptrons[self.Nperceptrons-j].Nneurons,'BPN') for j in xrange(2,self.Nperceptrons))
        
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

    def learn(self,invec,target,miniBatch=0):
        """
        ###########
        Feedforward
        ###########
        """
        self.eval(invec)

        """
        ###############
        Backpropagation
        ###############
        """
        self.copyFFtoBPweights()
        self.inputDcDz = self.FFperceptrons[self.Nperceptrons-1].dcFunDz(target)
        self.evalBP(self.inputDcDz)

        """
        ######
        Update
        ######
        """
        ### Last layer ###
        self.FFperceptrons[self.Nperceptrons-1].adapt([ N.afun for N in self.FFperceptrons[self.Nperceptrons-2].neurons ],self.inputDcDz,miniBatch)

        ### Intermediate layers ###
        if self.Nperceptrons > 2:
            for j in xrange(self.Nperceptrons-2):
                self.FFperceptrons[self.Nperceptrons-j-2].adapt([ N.afun for N in self.FFperceptrons[self.Nperceptrons-j-3].neurons ],
                                                                [ N.afun for N in self.BPperceptrons[j].neurons ],miniBatch)

        ### First layer ###
        self.FFperceptrons[0].adapt(invec,[ N.afun for N in self.BPperceptrons[self.Nperceptrons-2].neurons ],miniBatch)

        if self.FFperceptrons[0].neurons[0].regular != 0:
            return self.FFperceptrons[self.Nperceptrons-1].cFun(target) + self.FFperceptrons[0].neurons[0].regular/2 * sum(P.sum2W() for P in self.FFperceptrons)
        else:
            return self.FFperceptrons[self.Nperceptrons-1].cFun(target)

    def speed(self,who):
        return sqrt(sum(d * d for d in self.inputDcDz)) if who == 0 else self.BPperceptrons[who-1].speed()
    
    def printParams(self):
        print '\n\n===>>> Feedforward Neural Net ===>>>'
        for j,P in enumerate(self.FFperceptrons):
            print 'Perceptron[', j, ']'
            P.printParams()

        print '\n<<<=== Backpropagation Neural Net <<<==='
        for j,P in enumerate(self.BPperceptrons):
            print 'Perceptron[', j, ']'
            P.printParams()

    def reset(self):
        for P in self.FFperceptrons:
            P.reset()

        for P in self.BPperceptrons:
            P.reset()

    def scramble(self,who):
        if -1 in who.keys():
            who = { j:[-1] for j in xrange(self.Nperceptrons) }

        for j in who:
            self.FFperceptrons[j].scramble(who[j])

    def fixAllBut(self,who):
        """
        #########################################################################
        e.g.:
        NN.fixAllBut({2:[0,3]}) --> fix all but neurons 0 and 3 from perceptron 2
        #########################################################################
        """
        for j,P in enumerate(self.FFperceptrons):
            P.fixAllBut([-1]) if j not in who else P.fixAllBut(who[j])

    def release(self,who):
        """
        ###################################################################
        e.g.:
        NN.release({-1:[]})   --> release everything
        NN.release({1:[-1]})  --> release perceptron 1
        NN.release({2:[0,3]}) --> release neurons 0 and 3 from perceptron 2
        ###################################################################
        """
        if -1 in who.keys():
            who = { j:[-1] for j in xrange(self.Nperceptrons) }

        for j in who:
            self.FFperceptrons[j].release(who[j])

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

    def aFunMin(self,indx):
        return self.FFperceptrons[indx].aFunMin(0)

    def aFunMax(self,indx):
        return self.FFperceptrons[indx].aFunMax(0)

    def remove(self,who):
        """
        #################################################################
        e.g.:
        NN.remove({-1:[]})   --> remove everything
        NN.remove({1:[-1]})  --> remove perceptron 1
        NN.remove({2:[0,3]}) --> remove neurons 0 and 3 from perceptron 2
        #################################################################
        """
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
        """
        #####################################################################
        e.g.:
        NN.add({1:3})     --> add new perceptron with 3 neurons in position 1
        NN.add({2:[0,3]}) --> add neurons in position 0 and 3 in perceptron 2
        #####################################################################
        """
        for j in who:
            if type(who[j]) is list:
                self.FFperceptrons[j].addN(who[j])
            else:
                self.Nperceptrons += 1
                self.FFperceptrons.insert(j,Perceptron(who[j],self.FFperceptrons[j].neurons[0].Nvars))

            if j < self.Nperceptrons-1:
                self.FFperceptrons[j+1].addW(who[j])
                
        self.initBPnetwork()

    def setLearnRate(self,val):
        for P in self.FFperceptrons:
            P.setLearnRate(val)

    def save(self,name):
        f   = open(name,'w')
        out = []
        
        f.write('# === M.E.D. Neural Network === #\n')
        f.write('# {0:20s} {1:20s} {2:40s} {3:20s} {4:20s} {5:20s}\n'.format('Nvars','Nperceptrons','Nneurons','Learn rate','Regularization','RMSprop decay'))

        out.append(str(self.FFperceptrons[0].neurons[0].Nvars))
        out.append(str(self.Nperceptrons))
        out.append('[')
        for P in self.FFperceptrons:
            out[-1] += ' ' + str(P.Nneurons)
        out[-1] += ' ]'
        out.append(str(self.FFperceptrons[0].neurons[0].learnRate))
        out.append(str(self.FFperceptrons[0].neurons[0].regular))
        out.append(str(self.FFperceptrons[0].neurons[0].rmsPrDecay))
        f.write('  {0:20s} {1:20s} {2:40s} {3:20s} {4:20s} {5:20s}\n\n'.format(out[0],out[1],out[2],out[3],out[4],out[5]))

        for j,P in enumerate(self.FFperceptrons):
            f.write('Perceptron[ ' + str(j) + ' ]\n')
            P.save(f)

        f.close()

    def read(self,name):
        f     = open(name,'r')
        line  = f.readline()
        lele  = line.split()
        leled = []
        lelef = []
        
        while len(lele) == 0 or (len(lele) > 0 and '#' in lele[0]):
            line = f.readline()
            lele = line.split()

        for ele in lele:
            if ele.isdigit():
                leled.append(int(ele))
            elif ele.replace('.','').replace('-','').isdigit():
                lelef.append(float(ele))

        self.__init__(leled[0],leled[1],leled[2:])
        self.FFperceptrons[0].neurons[0].learnRate  = lelef[0]
        self.FFperceptrons[0].neurons[0].regular    = lelef[1]
        self.FFperceptrons[0].neurons[0].rmsPrDecay = lelef[2]

        for P in self.FFperceptrons:
            P.read(f)

        f.close()

    def saveHypPar(self,*args):
        """
        ###################################
        Calling sequence:
        args[0] = file name
        args[1] = number of runs
        args[2] = mini batch
        args[3] = learn rate starting value
        args[4] = learn rate ending value
        args[5] = learn rate decay time
        args[6] = scramble True/False
        args[7] = scrambled neurons
        ###################################
        """
        f = open(args[0],'a')

        f.write('\n# Hyper-parameter information\n')
        f.write('  Number of runs: {:_^20d}\n'.format(args[1]))
        f.write('  Mini-batch: {:_^20d}\n'.format(args[2]))

        f.write('\n  # Learn rate parameters\n')
        f.write('  Learn rate start: {:_^20f}\n'.format(args[3]))
        f.write('  Learn rate end: {:_^20f}\n'.format(args[4]))
        f.write('  Learn rate tau: {:_^20d}\n'.format(args[5]))

        f.write('\n  # Scramble parameters\n')
        f.write('  Neuron to scramble: {:_^20s}\n'.format(str(args[7]))) if args[6] == True else f.write('  Neuron to scramble: {:_^20s}\n'.format(''))

        f.close()

    def readHypPar(self,name):
        """
        ############################
        Return sequence:
        1. number of runs
        2. mini batch
        3. learn rate starting value
        4. learn rate ending value
        5. learn rate decay time
        6. scrambled neurons
        ############################
        """
        f    = open(name,'r')
        line = f.readline()

        while 'Hyper-parameter' not in line:
            line = f.readline()

        it = 0
        while f.tell() != os.fstat(f.fileno()).st_size:
            line = f.readline()
            lele = line.split('_')
            for ele in lele:
                if ele.isdigit():
                    yield int(ele)
                    it += 1
                elif ele.replace('.','').replace('-','').isdigit():
                    yield float(ele)
                    it += 1
                elif '{' and '}' in ele:
                    yield ast.literal_eval(ele)
                    it += 1

        f.close()
        
        while it < 6:
            yield
            it += 1
