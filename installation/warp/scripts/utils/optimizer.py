"""This file contains several optimizers, including:
  Spsa: Simultaneaous Perturbation Stochastic Approximation
  Evolution: Genetic based algorithm
  ParticleSwarm: Particle swarm optimizer
  Simpleoptimizer: simple search
"""

from numpy import *
from numpy import random

try:
    # --- Import threading for the ParticleSwarm optimizer
    import threading
except ImportError:
    pass

# --- Convenience function
def convertinputtoarray(p,n,default,name):
    """
Given a function input value, convert it to an array.
 - p: the input value
 - n: the length of the resulting array
 - default: default value
    """
    if p is None:
        result = default*ones(n)
    elif callable(p):
        result = p
    else:
        try:
            if len(p) == n:
                result = array(p)
            else:
                raise Exception("ERROR: %s does not have the correct length"%name)
        except TypeError: # p has no len (i.e. it is a scalar)
            result = ones(n)*p
    return result

class Spsa:
    """
Implements the Simultaneaous Perturbation Stochastic Approximation
minimization algorithm. To use, create an instance of the Spsa class
and then call the iter method.

opt = Spsa(...)
opt.iter(...)

Note that the parameters are all varied by approximately the same amount and so
should be normalized to be of the same order of magnitude.

Constructor arguments:
  - nparams: Number of params to vary
  - params: Initial values of the parameters
  - func: Method which does any calculations given the varying parameters.
          It takes a single argument, the list of parameters.
  - lossfunc: Method which calculates the loss and returns it.
  - c1: Amount by which params are initially varied to calculate gradient
  - a1: Amount by which the params are changed, scaled by the gradient (the gradient
        is the change in loss when params are changed by c1, divided by c1).
        To check the scaling, call the method gradloss, which returns the gradient.
  - a2=100.: Scale (in iteration numbers) over which c and a are decreased
  - paramsmin=-1.e+36: Min value of the parameters. This can be a function
                       which takes the params as its single argument and returns
                       the parameter mins.
  - paramsmax=+1.e+36: Max value of the parameters. This can be a function
                       which takes the params as its single argument and returns
                       the parameter maxes.
  - paramsave=None: average value of each parameter, used to scale params
                    internally to improve performance
  - paramsrms=None: RMS values of each parameter, used to scale params
                    internally to improve performance
  - verbose=False: when true, print diagnostics
  - errmax=+1.e36: Maximum acceptable value of the error. If the error
                   is greater, the iteration is skipped.
  - saveparamhist=False: When true, saves the history of the parameters in the
                         attribute hparam. Note that the history of the loss is
                         always saved in hloss.
  - picklehistfile=None: When given, save the history of the loss and the
                         parameters in the given pickle file.

Methods:
  iter: carry out the iterations.
  gradloss: return the gradient of the loss function
    """

    def __init__(self,nparams,params,func,lossfunc,c1,a1,a2=100.,
                 paramsmin=None,paramsmax=None,paramsave=None,paramsrms=None,
                 verbose=False,errmax=1.e36,
                 saveparamhist=False,picklehistfile=None):
        """
    Creates an instance of the Spsa class.
        """
        self.nparams = nparams
        self.params = array(params)
        self.func = func
        self.loss = lossfunc
        self.k = 1
        self.a1 = a1
        self.a2 = a2
        self.c1 = c1
        self.verbose = verbose
        self.hloss = []
        self.errmax = errmax
        self.saveparamhist = saveparamhist
        self.paramsmin = convertinputtoarray(paramsmin,self.nparams,-1.e+36,'paramsmin')
        self.paramsmax = convertinputtoarray(paramsmax,self.nparams,+1.e+36,'paramsmax')
        self.paramsave = paramsave
        self.paramsrms = paramsrms
        self.params = self.scaledparams(self.params)
        if self.saveparamhist: self.hparam = []
        self.picklehistfile = picklehistfile
    def ak(self):
        return self.a1/(self.k+self.a2)**0.602
    def ck(self):
        return self.c1/self.k**.101
    def scaledparams(self,params=None):
        if params is None: params = self.params
        if self.paramsave is not None: params = params - self.paramsave
        if self.paramsrms is not None: params = params/self.paramsrms
        return params
    def unscaledparams(self,params=None):
        if params is None: params = self.params
        if self.paramsrms is not None: params = params*self.paramsrms
        if self.paramsave is not None: params = params + self.paramsave
        return params
    def gradloss(self):
    # --- Calculated the approximated gradient of the loss function.
        deltak = (2*random.random(self.nparams)).astype(long) - .5
        nextparams = self.constrainparams(self.params + self.ck()*deltak)
        nextparams = self.unscaledparams(nextparams)
        self.func(nextparams)
        lplus = self.loss()
        nextparams = self.constrainparams(self.params - self.ck()*deltak)
        nextparams = self.unscaledparams(nextparams)
        self.func(nextparams)
        lminus = self.loss()
        return (lplus - lminus)/(2.*self.ck()*deltak)
    def printerror(self,err):
        print "iterations = %d  error = %e %e" %(self.k,self.loss(),self.loss()/err)
    def printparams(self):
        pp = self.unscaledparams(self.params)
        for i in range(self.nparams): print '%15.12e'%pp[i]
    def getparamsmin(self,params):
        """Returns the min limit of parameters."""
        if callable(self.paramsmin):
            return self.paramsmin(self.unscaledparams(params))
        return self.paramsmin
    def getparamsmax(self,params):
        """Returns the max limit of parameters."""
        if callable(self.paramsmax):
            return self.paramsmax(self.unscaledparams(params))
        return self.paramsmax
    def constrainparams(self,params):
        "Makes sure all params are within bounds"
        params = self.unscaledparams(params)
        params = maximum(params,self.getparamsmin(params))
        params = minimum(params,self.getparamsmax(params))
        params = self.scaledparams(params)
        return params

    def iter(self,err=1.e-9,imax=100,verbose=None,kprint=10,kprintlogmax=3):
        """
Function to do iterations.
  - err=err=1.e-9: Convergence criteria
  - imax=10000: Maximum number of iterations
  - verbose=None: can override main value of verbose
  - kprint=10: exponential scale for printing out current status
  - kprintlogmax=3: max value of exponential scale
        """
        if verbose is None: verbose = self.verbose
        i = 0
        while (self.loss()>err and i < imax):
            i = i + 1
            # --- calculate new parameters
            #self.params = self.params - self.ak()*self.gradloss()
            dp = self.gradloss()
            if verbose:
                print "ak = %f" % self.ak()
                print "gradloss = " + repr(dp)
                print "params = " + repr(self.unscaledparams(self.params))
            oldparams = self.params + 0.
            self.params = self.constrainparams(self.params - self.ak()*dp)
            if verbose: print "new params = " + repr(self.unscaledparams(self.params))
            # --- Calculate function with new params
            if self.saveparamhist:
                self.hparam.append(self.unscaledparams(self.params))
            self.func(self.unscaledparams(self.params))
            # --- Check if loss it too great.
            if self.loss() > self.errmax:
                self.params = oldparams
                if verbose: print "Skipping"
                continue
            # --- Save the latest value of the loss function.
            latestloss = self.loss()
            self.hloss.append(latestloss)
            # --- Increment the counter
            self.k = self.k + 1
            # --- Print out loss value
            if kprint > 1:
                klog = int(log(self.k)/log(kprint))
            else:
                klog = kprintlogmax
            klog = min(klog,kprintlogmax)
            if ((self.k%(kprint**klog)) == 0):
                self.printerror(err)
                if klog == kprintlogmax: self.printparams()
            if self.picklehistfile is not None:
                with open(self.picklehistfile,'ab') as ff:
                    cPickle.dump([self.k,self.unscaledparams(),latestloss],ff,-1)
        # --- Print out the resulting params
        self.printerror(err)
        self.printparams()



###########################################################################
###########################################################################
###########################################################################
# Performs global optimization using genetic evolution algorithms.
# Algorithm taken from Dr Dobb's Journal, April 1997, K Price, R. Storn
###########################################################################

class Evolution:
    """
Differential Evolution
  Performs global optimization using genetic evolution algorithms.
  Algorithm taken from Dr Dobb's Journal, April 1997, K Price, R. Storn

Input:
  - npop: size of population (must be greater than 3)
  - nparams: number of parameters
  - evaluate(params): is function, given a set of parameters, returns a score
  - params: intial set of parameters
  - deltas=0.01: fractional variation of parameters to fill initial population
  - shifts=0.: absolute variation of parameters to fill initial population
  - crossover=0.5: fraction of crossovers, in range [0,1)
  - f=0.7: differential factor, in range (0,1.2]
  - paramsmin=-1.e+36: Min value of the parameters. This can be a function
                       which takes the params as its single argument.
  - paramsmax=+1.e+36: Max value of the parameters. This can be a function
                       which takes the params as its single argument.

Methods:
  evolve: carry out the iteration
  best_params: returns parameters which give the lowest score

    """
    def __init__(self,npop,nparams,evaluate,params,deltas=None,shifts=None,
                 crossover=.5,f=.7,paramsmin=None,paramsmax=None):
        """
Differential Evolution
        """
        self.npop = npop
        self.nparams = nparams
        self.deltas = deltas
        self.shifts = shifts
        self.initparams = params
        self.crossover = crossover
        self.f = f
        self.evaluate = evaluate
        if (crossover < 0. or crossover > 1.):
            print "Warning: crossover outside of the range [0,1)"
        if (f < 0. or f > 1.2):
            print "Warning: differential factor f outside of the range (0,1.2]"
        if (npop < 4):
            raise Exception("Error: number of populations, npop, must be greater than 3")
        self.trial = zeros(nparams,'d')
        self.x1 = zeros((npop,nparams),'d')
        self.x2 = zeros((npop,nparams),'d')
        self.cost = zeros(npop,'d')
        self.count = 0
        self.paramsmin = convertinputtoarray(paramsmin,self.nparams,-1.e+36,'paramsmin')
        self.paramsmax = convertinputtoarray(paramsmax,self.nparams,+1.e+36,'paramsmax')
        self.linitialized = False
    def best_params(self):
        "Function to return best set of parameters so far"
        imin = 0
        costmin = self.cost[imin]
        for i in range(1,self.npop):
            if self.cost[i] <  costmin:
                imin = i
                costmin = self.cost[i]
        return self.x1[imin,:]
    def printbestcost(self):
        print "Generation %d, best cost %f worst cost %f"% \
              (self.count,min(self.cost),max(self.cost))
    def getparamsmin(self,params):
        """Returns the min limit of parameters."""
        if callable(self.paramsmin): return self.paramsmin(params)
        return self.paramsmin
    def getparamsmax(self,params):
        """Returns the max limit of parameters."""
        if callable(self.paramsmax): return self.paramsmax(params)
        return self.paramsmax
    def constrainparams(self,params):
        "Makes sure all params are within bounds"
        params = maximum(params,self.getparamsmin(params))
        params = minimum(params,self.getparamsmax(params))
        return params

    def evolve_reset(self):
        "Reset cost function.  Used for example if cost function is changed."
        for i in range(npop):
            cost[i] = evaluate(x1[i,:])

    def evolve_init(self):
        """
Function to initialize the population.
Picks parameters randomly distributed by deltas and shifts about a base
sample set of parameters.
  - sample: is the initial set of parameters
  - delta=0.01: is the fractional variation about the sample
                It can either be a scalar or an array the same size as sample.
  - shifts=0.: is the absolute variation about the sample.
               It can either be a scalar or an array the same size as sample.
        """
        if self.linitialized: return
        self.linitialized = True

        sample = self.initparams
        deltas = convertinputtoarray(self.deltas,self.nparams,.01,'deltas')
        shifts = convertinputtoarray(self.shifts,self.nparams,.0,'shifts')
        self.deltas = deltas
        self.shifts = shifts

        for i in range(self.npop):
            # --- The first one uses the sample as given.
            trial = sample
            if i > 0:
                # --- The rest use a perturbation from the sample.
                trial = (trial*(1.+2.*(random.random(self.nparams)-.5)*deltas)
                              + 2.*(random.random(self.nparams)-.5)*shifts)
            self.x1[i,:] = self.constrainparams(trial)
            self.cost[i] = self.evaluate(self.x1[i,:])

    def evolve(self,gen_max=1,nprint=100):
        """
Do the optimization
  - gen_max=1: number of generations to run through
  - nprint=100: base frequency to print cost
        """

        self.evolve_init()

        self.score = self.cost[0]

        # --- Loop over the generations
        for count in range(gen_max):
            self.count = self.count + 1

            # --- loop through population
            for i in range(self.npop):

                # Mutate/Recombine

                # --- Randomly pick three vectors different from each other and 'i'.
                a = i
                b = i
                c = i
                while (a == i):                     a = int(random.random()*self.npop)
                while (b == i or b == a):           b = int(random.random()*self.npop)
                while (c == i or c == a or c == b): c = int(random.random()*self.npop)

                # --- Randomly pick the first parameter
                j = int(random.random()*self.nparams)

                # --- Load parameters into trial, performing binomial trials
                for k in range(self.nparams):
                    if (random.random() < self.crossover or k == self.nparams-1):
                        # --- Source for trial is a random vector plus weighted differential
                        # --- The last parameter always comes from noisy vector
                        self.trial[j] = self.x1[c,j] + self.f*(self.x1[a,j] - self.x1[b,j])
                    else:
                        # --- Trial parameter come from x1 itself
                        self.trial[j] = self.x1[i,j]
                    # --- get next parameter
                    j = (j+1)%self.nparams

                # Evaluate/Select

                # --- Evaluate trial function
                self.trial = self.constrainparams(self.trial)
                self.score = self.evaluate(self.trial)

                if (self.score <= self.cost[i]):
                    # --- If trial improves on x1, move trial to secondary array
                    # --- and save the new score
                    self.x2[i,:] = self.trial
                    self.cost[i] = self.score
                else:
                    # --- otherwise move x1 to secondary array
                    self.x2[i,:] = self.x1[i,:]

            # --- End of population loop, so copy new parameters into x1
            self.x1[...] = self.x2[...]

            # --- Print out loss function
            if (self.count <= nprint):
                self.printbestcost()
            elif ((self.count>nprint) and (self.count<=nprint**2) and
                  ((self.count%nprint)==0)):
                self.printbestcost()
            elif ( (self.count%(nprint**2)) == 0):
                self.printbestcost()
                print self.best_params()

        self.printbestcost()
        print self.best_params()

    ### --- threaded code below --- ###
    def initializememberthread(self,i):
        # --- The first one uses the initial params as given.
        trial = self.initparams
        if i > 0:
            # --- The rest use a perturbation from the sample.
            trial = (trial*(1.+2.*(random.random(self.nparams)-.5)*self.deltas)
                          + 2.*(random.random(self.nparams)-.5)*self.shifts)

        self.x1[i,:] = self.constrainparams(trial)
        self.threadthrottle.acquire()
        self.cost[i] = self.evaluate(self.x1[i,:])
        self.threadthrottle.release()

    def evolve_initthread(self):
        """
Function to initialize the population.
Picks parameters randomly distributed by deltas and shifts about a base
sample set of parameters.
  - sample: is the initial set of parameters
  - delta=0.01: is the fractional variation about the sample
                It can either be a scalar or an array the same size as sample.
  - shifts=0.: is the absolute variation about the sample.
               It can either be a scalar or an array the same size as sample.
        """
        if self.linitialized: return

        sample = self.initparams
        deltas = convertinputtoarray(self.deltas,self.nparams,.01,'deltas')
        shifts = convertinputtoarray(self.shifts,self.nparams,.0,'shifts')
        self.deltas = deltas
        self.shifts = shifts

        initthreads = []
        for i in range(self.npop):
            print "starting init thread ",i
            initthreads.append(
              threading.Thread(target=self.initializememberthread,
                               name='init%d'%i,
                               args=(i,)))
            initthreads[-1].start()

        # --- Wait for the threads to finish
        for t in initthreads:
            t.join()
        #while threading.active_count() > 1:
        #  print threading.active_count()

        self.linitialized = True

    def evolvememberthread(self,i):

        # Mutate/Recombine

        # --- Randomly pick three vectors different from each other and 'i'.
        a = i
        b = i
        c = i
        while (a == i):                     a = int(random.random()*self.npop)
        while (b == i or b == a):           b = int(random.random()*self.npop)
        while (c == i or c == a or c == b): c = int(random.random()*self.npop)

        # --- Randomly pick the first parameter
        j = int(random.random()*self.nparams)

        # --- Load parameters into trial, performing binomial trials
        trial = zeros_like(self.x1[i,:])
        for k in range(self.nparams):
            if (random.random() < self.crossover or k == self.nparams-1):
                # --- Source for trial is a random vector plus weighted differential
                # --- The last parameter always comes from noisy vector
                trial[j] = self.x1[c,j] + self.f*(self.x1[a,j] - self.x1[b,j])
            else:
                # --- Trial parameter come from x1 itself
                trial[j] = self.x1[i,j]
            # --- get next parameter
            j = (j+1)%self.nparams

        # Evaluate/Select

        # --- Evaluate trial function
        trial = self.constrainparams(trial)
        self.threadthrottle.acquire
        score = self.evaluate(trial)
        self.threadthrottle.release

        if (score <= self.cost[i]):
            # --- If trial improves on x1, move trial to secondary array
            # --- and save the new score
            self.x2[i,:] = trial
            self.cost[i] = score
        else:
            # --- otherwise move x1 to secondary array
            self.x2[i,:] = self.x1[i,:]

    def evolvethread(self,gen_max=1,nprint=100,maxthreads=1000000):
        """
Do the optimization
  - gen_max=1: number of generations to run through
  - nprint=100: base frequency to print cost
  - maxthreads=1000000: The number of threads that will be started will be the
                        minimum of maxthreads and the size of the population
        """

        self.threadthrottle = threading.Semaphore(maxthreads)

        self.evolve_initthread()

        # --- Loop over the generations
        for count in range(gen_max):
            self.count = self.count + 1

            # --- loop through population
            iterthreads = []
            for i in range(self.npop):
                print "starting thread ",i
                iterthreads.append(
                  threading.Thread(target=self.evolvememberthread,
                                   name='iter%d'%i,
                                   args=(i,)))
                iterthreads[-1].start()

            # --- Wait for the threads to finish
            for t in iterthreads:
                t.join()
            #while threading.active_count() > 1:
            #  print threading.active_count()

            # --- End of population loop, so copy new parameters into x1
            self.x1[...] = self.x2[...]

            # --- Print out loss function
            if (self.count <= nprint):
                self.printbestcost()
            elif ((self.count>nprint) and (self.count<=nprint**2) and
                  ((self.count%nprint)==0)):
                self.printbestcost()
            elif ( (self.count%(nprint**2)) == 0):
                self.printbestcost()
                print self.best_params()

        self.printbestcost()
        print self.best_params()


###########################################################################

class ParticleSwarm:
    """
Particle swarm optimization
  Performs global optimization
  See for instance http://en.wikipedia.org/wiki/Particle_swarm_optimization

Input:
  - npop: size of population. A larger population will give better convergence
          though probably more function evaluations.
  - evaluate(params): is function, given a set of parameters, returns a score.
                      Its input argument is a list of the parameters.
  - initparams: intial set of parameters
  - deltas=0.01: fractional variation of parameters to fill initial population
  - shifts=0.: absolute variation of parameters to fill initial population
  - decel=1.: Scale factor on the particle "velocity"
  - cognitiveweight=2.: How much weight to give to the particle's best case
  - socialweight=2.: How much weight to give to the neighborhood's best case
  - globalweight=0.: How much weight to give to the global best case
  - neighbors=2: Number of neighbors to include in the neighborhoods
  - decelfinal=0.4: Final value that decel is reduced to during the iterations
  - coolingrate=100.: Number of iterations over which the decel is reduced
                      from its initial value to its final value.
  - paramsmin=-1.e+36: Min value of the parameters. This can be a function
                       which takes the params as its single argument.
  - paramsmax=+1.e+36: Max value of the parameters. This can be a function
                       which takes the params as its single argument.

Methods:
  swarm(niters=1,nprint=100): do the minimization
  swarmthread(niters=1,nprint=100,maxthreads=1000000): do the minimization,
      carrying out the function evaluations in separate threads

    """
    def __init__(self,npop,evaluate,initparams,deltas=None,shifts=None,
                 decel=1.,cognitiveweight=2.,socialweight=2.,globalweight=0.,
                 neighbors=2,decelfinal=0.4,coolingrate=100.,
                 paramsmin=None,paramsmax=None):
        self.npop = npop
        self.evaluate = evaluate
        self.initparams = initparams
        self.deltas = deltas
        self.shifts = shifts
        self.nparams = len(initparams)
        self.decel = decel
        self.cognitiveweight = cognitiveweight
        self.socialweight = socialweight
        self.globalweight = globalweight
        self.neighbors = neighbors
        self.decelfinal = decelfinal
        self.coolingrate = coolingrate

        self.globalbestparams = zeros(self.nparams,'d')
        self.globalbestcost = inf

        # --- Total iteration count
        self.count = 0

        self.paramsmin = convertinputtoarray(paramsmin,self.nparams,-1.e+36,'paramsmin')
        self.paramsmax = convertinputtoarray(paramsmax,self.nparams,+1.e+36,'paramsmax')

        self.setup_neighborhoods()
        self.linitialized = False

    def reset(self):
        """Resets the optimizer to start with a new population, starting from
the previous best global parameters. The deceleration is also reset."""
        self.initparams = self.globalbestparams.copy()
        self.linitialized = False
        self.count = 0

    def getdecel(self):
        """Handles the cooling rate on decel factor"""
        if self.count <= self.coolingrate:
            decel = (self.decel - (self.decel - self.decelfinal)*
                                  (self.count-1)/self.coolingrate)
        else:
            decel = self.decelfinal
        return decel

    def printbestcost(self):
        print "Generation %d, global best cost %e, best cost %e, worst cost %e"% \
              (self.count,self.globalbestcost,min(self.bestcost),max(self.bestcost))
        print "Global best params ",self.globalbestparams
    def getparamsmin(self,params):
        """Returns the min limit of parameters."""
        if callable(self.paramsmin): return self.paramsmin(params)
        return self.paramsmin
    def getparamsmax(self,params):
        """Returns the max limit of parameters."""
        if callable(self.paramsmax): return self.paramsmax(params)
        return self.paramsmax
    def constrainparams(self,params):
        "Makes sure all params are within bounds"
        params = maximum(params,self.getparamsmin(params))
        params = minimum(params,self.getparamsmax(params))
        return params

    def setup_neighborhoods(self):
        """Setup the neighbors for each particle.
        This uses a simple circle neighborhood"""
        self.neighborhoods = []
        halfneigh = self.neighbors//2
        for i in range(self.npop):
            self.neighborhoods.append([])
            for n in range(-halfneigh,halfneigh+1):
                self.neighborhoods[-1].append((i+n)%self.npop)

    def findbestneighbor(self,neighborhood):
        """Searches through the given neighborhood, finding the case with the best cost."""
        bestparams = self.bestparams[neighborhood[0]]
        bestcost = self.bestcost[neighborhood[0]]
        for n in neighborhood[1:]:
            if self.bestcost[n] < bestcost:
                bestparams = self.bestparams[n]
                bestcost = self.bestcost[n]
        return bestparams

    def init_population(self):
        """
Initializes the population.
Picks parameters randomly distributed by deltas and shifts about the initial
set of parameters.
  - initparams: is the initial set of parameters
  - delta=0.01: is the fractional variation
                It can either be a scalar or an array
  - shifts=0.: is the absolute variation
               It can either be a scalar or an array
        """
        if self.linitialized: return
        deltas = convertinputtoarray(self.deltas,self.nparams,.01,'deltas')
        shifts = convertinputtoarray(self.shifts,self.nparams,.0,'shifts')
        self.deltas = deltas
        self.shifts = shifts

        self.trial = zeros((self.npop,self.nparams),'d')
        self.velocity = zeros((self.npop,self.nparams),'d')
        self.cost = zeros(self.npop,'d')
        self.bestparams = zeros((self.npop,self.nparams),'d')
        self.bestcost = zeros(self.npop,'d')

        for i in range(self.npop):
            self.initializeparticle(i)

        self.linitialized = True

    def initializeparticle(self,i):

        # --- Start with the initial parameters
        trial = self.initparams
        if i > 0:
            # --- All but the first particle add a random perturbation from
            # --- the initial parameters.
            trial = (trial*(1.+2.*(random.random(self.nparams)-.5)*self.deltas)
                         + 2.*(random.random(self.nparams)-.5)*self.shifts)

        self.trial[i,:] = self.constrainparams(trial)
        self.cost[i] = self.evaluate(self.trial[i,:])
        self.bestparams[i,:] = self.trial[i,:]
        self.bestcost[i] = self.cost[i]

        if self.bestcost[i] < self.globalbestcost:
            self.globalbestcost = self.bestcost[i]
            self.globalbestparams[:] = self.bestparams[i,:]

    def updateparticle(self,i):

        # --- Update the velocity and trial
        r1 = random.random(self.nparams)
        r2 = random.random(self.nparams)
        r3 = random.random(self.nparams)
        bestneighbor = self.findbestneighbor(self.neighborhoods[i])
        if all(bestneighbor==self.bestparams[i,:]): r3 = 0.
        self.velocity[i,:] = (self.getdecel()*self.velocity[i,:] +
              r1*self.globalweight*(self.globalbestparams - self.trial[i,:]) +
              r2*self.cognitiveweight*(self.bestparams[i,:] - self.trial[i,:]) +
              r3*self.socialweight*(bestneighbor - self.trial[i,:]))
        self.trial[i,:] += self.velocity[i,:]
        self.trial[i,:] = self.constrainparams(self.trial[i,:])

        # --- Evaluate trial function
        self.cost[i] = self.evaluate(self.trial[i,:])

        if self.cost[i] < self.bestcost[i]:
            self.bestparams[i,:] = self.trial[i,:]
            self.bestcost[i] = self.cost[i]

        if self.cost[i] < self.globalbestcost:
            self.globalbestparams[:] = self.trial[i,:]
            self.globalbestcost = self.cost[i]

    def swarm(self,niters=1,nprint=100):
        """
Do the optimization
  - niters=1: number of iterations to run through
  - nprint=100: base frequency to print cost
        """

        self.init_population()

        for count in range(niters):
            self.count += 1

            # --- loop through the particles
            for i in range(self.npop):
                self.updateparticle(i)

            # --- Print out loss function
            if nprint > 0:
                if (self.count <= nprint):
                    self.printbestcost()
                elif ((self.count>nprint) and (self.count<=nprint**2) and
                      ((self.count%nprint)==0)):
                    self.printbestcost()
                elif ( (self.count%(nprint**2)) == 0):
                    self.printbestcost()

        if nprint > 0:
            self.printbestcost()
            print self.globalbestparams

    ##### --- threaded code below here --- #####
    def findbestneighborthread(self,neighborhood):
        self.bestparamslock.acquire()
        bestparams = self.bestparams[neighborhood[0]]
        bestcost = self.bestcost[neighborhood[0]]
        for n in neighborhood[1:]:
            if self.bestcost[n] < bestcost:
                bestparams = self.bestparams[n]
                bestcost = self.bestcost[n]
        self.bestparamslock.release()
        return bestparams

    def initializeparticlethread(self,i):

        # --- Start with the initial parameters
        trial = self.initparams
        if i > 0:
            # --- All but the first particle add a random perturbation from
            # --- the initial parameters.
            trial = (trial*(1.+2.*(random.random(self.nparams)-.5)*self.deltas)
                          +2.*(random.random(self.nparams)-.5)*self.shifts)

        self.trial[i,:] = self.constrainparams(trial)
        self.threadthrottle.acquire()
        self.cost[i] = self.evaluate(self.trial[i,:])
        self.threadthrottle.release()
        self.bestparams[i,:] = self.trial[i,:]
        self.bestcost[i] = self.cost[i]

        self.globalbestparamslock.acquire()
        if self.bestcost[i] < self.globalbestcost:
            self.globalbestcost = self.bestcost[i]
            self.globalbestparams[:] = self.bestparams[i,:]
        self.globalbestparamslock.release()

    def updateparticlethread(self,i):

        # --- Update the velocity and trial
        r1 = random.random(self.nparams)
        r2 = random.random(self.nparams)
        r3 = random.random(self.nparams)
        bestneighbor = self.findbestneighbor(self.neighborhoods[i])
        if all(bestneighbor==self.bestparams[i,:]): r3 = 0.
        self.velocity[i,:] = (self.getdecel()*self.velocity[i,:] +
              r1*self.globalweight*(self.globalbestparams - self.trial[i,:]) +
              r2*self.cognitiveweight*(self.bestparams[i,:] - self.trial[i,:]) +
              r3*self.socialweight*(bestneighbor - self.trial[i,:]))
        self.trial[i,:] += self.velocity[i,:]

        # --- Evaluate trial function
        self.trial[i,:] = self.constrainparams(self.trial[i,:])
        self.threadthrottle.acquire()
        self.cost[i] = self.evaluate(self.trial[i,:])
        self.threadthrottle.release()

        if self.cost[i] < self.bestcost[i]:
            self.bestparamslock.acquire()
            self.bestparams[i,:] = self.trial[i,:]
            self.bestcost[i] = self.cost[i]
            self.bestparamslock.release()

        self.globalbestparamslock.acquire()
        if self.cost[i] < self.globalbestcost:
            self.globalbestparams[:] = self.trial[i,:]
            self.globalbestcost = self.cost[i]
        self.globalbestparamslock.release()

    def init_populationthread(self):
        """
Function to initialize the population.
Picks parameters randomly distributed by deltas and shifts about the initial
set of parameters.
  - initparams: is the initial set of parameters
  - delta=0.01: is the fractional variation
                It can either be a scalar or an array
  - shifts=0.: is the absolute variation
               It can either be a scalar or an array
        """
        if self.linitialized: return
        deltas = convertinputtoarray(self.deltas,self.nparams,.01,'deltas')
        shifts = convertinputtoarray(self.shifts,self.nparams,.0,'shifts')
        self.deltas = deltas
        self.shifts = shifts

        self.trial = zeros((self.npop,self.nparams),'d')
        self.velocity = zeros((self.npop,self.nparams),'d')
        self.cost = zeros(self.npop,'d')
        self.bestparams = zeros((self.npop,self.nparams),'d')
        self.bestcost = zeros(self.npop,'d')

        initthreads = []
        for i in range(self.npop):
            initthreads.append(
              threading.Thread(target=self.initializeparticlethread,
                               name='init%d'%i,
                               args=(i,)))
            initthreads[-1].start()

        # --- Wait for the threads to finish
        for t in initthreads:
            t.join()
        #while threading.active_count() > 1:
        #  print threading.active_count()

        self.linitialized = True

    def swarmthread(self,niters=1,nprint=100,maxthreads=1000000):
        """
Do the optimization
  - niters=1: number of iterations to run through
  - nprint=100: base frequency to print cost
  - maxthreads=1000000: The number of threads that will be started will be the
                        minimum of maxthreads and the size of the population
        """

        self.bestparamslock = threading.RLock()
        self.globalbestparamslock = threading.RLock()
        self.threadthrottle = threading.Semaphore(maxthreads)

        self.init_populationthread()

        for count in range(niters):
            self.count += 1

            # --- loop through the particles, starting a thread for each one
            iterthreads = []
            for i in range(self.npop):
                iterthreads.append(
                  threading.Thread(target=self.updateparticlethread,
                                   name='iter%d'%i,
                                   args=(i,)))
                iterthreads[-1].start()

            # --- Wait for the threads to finish
            for t in iterthreads:
                t.join()
            #while threading.active_count() > 1:
            #  print threading.active_count()

            # --- Print out loss function
            if (self.count <= nprint):
                self.printbestcost()
            elif ((self.count>nprint) and (self.count<=nprint**2) and
                  ((self.count%nprint)==0)):
                self.printbestcost()
            elif ( (self.count%(nprint**2)) == 0):
                self.printbestcost()

        self.printbestcost()
        print self.globalbestparams


##############################################################################
class Simpleoptimizer:
    """
Simple optimization over each parameter. Simple means simplistic algorithm
not ease of use. This is probably not very robust.
    """
    def pxpone(self,id):
        return self.params[id] + self.x[id,+1]*abs(self.params[id])
    def pxmone(self,id):
        return self.params[id] + self.x[id,-1]*abs(self.params[id])
    def paramspone(self,id):
        result = self.params + 0.
        result[id] = self.pxpone(id)
        return result
    def paramsmone(self,id):
        result = self.params + 0.
        result[id] = self.pxmone(id)
        return result
    def checklimits(self,id):
        if self.pxmone(id) < self.paramsmin[id]:
            self.x[id,-1] = self.paramsmin[id]/self.params[id] - 1. + 1.e-14
        if self.pxpone(id) > self.paramsmax[id]:
            self.x[id,+1] = self.paramsmax[id]/self.params[id] - 1. - 1.e-14
        self.x[id,+1] = min(self.x[id,+1],+10.*self.vary)
        self.x[id,-1] = max(self.x[id,-1],-10.*self.vary)
        self.x[id,+1] = max(self.x[id,+1],+1.e-14)
        self.x[id,-1] = min(self.x[id,-1],-1.e-14)
        #if self.params[id] < self.paramsmin[id] or self.params[id] > self.paramsmax[id]:
            #print id
            #raise Exception("Params out of bounds")
    def __init__(self,params,func,loss,vary=0.01,paramsmin=None,paramsmax=None,
                 maxxdisparity=1.e5 ):
        self.nparams = len(params)
        self.params = params
        self.func = func
        self.loss = loss
        self.vary = vary
        self.maxxdisparity = maxxdisparity
        self.paramsmin = convertinputtoarray(paramsmin,self.nparams,-1.e+36,'paramsmin')
        self.paramsmax = convertinputtoarray(paramsmax,self.nparams,+1.e+36,'paramsmax')
        for i in range(self.nparams):
            if not (self.paramsmin[i] < self.params[i] < self.paramsmax[i]):
                raise Exception("ERROR: Starting value is outside the parameter limits")
                return
        self.x = zeros((self.nparams,3),'d')
        self.f = zeros((self.nparams,3),'d')
        self.func(self.params)
        self.f[:,0] = self.loss()
        for i in range(self.nparams):
            self.x[i,-1] =  - vary
            self.checklimits(i)
            self.func(self.paramsmone(i))
            self.f[i,-1] = self.loss()
            self.x[i,+1] = + vary
            self.checklimits(i)
            self.func(self.paramspone(i))
            self.f[i,+1] = self.loss()
    def minimize1d(self,id,niters=1):
        for ii in range(niters):
            self.func(self.paramspone(id))
            self.f[id,+1] = self.loss()
            self.func(self.params)
            self.f[id,0] = self.loss()
            self.func(self.paramsmone(id))
            self.f[id,-1] = self.loss()
            #print self.params
            print "loss = %e"%(self.f[id,0])
            if self.f[id,-1] < self.f[id,0] < self.f[id,1]:
                self.x[id,+1] = - self.x[id,-1] + self.x[id,+1]
                self.params[id] = self.pxmone(id)
                self.x[id,-1] = 2.*self.x[id,-1]
                self.checklimits(id)
                self.f[id,0] = self.f[id,-1]
                self.func(self.paramsmone(id))
                self.f[id,-1] = self.loss()
            elif self.f[id,-1] > self.f[id,0] > self.f[id,1]:
                self.x[id,-1] = - self.x[id,+1] + self.x[id,-1]
                self.params[id] = self.params[id]*(1.+self.x[id,+1])
                self.x[id,+1] = 2.*self.x[id,+1]
                self.checklimits(id)
                self.f[id,0] = self.f[id,+1]
                self.func(self.paramspone(id))
                self.f[id,+1] = self.loss()
            else:
                slopem1 = (self.f[id,-1] - self.f[id,0])/abs(self.x[id,-1])
                slopep1 = (self.f[id,+1] - self.f[id,0])/abs(self.x[id,+1])
                if slopem1 < slopep1:
                    self.x[id,-1] = self.x[id,-1]/2.
                    self.params[id] = self.pxmone(id)
                    self.func(self.params)
                    self.f[id,0] = self.loss()
                    self.func(self.paramspone(id))
                    self.f[id,+1] = self.loss()
                else:
                    self.x[id,+1] = self.x[id,+1]/2.
                    self.params[id] = self.pxpone(id)
                    self.func(self.params)
                    self.f[id,0] = self.loss()
                    self.func(self.paramsmone(id))
                    self.f[id,-1] = self.loss()
    def rebalancex(self):
        # --- Make sure that the difference in the sizes of x doesn't become
        # --- too large.
        xmin = min(min(abs(self.x[:,+1])),min(abs(self.x[:,-1])))
        xmax = max(max(abs(self.x[:,+1])),max(abs(self.x[:,-1])))
        if xmax/xmin > self.maxxdisparity:
            newx = 0.5*(ave(abs(self.x[:,-1]))+ave(abs(self.x[:,+1])))
            self.x[:,-1] = - newx
            self.x[:,+1] = + newx
            for i in range(self.nparams):
                self.checklimits(i)
                self.func(self.paramsmone(i))
                self.f[i,-1] = self.loss()
                self.func(self.paramspone(i))
                self.f[i,+1] = self.loss()
    def minimize(self,niters=1,nsubiters=1,tol=1.e-10):
        self.ii = 0
        while self.ii < niters and max(self.f[:,0]) > tol:
            self.rebalancex()
            self.ii = self.ii + 1
            for i in range(self.nparams):
                self.minimize1d(i,nsubiters)
            print "Iteration number %d" %(self.ii)
            print "Current values = "
            print self.params
            print "loss = %e %e"%(self.f[-1,0],self.f[-1,0]/tol)
