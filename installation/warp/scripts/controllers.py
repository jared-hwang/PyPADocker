"""
Controller operations
=====================

These are the functions which allow installing user created functions so that
they are called at various places along the time step.

For each controller, the following three functions are defined.
 - install___: Installs a function to be called at that specified time
 - uninstall___: Uninstalls the function (so it won't be called anymore)
 - isinstalled___: Checks if the function is installed

The install can be done using a decorator, which has the prefix "callfrom". See example below.

The functions all take a function or instance method as an argument. Note that
if an instance method is used, an extra reference to the method's object is saved.

Functions can be called at the following times:
 - :py:func:`aftergenerate <installaftergenerate>`: immediately after the generate is complete
 - :py:func:`beforefs <installbeforefs>`: before the field solve
 - :py:func:`afterfs <installafterfs>`: after the field solve
 - :py:func:`beforeloadrho <installbeforeloadrho>`: before the rho is deposited, at the beginning of loadrho
 - :py:func:`afterloadrho <installafterloadrho>`: after the rho is deposited, at the end of loadrho
 - :py:func:`othereuser <installothereuser>`: during execution of electric fields gathering
 - :py:func:`beforestep <installbeforestep>`: before the time step
 - :py:func:`afterstep <installafterstep>`: after the time step
 - :py:func:`beforescraper <installbeforescraper>`: just before the particle boundary conditions are applied
 - :py:func:`particlescraper <installparticlescraper>`: just after the particle boundary conditions are applied
                    but before lost particles are processed
 - :py:func:`afterscraper <installafterscraper>`: just after the particle boundary conditions are applied
 - :py:func:`particleloader <installparticleloader>`: at the time that the standard particle loader is called
 - :py:func:`addconductor <installaddconductor>`: at the start of the multigrid solver (to installconductors)
 - :py:func:`beforeplot <installbeforeplot>`: before a plot (actually after a frame advance)
 - :py:func:`afterplot <installafterplot>`: after a plot (acutally just before a frame advance)
 - :py:func:`plseldom <installplseldom>`: during a special time step, when position and velocity are
             synchronized, specified by itplseldom or zzplseldom
 - :py:func:`plalways <installplalways>`: during a special time step, when position and velocity are
             synchronized, specified by itplalways or zzplalways
 - :py:func:`userinjection <installuserinjection>`: called when particle injection happens, after the position
                  advance and before loadrho is called, allowing a user defined
                  particle distribution to be injected each time step
 - :py:func:`userinjection2 <installuserinjection2>`: called when the 2nd part of particlesinjection happens, after the field solve,
                  allowing a user defined injected particle velocity synchronization
 - :py:func:`userparticlesinjection <installuserparticlesinjection>`: allows directly specifying the particles to be injected
 - :py:func:`userappliedfields <installuserappliedfields>`: allows directly specifying any fields to be applied to the particles
                                                            during the advance

To use a decorator, the syntax is as follows. This will install the function myplots to be called after each step.

@callfromafterstep
def myplots():
  ppzx()

This is equivalent to the following:

def myplots():
  ppzx()
installafterstep(myplots)

"""
from __future__ import generators


def controllersdoc():
    import controllers
    print controllers.__doc__

from warp import *
import types
import copy
import time

class ControllerFunction:
    """
    Class to handle the function lists.

    Note that for functions passed in that are methods of a class instance,
    a full reference of the instance is saved. This extra reference means
    that the object will not actually deleted if the user deletes the
    original reference.  This is good since the user does not need to keep
    the reference to the object (for example it can be created using a local
    variable in a function). It is also good since this allows the installed
    method to transfer across a dump and restart. It may be bad if the user
    thinks an object was deleted, but it actually isn't since it had (unkown
    to the user) installed a method in one of the controllers.

    This class also provides what is effectively a picklable function
    reference. Though it is not complete, since in some cases, functions
    won't be restorable. For example, functions typed in interactively cannot
    be restored since the source is not saved anywhere.
    """

    def __init__(self,name=None,lcallonce=0):
        self.funcs = []
        self.time = 0.
        self.timers = {}
        self.name = name
        self.lcallonce = lcallonce

    def __call__(self,*args,**kw):
        "Call all of the functions in the list"
        tt = self.callfuncsinlist(*args,**kw)
        self.time = self.time + tt
        if self.lcallonce: self.funcs = []

    def clearlist(self):
        self.funcs = []

    def __getstate__(self):
        """
        The instance is picklable. Only the names of functions are saved. A full
        reference to a method's object is saved. The names of functions replace
        the funcs attribute in the dictionary returned. Note that nothing
        special is needed on a restore since the function names will
        automatically be converted back into functions the first time they are
        called (so there is no __setstate__). For methods, the name is saved since
        instancemethods cannot be pickled.
        The ControllerFunctionContainer class below ensures that top level
        controllers are restored properly.
        """
        dict = self.__dict__.copy()
        del dict['funcs']
        funcnamelist = []
        for f in self.getpicklablefuncs():
            funcnamelist.append(f)
        dict['funcs'] = funcnamelist
        return dict

    def __nonzero__(self):
        "Returns True of functions are installed, otherwise false"
        return self.hasfuncsinstalled()

    def __len__(self):
        return len(self.funcs)

    def hasfuncsinstalled(self):
        "Checks if there are any functions installed"
        return len(self.funcs) > 0

    def _getmethodobject(self,func):
        return func[0]

    def getpicklablefuncs(self):
        """Returns the functions in a form that is picklable.
        For functions, returns the name.
        For instance methods, returns as is, the list containing the instance and method name.
        """
        for f in self.funcs:
            if isinstance(f,list):
                result = f
            elif isinstance(f,basestring):
                import __main__
                if f in __main__.__dict__:
                    result = f
                else:
                    continue
            elif isinstance(f,FunctionType):
                result = f.__name__
            else:
                result = f
            yield result

    def controllerfunclist(self):
        funclistcopy = copy.copy(self.funcs)
        for f in funclistcopy:
            if isinstance(f,list):
                object = self._getmethodobject(f)
                if object is None:
                    self.funcs.remove(f)
                    continue
                result = getattr(object,f[1])
            elif isinstance(f,basestring):
                import __main__
                if f in __main__.__dict__:
                    result = __main__.__dict__[f]
                    # --- If the function with the name is found, then replace the
                    # --- name in the list with the function.
                    self.funcs[self.funcs.index(f)] = result
                else:
                    continue
            else:
                result = f
            if not callable(result):
                print "\n\nWarning: a controller was found that is not callable."
                if self.name is not None:
                    print "For %s"%self.name
                print "Only callable objects can be installed."
                print "It is possible that the callable's name has been overwritten"
                print "by something not callable. This can happen during restart"
                print "if a function name had later been used as a variable name."
                print self.name
                if isinstance(f,basestring):
                    print "The name of the controller is ",f
                print "\n\n"
                continue
            yield result

    def installfuncinlist(self,f):
        if isinstance(f,types.MethodType):
            # --- If the function is a method of a class instance, then save a full
            # --- reference to that instance and the method name.
            finstance = f.im_self
            fname = f.__name__
            self.funcs.append([finstance,fname])
        elif callable(f):
            # --- If a function had already been installed by name, then skip the install.
            # --- This is problematic, since no warning message is given, but it is unlikely
            # --- to arise under normal circumstances.
            # --- The purpose of this check is to avoid redundant installation of functions
            # --- during a restore from a dump file. Without the check, functions that had been
            # --- installed via a decorator would be installed an extra time since the source
            # --- of the function contains the decoration (which is activated when the source
            # --- is exec'd).
            if f.__name__ not in self.funcs:
                self.funcs.append(f)
        else:
            self.funcs.append(f)

    def uninstallfuncinlist(self,f):
        # --- An element by element search is needed
        # --- f can be a function or method object, or a name (string).
        # --- Note that method objects can not be removed by name.
        funclistcopy = copy.copy(self.funcs)
        for func in funclistcopy:
            if f == func:
                self.funcs.remove(f)
                return
            elif isinstance(func,list) and isinstance(f,types.MethodType):
                object = self._getmethodobject(func)
                if f.im_self is object and f.__name__ == func[1]:
                    self.funcs.remove(func)
                    return
            elif isinstance(func,basestring):
                if f.__name__ == func:
                    self.funcs.remove(func)
                    return
            elif isinstance(f,basestring):
                if isinstance(func,basestring): funcname = func
                elif isinstance(func,list): funcname = None
                else:                        funcname = func.__name__
                if f == funcname:
                    self.funcs.remove(func)
                    return
        raise Exception('Warning: no such function had been installed')

    def isinstalledfuncinlist(self,f):
        # --- An element by element search is needed
        funclistcopy = copy.copy(self.funcs)
        for func in funclistcopy:
            if f == func:
                return 1
            elif isinstance(func,list) and isinstance(f,types.MethodType):
                object = self._getmethodobject(func)
                if f.im_self is object and f.__name__ == func[1]:
                    return 1
            elif isinstance(func,basestring):
                if f.__name__ == func:
                    return 1
        return 0

    def callfuncsinlist(self,*args,**kw):
        printentering(self.name, 2)
        bb = time.time()
        for f in self.controllerfunclist():
            #barrier()
            t1 = time.time()
            printentering(str(f), 3)
            f(*args,**kw)
            printexiting(str(f), 3)
            #barrier()
            t2 = time.time()
            # --- For the timers, use the function (or method) name as the key.
            # --- This is done since instancemethods cannot be pickled.
            self.timers[f.__name__] = self.timers.get(f.__name__,0.) + (t2 - t1)
        aa = time.time()
        printexiting(self.name, 2)
        return aa - bb

class PicklableFunction:
    """
    Class allowing pickable functions

    This class provides a picklable function reference. Though it is not
    complete, since in some cases, functions won't be restorable. For example,
    functions typed in interactively cannot be restored since the source is not
    saved anywhere.

    Note that for functions passed in that are methods of a class instance,
    a full reference of the instance is saved. This extra reference means
    that the object will not actually deleted if the user deletes the
    original reference.  This is good since the user does not need to keep
    the reference to the object (for example it can be created using a local
    variable in a function). It is also good since this allows the installed
    method to transfer across a dump and restart. It may be bad if the user
    thinks an object was deleted, but it actually isn't since it had (unkown
    to the user) installed a method in one of the controllers.

    """

    def __init__(self,func):
        self.installfunc(func)
        self.time = 0.

    def __call__(self,*args,**kw):
        "Call the function"
        return self.callfunc(*args,**kw)

    def __getstate__(self):
        """
        The instance is picklable. Only the names of functions are saved. A full
        reference to a method's object is saved. The names of function replace
        the func attribute in the dictionary returned. Note that nothing
        special is needed on a restore since the function names will
        automatically be converted back into functions the first time they are
        called (so there is no __setstate__). For methods, the name is saved since
        instancemethods cannot be pickled.
        The ControllerFunctionContainer class below ensures that top level
        controllers are restored properly.
        """
        dict = self.__dict__.copy()
        dict['func'] = self.getpicklablefunc()
        return dict

    def _getmethodobject(self,func):
        return func[0]

    def getpicklablefunc(self):
        """Returns the function in a form that is picklable.
        For functions, returns the name.
        For instance methods, returns as is, the list containing the instance and method name.
        """
        if isinstance(self.func,list):
            result = self.func
        elif isinstance(self.func,basestring):
            import __main__
            if self.func in __main__.__dict__:
                result = self.func
            else:
                # --- The function couldn't be found
                result = None
        elif isinstance(self.func,FunctionType):
            result = self.func.__name__
        else:
            result = self.func
        return result

    def callablefunc(self):
        if isinstance(self.func,list):
            object = self._getmethodobject(self.func)
            if object is None:
                raise Exception("The method's class instance could not be found")
            result = getattr(object,self.func[1])
        elif isinstance(self.func,basestring):
            import __main__
            if self.func in __main__.__dict__:
                result = __main__.__dict__[self.func]
                # --- If the function with the name is found, then replace the
                # --- name in the list with the function.
                self.func = result
            else:
                raise Exception("The function couldn't be found")
        else:
            result = self.func
        if not callable(result):
            print "\n\nWarning: a controller was found that is not callable."
            print "Only callable objects can be installed."
            print "It is possible that the callable's name has been overwritten"
            print "by something not callable. This can happen during restart"
            print "if a function name had later been used as a variable name."
            if isinstance(self.func,basestring):
                print "The name of the controller is ",self.func
            print "\n\n"
            raise Exception("The function is not callable")
        return result

    def installfunc(self,f):
        if isinstance(f,types.MethodType):
            # --- If the function is a method of a class instance, then save a full
            # --- reference to that instance and the method name.
            finstance = f.im_self
            fname = f.__name__
            self.func = [finstance,fname]
        elif callable(f):
            self.func = f
        else:
            self.func = f

    def isfunc(self,f):
        if f == self.func:
            return 1
        elif isinstance(self.func,list) and isinstance(f,types.MethodType):
            object = self._getmethodobject(self.func)
            if f.im_self is object and f.__name__ == self.func[1]:
                return 1
        elif isinstance(self.func,basestring):
            if f.__name__ == self.func:
                return 1
        return 0

    def callfunc(self,*args,**kw):
        t1 = time.time()
        f = self.callablefunc()
        result = f(*args,**kw)
        t2 = time.time()
        self.time += t2 - t1
        return result

#=============================================================================

# --- Now create the actual instances.
aftergenerate = ControllerFunction('aftergenerate')
beforefs = ControllerFunction('beforefs')
afterfs = ControllerFunction('afterfs')
beforeloadrho = ControllerFunction('beforeloadrho')
afterloadrho = ControllerFunction('afterloadrho')
othereuser = ControllerFunction('othereuser')
beforescraper = ControllerFunction('beforescraper')
afterscraper = ControllerFunction('afterscraper')
callscraper = ControllerFunction('callscraper')
callparticleloader = ControllerFunction('callparticleloader')
calladdconductor = ControllerFunction('calladdconductor')
callbeforestepfuncs = ControllerFunction('callbeforestepfuncs')
callafterstepfuncs = ControllerFunction('callafterstepfuncs')
callbeforeplotfuncs = ControllerFunction('callbeforeplotfuncs')
callafterplotfuncs = ControllerFunction('callafterplotfuncs')
callplseldomfuncs = ControllerFunction('callplseldomfuncs')
callplalwaysfuncs = ControllerFunction('callplalwaysfuncs')
callafterrestartfuncs = ControllerFunction('callafterrestartfuncs',lcallonce=1)
userinjection = ControllerFunction('userinjection')
userinjection2 = ControllerFunction('userinjection2')
generateuserparticlesforinjection = ControllerFunction('generateuserparticlesforinjection')
userappliedfields = ControllerFunction('userappliedfields')

#=============================================================================
class ControllerFunctionContainer:
    """
  This is a somewhat kludgy fix to how to get any saved functions restored.
  A single instance of this class is created and this instance is what is saved
  in a dump. This instance will have a list of the controllers, so the
  controllers will be saved, but not as top level python variables.
  Upon restoration, this container will go through each of the saved controllers
  and reinstall the functions saved therein. This installs the functions in the
  original set of controllers created when this module was first imported.
  Anything that may have already been installed will by default be deleted.
    """
    # --- This needs to be a class attribute, so that it can be set before
    # --- and instance is being restored and __setstate__ called.
    clearfunctionlists = True
    def __init__(self,clist):
        self.clist = clist
    def __setstate__(self,dict):
        import controllers
        import __main__
        self.__dict__.update(dict)
        for c in self.clist:
            origcontroller = controllers.__dict__[c.name]
            # --- NOTE: If clearfunctionlists was always True, the rest of the clist loop could
            # --- be simplified to one of these single statements:
            # --- origcontroller.funcs = c.funcs
            # --- OR
            # --- controllers.__dict__[c.name] = c
            # --- With this second version, the resetting of self.clist below would not be needed.
            if ControllerFunctionContainer.clearfunctionlists:
                origcontroller.clearlist()
            for f in c.funcs:
                if isinstance(f,basestring):
                    # --- Check if f is already in the original list of functions,
                    # --- and skip it if it is. Both the function name (f) and the
                    # --- actual function in main are checked.
                    # --- This will be the case if, for example, the user execs the
                    # --- original input file, which sets up some functions, before
                    # --- doing the restart.
                    origfuncs = origcontroller.funcs
                    try:
                        ffunc = __main__.__dict__[f]
                    except KeyError:
                        ffunc = None
                    if (f not in origfuncs and ffunc not in origfuncs):
                        origcontroller.installfuncinlist(f)
                else:
                    # --- Otherwise, f is a method, so it can be directly installed.
                    # --- A check is still made to ensure it isn't installed twice.
                    # --- The check is only needed temporarily until the classes
                    # --- are fixed to not resinstall in the getstate.
                    ffunc = getattr(f[0],f[1])
                    if not origcontroller.isinstalledfuncinlist(ffunc):
                        origcontroller.installfuncinlist(ffunc)
        # --- The clist is obtained from the original instance in the controllers
        # --- module so that the list contains references to the original
        # --- controller instances. This is needed, since in the next dump,
        # --- this instance will be written out and must contain an updated
        # --- list of controllers.
        self.clist = controllers.__dict__['controllerfunctioncontainer'].clist

    def printtimers(self,tmin=1.,lminmax=0.,ff=None):
        """Prints timings of install functions.
     - tmin=1.: only functions with time greater than tmin will be printed
        """
        if ff is None: ff = sys.stdout
        for c in self.clist:
            for fname,time in c.timers.items():
                vlist = array(gather(time))
                if me > 0: continue
                vsum = sum(vlist)
                if vsum <= tmin: continue
                vrms = sqrt(max(0.,ave(vlist**2) - ave(vlist)**2))
                ff.write('%20s %s %10.4f  %10.4f %10.4f'%(c.name,fname,vsum,vsum/npes,vrms))
                if lminmax:
                    vmin = min(vlist)
                    vmax = max(vlist)
                    ff.write('  %10.4f  %10.4f'%(vmin,vmax))
                if top.it > 0:
                    ff.write('   %10.4f'%(vsum/npes/(top.it)))
                ff.write('\n')

# --- This is primarily needed by warp.py so that these objects can be removed
# --- from the list of python objects which are not written out.
controllerfunctioncontainer = ControllerFunctionContainer(
                               [aftergenerate,beforefs,afterfs,
                                beforeloadrho,afterloadrho,othereuser,
                                beforescraper,afterscraper,callscraper,
                                callparticleloader,calladdconductor,
                                callbeforestepfuncs,callafterstepfuncs,
                                callbeforeplotfuncs,callafterplotfuncs,
                                callplseldomfuncs,callplalwaysfuncs,
                                callafterrestartfuncs,
                                userinjection,userinjection2,
                                generateuserparticlesforinjection,
                                userappliedfields])


#=============================================================================
# ----------------------------------------------------------------------------
def callfromaftergenerate(f):
    installaftergenerate(f)
    return f
def installaftergenerate(f):
    "Adds a function to the list of functions called after a generate"
    aftergenerate.installfuncinlist(f)
def uninstallaftergenerate(f):
    "Removes the function from the list of functions called after a generate"
    aftergenerate.uninstallfuncinlist(f)
def isinstalledaftergenerate(f):
    "Checks if the function is called after a generate"
    return aftergenerate.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfrombeforefs(f):
    installbeforefs(f)
    return f
def installbeforefs(f):
    "Adds a function to the list of functions called before a field-solve"
    beforefs.installfuncinlist(f)
    w3d.lbeforefs = true
def uninstallbeforefs(f):
    "Removes the function from the list of functions called before a field-solve"
    beforefs.uninstallfuncinlist(f)
    if not beforefs.hasfuncsinstalled(): w3d.lbeforefs = false
def isinstalledbeforefs(f):
    "Checks if the function is called before a field-solve"
    return beforefs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterfs(f):
    installafterfs(f)
    return f
def installafterfs(f):
    "Adds a function to the list of functions called after a field-solve"
    afterfs.installfuncinlist(f)
    w3d.lafterfs = true
def uninstallafterfs(f):
    "Removes the function from the list of functions called after a field-solve"
    afterfs.uninstallfuncinlist(f)
    if not afterfs.hasfuncsinstalled(): w3d.lafterfs = false
def isinstalledafterfs(f):
    "Checks if the function is called after a field-solve"
    return afterfs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfrombeforeloadrho(f):
    installbeforeloadrho(f)
    return f
def installbeforeloadrho(f):
    "Adds a function to the list of functions called before a load rho"
    beforeloadrho.installfuncinlist(f)
    w3d.lbeforelr = true
def uninstallbeforeloadrho(f):
    "Removes the function from the list of functions called before a load rho"
    beforeloadrho.uninstallfuncinlist(f)
    if not beforeloadrho.hasfuncsinstalled(): w3d.lbeforelr = false
def isinstalledbeforeloadrho(f):
    "Checks if the function is called before a load rho"
    return beforeloadrho.isinstalledfuncinlist(f)

# --- This are defined for backwards compatibility
installbeforelr = installbeforeloadrho
uninstallbeforelr = uninstallbeforeloadrho
isinstalledbeforelr = isinstalledbeforeloadrho

# ----------------------------------------------------------------------------
def callfromafterloadrho(f):
    installafterloadrho(f)
    return f
def installafterloadrho(f):
    "Adds a function to the list of functions called after a load rho"
    afterloadrho.installfuncinlist(f)
    w3d.lafterloadrho = true
def uninstallafterloadrho(f):
    "Removes the function from the list of functions called after a load rho"
    afterloadrho.uninstallfuncinlist(f)
    if not afterloadrho.hasfuncsinstalled(): w3d.lafterloadrho = false
def isinstalledafterloadrho(f):
    "Checks if the function is called after a load rho"
    return afterloadrho.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromothereuser(f):
    installothereuser(f)
    return f
def installothereuser(f):
    "Adds a function to the list of functions called during the electric fields gathering"
    othereuser.installfuncinlist(f)
    w3d.lothereuser = true
def uninstallothereuser(f):
    "Removes the function from the list of functions called during the electric fields gathering"
    othereuser.uninstallfuncinlist(f)
    if not othereuser.hasfuncsinstalled(): w3d.lothereuser = false
def isinstalledothereuser(f):
    "Checks if the function is called during the electric fields gathering"
    return othereuser.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfrombeforescraper(f):
    installbeforescraper(f)
    return f
def installbeforescraper(f):
    "Adds a function to the list of functions called before scraping particles"
    beforescraper.installfuncinlist(f)
    w3d.lbeforescraper = true
def uninstallbeforescraper(f):
    "Removes the function from the list of functions called before scraping particles"
    beforescraper.uninstallfuncinlist(f)
    if not beforescraper.hasfuncsinstalled(): w3d.lbeforescraper = false
def isinstalledbeforescraper(f):
    "Checks if the function is called before scraping particles"
    return beforescraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterscraper(f):
    installafterscraper(f)
    return f
def installafterscraper(f):
    "Adds a function to the list of functions called after scraping particles"
    afterscraper.installfuncinlist(f)
    w3d.lafterscraper = true
def uninstallafterscraper(f):
    "Removes the function from the list of functions called after scraping particles"
    afterscraper.uninstallfuncinlist(f)
    if not afterscraper.hasfuncsinstalled(): w3d.lafterscraper = false
def isinstalledafterscraper(f):
    "Checks if the function is called after scraping particles"
    return afterscraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromparticlescraper(f):
    installparticlescraper(f)
    return f
def installparticlescraper(f):
    "Adds a function to the list of functions called to scrape particles"
    callscraper.installfuncinlist(f)
    w3d.lcallscraper = true
def uninstallparticlescraper(f):
    "Removes the function from the list of functions called to scrape particles"
    callscraper.uninstallfuncinlist(f)
    if not callscraper.hasfuncsinstalled(): w3d.lcallscraper = false
def isinstalledparticlescraper(f):
    "Checks if the function is called to scrape particles"
    return callscraper.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromparticleloader(f):
    installparticleloader(f)
    return f
def installparticleloader(f):
    "Adds a function to the list of functions called to load particles"
    callparticleloader.installfuncinlist(f)
    w3d.lcallparticleloader = true
def uninstallparticleloader(f):
    "Removes the function from the list of functions called to load particles"
    callparticleloader.uninstallfuncinlist(f)
    if not callparticleloader.hasfuncsinstalled():
        w3d.lcallparticleloader = false
def isinstalledparticleloader(f):
    "Checks if the function is called to load particles"
    return callparticleloader.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromaddconductor(f):
    installaddconductor(f)
    return f
def installaddconductor(f):
    "Adds a function to the list of functions called to add conductors"
    calladdconductor.installfuncinlist(f)
    f3d.laddconductor = true
def uninstalladdconductor(f):
    "Removes the function from the list of functions called to add conductors"
    calladdconductor.uninstallfuncinlist(f)
    if not calladdconductor.hasfuncsinstalled(): f3d.laddconductor = false
def isinstalledaddconductor(f):
    "Checks if the function is called to add conductors"
    return calladdconductor.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfrombeforestep(f):
    installbeforestep(f)
    return f
def installbeforestep(f):
    "Adds a function to the list of functions called before a step"
    callbeforestepfuncs.installfuncinlist(f)
def uninstallbeforestep(f):
    "Removes the function from the list of functions called before a step"
    callbeforestepfuncs.uninstallfuncinlist(f)
def isinstalledbeforestep(f):
    "Checks if the function is called before a step"
    return callbeforestepfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterstep(f):
    installafterstep(f)
    return f
def installafterstep(f):
    "Adds a function to the list of functions called after a step"
    callafterstepfuncs.installfuncinlist(f)
def uninstallafterstep(f):
    "Removes the function from the list of functions called after a step"
    callafterstepfuncs.uninstallfuncinlist(f)
def isinstalledafterstep(f):
    "Checks if the function is called after a step"
    return callafterstepfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfrombeforeplot(f):
    installbeforeplot(f)
    return f
def installbeforeplot(f):
    "Adds a function to the list of functions called before a plot"
    beforeplotfuncs.installfuncinlist(f)
def uninstallbeforeplot(f):
    "Removes the function from the list of functions called before a plot"
    beforeplotfuncs.uninstallfuncinlist(f)
def isinstalledbeforeplot(f):
    "Checks if the function is called before a plot"
    return beforeplotfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterplot(f):
    installafterplot(f)
    return f
def installafterplot(f):
    "Adds a function to the list of functions called after a plot"
    callafterplotfuncs.installfuncinlist(f)
def uninstallafterplot(f):
    "Removes the function from the list of functions called after a plot"
    callafterplotfuncs.uninstallfuncinlist(f)
def isinstalledafterplot(f):
    "Checks if the function is called after a plot"
    return callafterplotfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromplseldom(f):
    installplseldom(f)
    return f
def installplseldom(f):
    "Adds a function to the list of functions controlled by itplseldom and zzplseldom"
    callplseldomfuncs.installfuncinlist(f)
def uninstallplseldom(f):
    "Removes the function from the list of functions controlled by itplseldom and zzplseldom"
    callplseldomfuncs.uninstallfuncinlist(f)
def isinstalledplseldom(f):
    "Checks if the function is controlled by itplseldom and zzplseldom"
    return callplseldomfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromplalways(f):
    installplalways(f)
    return f
def installplalways(f):
    "Adds a function to the list of functions controlled by itplalways and zzplalways"
    callplalwaysfuncs.installfuncinlist(f)
def uninstallplalways(f):
    "Removes the function from the list of functions controlled by itplalways and zzplalways"
    callplalwaysfuncs.uninstallfuncinlist(f)
def isinstalledplalways(f):
    "Checks if the function is controlled by itplalways and zzplalways"
    return callplalwaysfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromafterrestart(f):
    installafterrestart(f)
    return f
def installafterrestart(f):
    "Adds a function to the list of functions called immediately after a restart"
    callafterrestartfuncs.installfuncinlist(f)
def uninstallafterrestart(f):
    "Removes the function from the list of functions called immediately after a restart"
    callafterrestartfuncs.uninstallfuncinlist(f)
def isinstalledafterrestart(f):
    "Checks if the function is called immediately after a restart"
    return callafterrestartfuncs.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromuserinjection(f):
    installuserinjection(f)
    return f
def installuserinjection(f):
    """
  Adds a user defined function that is to be called when particle
  injection happens, after the position advance and before loadrho is
  called, allowing a user defined particle distribution to be injected
  each time step"""
    userinjection.installfuncinlist(f)
def uninstalluserinjection(f):
    "Removes the function installed by installuserinjection"
    userinjection.uninstallfuncinlist(f)
def isinstalleduserinjection(f):
    "Checks if the function is called when particles injection happens"
    return userinjection.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromuserinjection2(f):
    installuserinjection2(f)
    return f
def installuserinjection2(f):
    """
  Adds a user defined function that is to be called when the second part of
  particle injection happens, after the field solve,
  allowing a user defined injected particle velocity synchronization"""
    userinjection2.installfuncinlist(f)
def uninstalluserinjection2(f):
    "Removes the function installed by installuserinjection2"
    userinjection2.uninstallfuncinlist(f)
def isinstalleduserinjection2(f):
    "Checks if the function is called when particles injection happens"
    return userinjection2.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromuserparticlesinjection(f):
    installuserparticlesinjection(f)
    return f
def installuserparticlesinjection(f):
    """
  Adds a user defined function that is to be called during injection which
  allows the user to specify the distribution of particles on the emitting
  surface. For expert use only.
  To use, the installed function should set w3d.npgrp to the number of
  particles to inject, call gchange("Setpwork3d") to allocate the arrays, and
  fill the arrays w3d.xt, yt, uxt, uyt, and uzt with the particle data. The
  particles start on the emitting surface and the code will advance them away
  from the surface. The function will be called once for each species each time
  step, with the variable w3d.inj_js set to the species being injected. Note
  that if no particles are to be injected, set w3d.npgrp=0 to avoid injection
  of bad particles."""
    w3d.l_inj_user_particles = true
    generateuserparticlesforinjection.installfuncinlist(f)
def uninstalluserparticlesinjection(f):
    "Removes the function installed by installuserparticlesinjection"
    generateuserparticlesforinjection.uninstallfuncinlist(f)
    if not generateuserparticlesforinjection.hasfuncsinstalled():
        w3d.l_inj_user_particles = false
def isinstalleduserparticlesinjection(f):
    "Checks if the function is called during injection"
    return generateuserparticlesforinjection.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
def callfromuserappliedfields(f):
    installuserappliedfields(f)
    return f
def installuserappliedfields(f):
    """
  Adds a user defined function which can specify E and B fields which are applied
  to the particles during the particle advance. This function is called at the
  start of padvnc3d. Note, that the following must be done in this install function:
   - lresetparticlee and/or lresetparticleb must be set to false
   - the E and/or B of the particles must be zeroed if they are otherwise unset
  """
    userappliedfields.installfuncinlist(f)
    w3d.luserappliedfields = true
def uninstalluserappliedfields(f):
    "Removes the function installed by installuserappliedfields"
    userappliedfields.uninstallfuncinlist(f)
    if not userappliedfields.hasfuncsinstalled():
        w3d.luserappliedfields = false
def isinstalleduserappliedfields(f):
    "Checks if the function is called when which applies fields"
    return userappliedfields.isinstalledfuncinlist(f)

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def fixcontrollersfromolddump():
    import __main__
    import controllers
    controllernames = ['aftergenerate','beforefs','afterfs','callscraper',
                       'calladdconductor','callbeforestepfuncs',
                       'callafterstepfuncs','callbeforeplotfuncs',
                       'callafterplotfuncs','callplseldomfuncs',
                       'callplalwaysfuncs']
    for cname in controllernames:
        if cname in __main__.__dict__:
            controller = __main__.__dict__[cname]
            if 'funcnamelist' in controller.__dict__:
                controllers.__dict__[controller.name] = controller
                controller.funcs = controller.funcnamelist
                del controller.funcnamelist
        else:
            print "Controller ",cname," not found"
