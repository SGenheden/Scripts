# Author: Samuel Genheden, 2012-2014


"""
Module containing classes for calculating quality metrics
"""


import itertools as it
import sys
import inspect

import numpy as np
import numpy.random as rnd

import scipy.stats.mstats as stat
import scipy.stats.stats as stats


class Qmetric :
  """
    Qmetric is the parent class for all quality metric classes

    Each class that inherit this class should implement the following methods
      eval     - should return the quality metric 
      h0       - should return the null hypothesis (optional)
      init     - should initialize additional member variables (optional)
      lower    - should return the lower limit (optional)
      name     - should return a string name of the quality metric
      upper    - should return the upper limit (optional)
    and the following class variable
      tag      - a short string on which the user can select the quality metric
  """

  def __init__(self,pred,prederr,expect,expecterr,**kwargs) : 
    """
    
    Attributes
    ----------
    pred : numpy array
      predicted values
    prederr : numpy array
      predicted errors
    expect : numpy array
      expected values
    expecterr : numpy array
      expected errors
    kwargs : dictionary
      dictionary of global variables
    """
    self.std = 0.0
    self.av = 0.0
    self.bias = 0.0
    self.unbiased = 0.0
    self.median = 0.0
    self.nlow = 0.0
    self.nhigh = 0.0
    self.dlow = 0.0
    self.dhigh = 0.0
    self.bootstrapped = []
    self.delta = []
    if "verbose" in kwargs :
      self.verbose = kwargs["verbose"]
    else :
      self.verbose = False
    if "splitstat" in kwargs :
      self.splitstat = kwargs["splitstat"]
    else :
      self.splitstat = True
    self.init(pred,prederr,expect,expecterr,**kwargs)
    self.biased = self.eval(pred,expect)

  def bootstrap(self,pred,expect) :
    """
    Calculate bootstrapped values
    
    Parameters
    ----------
    pred : numpy array
      the bootstrapped predicted values
    expect : numpy array
      the bootstrapped expected values
    """
    nboots = pred.shape[1]
    nval = pred.shape[0]
    if nboots < 1 : return
    self.bootstrapped = np.zeros(nboots)
    for i in range(nboots) :
      self.bootstrapped[i] = self.eval(pred[:,i],expect[:,i])
    self.delta = self.bootstrapped - self.biased
    self.std = np.std(self.bootstrapped,ddof=1)
    self.av = np.sum(self.bootstrapped)/nboots
    self.bias = np.sum(self.delta)/nboots
    self.unbiased = self.biased+self.bias
    self.median = np.median(self.bootstrapped)
    self.nlow = self.lower(self.biased - 1.96*self.std/np.sqrt(nval))
    self.nhigh = self.upper(self.biased + 1.96*self.std/np.sqrt(nval))
    self.dlow = self.lower(self.unbiased - stat.mquantiles(self.delta,prob=[0.95]))
    self.dhigh = self.upper(self.unbiased - stat.mquantiles(self.delta,prob=[0.05]))

  def eval(self,pred,expect) :
    """
    Evaluates this quality metric
    
    Parameters
    ----------
    pred : numpy array
      the predicted values
    expect : numpy array
      the expected values
      
    Returns
    -------
    float
      the quality metric
    """
    pass

  def h0(self,expect) :
    """
    Evaluates the null hypothesis
    
    Parameters
    ----------
    expect : numpy array
      the expected values
      
    Returns
    -------
    float or None
      the null hypothesis or if it is non existent
    """
    return None

  def init(self,pred,prederr,expect,expecterr,**kwargs) :
    """
    Initialise the object
    
    Attributes
    ----------
    pred : numpy array
      predicted values
    prederr : numpy array
      predicted errors
    expect : numpy array
      expected values
    expecterr : numpy array
      expected errors
    kwargs : dictionary
      dictionary of global variables
    """
    pass

  def lower(self,value) :
    """
    Returns the value it self or an lower bound of the quality metric
    """
    return value

  def name(self) :
    """
    Returns a short string abbreviation of the quality metric
    """
    pass 

  def plotDelta(self) :
    """
    Plot a histogram of the delta distribution
    """
    if len(self.delta) == 0 : return
    import matplotlib.pyplot as plt
    plt.clf()
    plt.hist(self.delta[:,idx],bins=75,normed=True)  
    plt.savefig(self.name()+"_delta.png",format="png")
 
  def plotSample(self) :
    """
    Plot a histogram of the sample distribution
    """
    if len(self.bootstrapped) == 0 : return
    import matplotlib.pyplot as plt
    plt.clf()
    plt.hist(self.bootstrapped,bins=75,normed=True)
    plt.savefig(self.name()+"_samples.png",format="png") 

  def __str__(self) :
    if self.verbose :
      return "%-10s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f"%(self.name(),self.biased,self.std,self.av,self.nlow,self.nhigh,self.bias,self.unbiased,self.median,self.dlow,self.dhigh)
    else :
      return "%-10s %9.4f %9.4f"%(self.name(),self.biased,self.std)
  """
    Returns the value it self or an upper bound of the quality metric
    """ 
  def upper(self,value) :
    return value


class QualityCollection :
  """
  Class for collection of Qmetric objects
  """
  
  def __init__(self,pred,prederr,expect,expecterr,metrics,**kwargs) :
    """
    
    """
    self.metriclst = []
    # Create all the Qmetric objects based on a key string and put them in a list
    QualityCollection._init_metrics()
    for str in metrics :
      try :
        q = QualityCollection._metrics[str](pred,prederr,expect,expecterr,**kwargs)
      except KeyError :
        print "Warning: Could not identify metric (%s)"%str
        print "Available metric keys are: ",QualityCollection.metrickeys()
      except Exception as e:
        print e.args[0]
        print "This metric (%s) will be ignored"%str
      else :
        self.metriclst.append(q)
    if len(self.metriclst) == 0 :
      raise Exception("Error: Trying to create a collection with zero metrics.")
    self.pred = pred
    self.prederr = prederr
    self.expect = expect
    self.expecterr = expecterr
    if "verbose" in kwargs :
      self.verbose = kwargs["verbose"]
    else :
      self.verbose = False
    if "nboots" in kwargs :
      self.nboots = kwargs["nboots"]
    else :
      self.nboots = 0
    rnd.seed(rnd.randint(107575))
 
  def bootstrap(self) :
    """
    Created bootstrapped predictions and expected values
    and calculate statistics for each of the Qmetric objects
    """
    if self.nboots < 1 : return
    nval = self.pred.shape[0]
    # Create bootstrapped values of the predicted values
    pboots = rnd.randn(nval,self.nboots)
    for i in range(self.nboots) :
      pboots[:,i] = np.multiply(pboots[:,i],self.prederr)+self.pred
    # Create bootstrapped values of the expected values or just copy them
    if len(self.expecterr) == nval :
      eboots = rnd.randn(nval,self.nboots)
      for i in range(self.nboots) :
        eboots[:,i] = np.multiply(eboots[:,i],self.expecterr)+self.expect
    else :
      eboots = np.zeros((nval,self.nboots))
      for i in range(self.nboots) :
        eboots[:,i] = eboots[:,i]+self.expect
    # Calculate bootstrapped quality metrics
    for m in self.metriclst :
      m.bootstrap(pboots,eboots)
      #m.plotBoot()

  def plotDistributions(self,plotDelta=True,plotSample=True) :
    """
    Plot the delta and sample distribution for the quality metrics
    
    Parameters
    ----------
    plotDelta : boolean, optional
      if to plot the delta distribution
    plotSample : boolean, optional
      if to plot the sample distribution
    """
    for m in self.metriclst :
      if plotDelta : m.plotDelta()
      if plotSample : m.plotSample()

  def printH0(self) :
    """
    Print the H0 for each of the Qmetric objects
    """
    for m in self.metriclst :
      h0 = m.h0(self.expect)
      if h0 != None : 
        print "\nNull hypothesis of %s = %9.4f"%(m.name(),h0)
 
  def printStat(self) :
    """
    Print statistics for each of the Qmetric objects
    """
    print ""
    if self.verbose :
      print "---------------------------------------------------------------------------------------------------------------"
      print "              Biased     Stdev      Mean      95%-confidence     Bias   Unbiased    Median  95%-confidence int."
      print "---------------------------------------------------------------------------------------------------------------"
    else :
      print "------------------------------"
      print "              Biased    Stdev"
      print "------------------------------"
    for m in self.metriclst :
      print m
     
  # Class variables and methods
     
  _metrickeys = []
  _metrics = {}
   
  def metrickeys(cls) :
    cls._init_metrics()
    return cls._metrickeys
    
  def metrics(cls) :
    cls._init_metrics()
    return cls._metrics[tag]
    
  def _init_metrics(cls) :
    """
    Initialises a dictionary of Qmetric sub classes and returns
    such a class based on a tag
    """
    if len(cls._metrickeys) == 0 :
      def pred(c) :
        return inspect.isclass(c) and c.__module__ == cls.__module__ and issubclass(c,Qmetric) and c.__name__ != "Qmetric"
      for c in  inspect.getmembers(sys.modules[__name__],pred) :
        cls._metrickeys.append(c[1].tag)
        cls._metrics[c[1].tag] = c[1]
        
  metrickeys = classmethod(metrickeys)      
  metrics = classmethod(metrics)
  _init_metrics = classmethod(_init_metrics)
    
#
# ----------------------------------------------------------------------------------------
# Below here is all classes that inherits Qmetric
# ----------------------------------------------------------------------------------------
#

class _Qauc(Qmetric) :
  """
  Class to calculate the area under the curve for an ROC
  """  
  
  tag = "auc"
  
  def eval(self,pred,expect) :
    maxv = max(pred)
    minv = min(pred)
    maxv = pred(1)
    dt = (maxv-minv) / 20.0
    l = len(pred)
    roc = []
    for i in range(21) :
      t = minv+i*dt
      tp = 0.0
      fp = 0.0
      fn = 0.0
      tn = 0.0
      for j in range(l) :
        if pred[j] < t and self.active[j] :
          tp = tp + 1.0
        elif pred[j] < t and not self.active[j] :
          fp = fp + 1.0
        elif pred[j] >= t and self.active[j] :
          fn = fn + 1.0
        elif pred[j] >= t and not self.active[j] :
          tn = tn + 1.0
      roc.append((fp/(fp+tn),tp/(tp+fn)))
    auc = 0.0
    for i in range(1,21) :
      auc = auc + (roc[i][0]-roc[i-1][0])*(roc[i][1]+roc[i-1][1])
    return auc*0.50

  def init(self,pred,prederr,expect,expecterr,**kwargs) :
    # Determine which items are active and which are inactive
    # This is stored in the list self.active 
    nmax = 0
    l = len(pred)
    for i in range(l) :
      npair = 0
      for j in range(i+1,l) :
        if expect[i] == expect[j] : npair = npair + 1
      if npair > nmax :
        nmax = npair
        nullval = expect[i]
    if npair == 0 :
      raise Exception("Error: Could not identify active items in Qauc, npair=0")
    self.active = []
    for e in expect :
      self.active(not(e==nullval))

  def lower(self,value) :
    return max(value,0.0)
  
  name = "AUC"
    
  def upper(self,value) :
    return min(value,1.0)

class _Qabsmedian(Qmetric) :
  """
  Class to calculate the absolute median of deviation
  """
  
  tag = "absmed"

  def eval(self,pred,expect) :
    return np.median(np.abs(pred-expect))
    
  def name(self) :
    return "Abs. med."

class _Qintercept(Qmetric) :
  """
  Class to calculate the intercept of a linear regression
  """
  
  tag = "inter"

  def eval(self,pred,expect) :
    z = np.polyfit(expect,pred,1)
    return z[1]
  
  def name(self) :
    return "intercept"

 
class _Qmad(Qmetric) :
  """
  Class to calculate the mean absolute deviation (or mean unsigned error)
  """
  
  tag = "mad"

  def eval(self,pred,expect) :
    return np.sum(np.abs(pred-expect)) / len(pred)
  
  def name(self) :
    return "MAD"

class _Qmadtr(Qmetric) :
  """
  Class to calculate the translated mean absolute deviation
  """
  
  tag = "madtr"

  def eval(self,pred,expect) :
    mse = np.sum(pred-expect) / len(pred)
    return np.sum(np.abs(pred-expect-mse)) / len(pred)
  
  def h0(self,expect) :
    av = np.sum(expect)/len(expect)
    return np.sum(np.abs(expect-av)) / len(expect)
  
  def name(self) :
    return "MADtr"

class _Qmedian(Qmetric) :
  """
  Class to calculate the median of deviation
  """
  
  tag = "med"

  def eval(self,pred,expect) :
    return np.median(pred-expect)
  
  def name(self) :
    return "Median"

class _Qmsd(Qmetric) :
  """
  Class to calculate the mean signed deviation
  """
  
  tag = "msd"

  def eval(self,pred,expect) :
    return np.sum(pred-expect) / len(pred)
  
  def name(self) :
    return "MSD"

class _Qmq(Qmetric) :
  """
  Class to calculate the mean quote
  """
  
  tag = "mq"

  def eval(self,pred,expect) :
    return np.sum(np.divide(pred,expect)) / len(pred)
  
  def name(self) :
    return "MQ"

class _Qpi(Qmetric) :
  """
  Class to calculate the predictive index
  """
  
  tag = "pi"

  def eval(self,pred,expect) :
    wsum = 0.0
    csum = 0.0
    r = range(len(pred))
    for i,j in it.product(r,r) :
      a = expect[j] - expect[i]
      b = pred[j] - pred[i]
      if abs(b) > 0.0 :
        if a/b < 0.0 :
          csum = csum - abs(a)
        else :
          csum = csum + abs(a)
      wsum = wsum + abs(a)
    return csum / wsum
    
  def lower(self,value) :
    return max(value,0.0)
  
  def name(self) :
    return "PI"

  def upper(self,value) :
    return min(value,1.0)

class _Qr(Qmetric) :
  """
  Class to calculate the linear correlation coefficient
  """
  
  tag = "r"

  def eval(self,pred,expect) :
    r = np.corrcoef(pred,expect)
    return r[1,0]
    
  def lower(self,value) :
    return max(value,-1.0)
  
  def name(self) :
    return "r"
    
  def upper(self,value) :
    return min(value,1.0)

class _Qr2(Qmetric) :
  """
  Class to calculate the square of the linear correlation coefficient
  """
  
  tag = "r2"

  def eval(self,pred,expect) :
    r = np.corrcoef(pred,expect)
    return r[1,0]**2
    
  def lower(self,value) :
    return max(value,0.0)
  
  def name(self) :
    return "r2"
    
  def upper(self,value) :
    return min(value,1.0)

class _Qregmad(Qmetric) :
  """
  Class to calculate the mean absolute deviation to a linear regression
  """
  
  tag = "regmad"

  def eval(self,pred,expect) :
    z = np.polyfit(pred,expect,1)
    return np.sum(np.abs(expect-(pred*z[0]+z[1]))) / len(pred)
  
  def name(self) :
    return "Reg. MAD"

class _Qrho(Qmetric) :
  """
  Class to calculate Spearman's rank correlation
  """
  
  tag = "rho"

  def eval(self,pred,expect) :
    rho = stats.spearmanr(pred, expect)
    return rho[0]
    
  def lower(self,value) :
    return max(value,-1.0)
  
  def name(self) :
    return "rho"
    
  def upper(self,value) :
    return min(value,1.0)

class _Qrmsd(Qmetric) :
  """
  Class to calculate the root mean squared deviation
  """
  
  tag = "rmsd"

  def eval(self,pred,expect) :
    return np.sqrt(np.sum((pred-expect)**2) / len(pred))
  
  def name(self) :
    return "RMSD"

class _Qq(Qmetric) :
  """
  Class to calculate Q-value
  """
  
  tag = "q"

  def eval(self,pred,expect) :
    sum1 = np.sum((pred-expect)**2)
    sum2 = np.sum(expect**2)
    return sum1 / sum2
  
  def name(self) :
    return "Q"

class _Qslope(Qmetric) :
  """
  Class to calculate the slope of a linear regression
  """
  
  tag = "slope"

  def eval(self,pred,expect) :
    z = np.polyfit(pred,expect,1)
    return z[0]
  
  def name(self) :
    return "slope"

class _Qtau(Qmetric) :
  """
  Class to calculate Kendall's tau for all pairs
  """     
  
  tag = "tau"
  
  def eval(self,pred,expect) :
    pred2 = []
    expect2 = []
    l = len(pred)
    for i in range(l) :
      for j in range(i+1,l) :
        pred2.append(pred[i]-pred[j])
        expect2.append(expect[i]-expect[j])
    return _kendall(np.array(pred2),np.array(expect2),self.take)

  def init(self,pred,prederr,expect,expecterr,**kwargs) :
    # Exclude any items with zero expected differences
    self.take = []
    l = len(pred)
    n = 0
    for i in range(l) :
      for j in range(i+1,l) :
        if abs(expect[i]-expect[j]) < 0.00000001 :
          self.take.append(False)
          n = n + 1
        else :
          self.take.append(True)
    if self.verbose : 
      ntot = l*(l-1)/2.0
      print "\ntaur initialization:"
      print "Excluded pairs due to zero exp. diff. = %d of %d"%(n,ntot)
      print "Included pairs = %d"%(ntot-n)

  def lower(self,value) :
    return max(value,0.0)
  
  def name(self) :
    return "tau"
    
  def upper(self,value) :
    return min(value,1.0)

class _Qtaur(Qmetric) :
  """
  Class to calculate Kendall's tau for given pairs
  """ 
  
  tag = "taur"
     
  def eval(self,pred,expect) :
    return _kendall(pred,expect,self.take)

  def init(self,pred,prederr,expect,expecterr,**kwargs) :
    # Exclude any items with zero expected differences
    self.take = []
    l = len(pred)
    n = 0
    for e in expect :
      if abs(e) < 0.00000001 :
        self.take.append(False)
        n = n + 1
      else :
        self.take.append(True)
    if self.verbose :
      print  "\ntaur initialization:"
      print "Excluded pairs due to zero exp. diff. = %d of %d"%(n,l)
      print "Included pairs = %d"%(l-n)

  def lower(self,value) :
    return max(value,0.0)
  
  def name(self) :
    return "taur"
    
  def upper(self,value) :
    return min(value,1.0)

class _Qtaurx(Qmetric) :
  """
  Class to calculate Kendall's tau for given pairs that are significant
  """  
  
  tag = "taurx"
    
  def eval(self,pred,expect) :
    return _kendall(pred,expect,self.take)

  def init(self,pred,prederr,expect,expecterr,**kwargs) :
    # Exclude any items with non-significant differences
    self.take = []
    if "level" in kwargs :
      level = kwargs["level"]
    else :
      level = 1.64500
    l = len(pred)
    npredexcl = 0
    nexpectexcl = 0
    nexcl = 0
    if self.verbose : print "\ntaurx initialization:"
    for i in range(l) :
      takethis = True
      se = level*prederr[i]
      if abs(pred[i]) < se :
         takethis = False
         npredexcl = npredexcl + 1
         if self.verbose : print "P %d %8.3f %8.3f"%(i+1,se,pred[i])
      if len(expecterr) > 0 :
        se = level*expecterr[i]
        if abs(expect[i]) < se :
          takethis = False
          nexpectexcl = nexpectexcl + 1
          if self.verbose : print "P %d %8.3f %8.3f"%(i+1,se,expect[i])
      self.take.append(takethis)
      if not takethis : nexcl = nexcl + 1
    if self.verbose : 
      ntot = l*(l-1)/2.0
      print "Excluded pairs based on experiment = %d"%(nexpectexcl)
      print "Excluded pairs based on predictions = %d"%(npredexcl)
      print "Excluded pairs in total = %d of %d"%(nexcl,ntot)
      print "Included pairs = %d"%(ntot-nexcl)
    
  def lower(self,value) :
    return max(value,0.0)
  
  def name(self) :
    return "taurx"
    
  def upper(self,value) :
    return min(value,1.0)

class _Qtaux(Qmetric) :
  """
  Class to calculate Kendall's tau for all significant pairs
  """  
  
  tag = "taux"
     
  def eval(self,pred,expect) :
    pred2 = []
    expect2 = []
    l = len(pred)
    for i in range(l) :
      for j in range(i+1,l) :
        pred2.append(pred[i]-pred[j])
        expect2.append(expect[i]-expect[j])
    return _kendall(np.array(pred2),np.array(expect2),self.take)

  def init(self,pred,prederr,expect,expecterr,**kwargs) :
    # Exclude any items with non-significant differences
    self.take = []
    if "level" in kwargs :
      level = kwargs["level"]
    else :
      level = 1.64500
    l = len(pred)
    npredexcl = 0
    nexpectexcl = 0
    nexcl = 0
    if self.verbose : print "\ntaux initialization:"
    for i in range(l) :
      for j in range(i+1,l) :
        takethis = True
        se=level*np.sqrt(prederr[i]**2+prederr[j]**2)
        if abs(pred[i]-pred[j]) < se :
          npredexcl = npredexcl + 1
          takethis = False
          if self.verbose : print "P %d %d %8.3f %8.3f %8.3f"%(i+1,j+1,se,pred[i],pred[j])
        if len(expecterr) > 0 :
          se=level*np.sqrt(expecterr[i]**2+expecterr[j]**2)
          if abs(expect[i]-expect[j]) < se :
            nexpectexcl = nexpectexcl + 1
            takethis = False
            if self.verbose : print "E %d %d %8.3f %8.3f %8.3f"%(i+1,j+1,se,expect[i],expect[j])
        self.take.append(takethis)
        if not takethis : nexcl = nexcl + 1
    if self.verbose : 
      ntot = l*(l-1)/2.0
      print "Excluded pairs based on experiment = %d"%(nexpectexcl)
      print "Excluded pairs based on predictions = %d"%(npredexcl)
      print "Excluded pairs in total = %d of %d"%(nexcl,ntot)
      print "Included pairs = %d"%(ntot-nexcl)

  def lower(self,value) :
    return max(value,0.0)
  
  def name(self) :
    return "taux"
    
  def upper(self,value) :
    return min(value,1.0)

def _kendall(l1,l2,take=[]) :
  """
  Driver for the calculation of Kendall's tau values
  """
  if len(take) == 0 :
    return np.sum(np.multiply(np.sign(l1),np.sign(l2))) / len(l1)
  else :
    r = range(len(l1))
    ksum = 0.0
    n = 0
    for i,j in it.izip(r,r) :
      if take[i] : 
        ksum = ksum + np.sign(l1[i])*np.sign(l2[j])
        n = n + 1
    return ksum / n

