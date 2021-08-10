"""Program that runs MetroMan on sets of reaches and writes a NetCDF file that
contains A0, n, and Q time series.
"""

# Standard imports
import os
from pathlib import Path
import sys

# Third party imports
from netCDF4 import Dataset
from numpy import linspace,reshape,diff,ones,array,empty,mean,zeros,putmask
import numpy as np

# Application imports
from metroman.CalcdA import CalcdA
from metroman.CalculateEstimates import CalculateEstimates
from metroman.GetCovMats import GetCovMats
from metroman.MetroManVariables import Domain,Observations,Chain,RandomSeeds,Experiment,Prior
from metroman.MetropolisCalculations import MetropolisCalculations
from metroman.ProcessPrior import ProcessPrior
from metroman.SelObs import SelObs

def get_reachids(reachtxt):
    """Extract and return a list of reach identifiers from txt file."""

    index = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
    with open(reachtxt) as reader:
        lines = reader.readlines()
    return lines[index].strip().split(',')

def get_domain_obs(nr):
    """Define and return domain and observations as a tuple."""

    DAll=Domain()
    DAll.nR=nr #number of reaches
    xkm=array([25e3, 50e3, 75e3, 25e3, 50e3, 75e3])
    DAll.xkm=xkm[0:nr] #reach midpoint distance downstream [m]
    L=array([50e3, 50e3, 50e3, 50e3, 50e3, 50e3])
    DAll.L=L[0:nr]  #reach lengths, [m]
    DAll.nt=15 #number of overpasses
    DAll.t=reshape(linspace(1,15,num=15),[1,15]) #time, [days]
    DAll.dt= reshape(diff(DAll.t).T*86400 * ones((1,DAll.nR)),(DAll.nR*(DAll.nt-1),1))

    AllObs=Observations(DAll)
    AllObs.sigS=1.7e-5
    AllObs.sigh=0.1
    AllObs.sigw=10

    return DAll, AllObs

def retrieve_obs(DAll, inputdir, reachids, AllObs):
    """ Retrieves data from SWOT and SoS files, populates observation object and
    returns qbar."""

    Qbar=empty([DAll.nR])
    i=0
    for reachid in reachids:
        swotfile=inputdir.joinpath('swot', f"{reachid}_SWOT.nc")
        sosfile=inputdir.joinpath('sos', f"{reachid}_SOS.nc")
        sos_dataset = Dataset(sosfile)
        swot_dataset = Dataset(swotfile)
        AllObs.h[i,:]=swot_dataset["reach/wse"][0:DAll.nt]
        AllObs.w[i,:]=swot_dataset["reach/width"][0:DAll.nt]
        AllObs.S[i,:]=swot_dataset["reach/slope2"][0:DAll.nt]      
        Qbar[i]=sos_dataset["reach/Qhat"][:]
        i += 1
    # Reshape observations
    AllObs.hv=reshape(AllObs.h, (DAll.nR*DAll.nt,1))
    AllObs.Sv=reshape(AllObs.S, (DAll.nR*DAll.nt,1))
    AllObs.wv=reshape(AllObs.w, (DAll.nR*DAll.nt,1))
    return Qbar

def set_up_experiment(DAll, Qbar):
    """Define and set parameters for experiment and return a tuple of 
    Chain, Random Seed, Expriment and Prior."""

    C=Chain()
    C.N=10000
    C.Nburn=2000
    R=RandomSeeds()
    R.Seed=9

    Exp=Experiment()
    Exp.tUse=array([1,	31])
    Exp.nOpt=5

    P=Prior(DAll)
    P.meanQbar=mean(Qbar)
    P.covQbar=0.5
    P.eQm=0.
    P.Geomorph.Use=False
    # this is for Laterals=false
    P.AllLats.q=zeros((DAll.nR,DAll.nt))
    return C, R, Exp, P

def process(DAll, AllObs, Exp, P, R, C):
    """Process observations and priors and return an estimate."""
    
    D,Obs,AllObs,DAll,Truth,Prior.Lats.q=SelObs(DAll,AllObs,Exp,[],Prior.AllLats)
    Prior.Lats.qv=reshape(Prior.Lats.q,(D.nR*(D.nt-1),1))
    Obs=CalcdA(D,Obs)
    AllObs=CalcdA(DAll,AllObs)
    ShowFigs=False
    DebugMode=False
    P,jmp=ProcessPrior(P,AllObs,DAll,Obs,D,ShowFigs,Exp,R,DebugMode)
    Obs,P2=GetCovMats(D,Obs,Prior)
    Obs.S[Obs.S<0]=putmask(Obs.S,Obs.S<0,0) #limit slopes to zero
    AllObs.S[AllObs.S<0]=putmask(AllObs.S,AllObs.S<0,0)
    C=MetropolisCalculations(P,D,Obs,jmp,C,R,DAll,AllObs,Exp.nOpt,DebugMode)
    Estimate,C=CalculateEstimates(C,D,Obs,P,DAll,AllObs,Exp.nOpt) 
    return Estimate

def write_output(outputdir, reachids, Estimate):
    """Write data from MetroMan run to NetCDF file in output directory."""

    setid = '-'.join(reachids) + "_metroman.nc"
    outfile = outputdir.joinpath(setid)
    dataset = Dataset(outfile, 'w', format="NETCDF4")
    dataset.set_id = setid    # TODO decide on how to identify sets
    dataset.valid = 'X'    # TODO decide what's valid if applicable
    dataset.createDimension("nr", len(reachids))
    dataset.createDimension("nt", len(Estimate.AllQ[0]))
    dataset.createDimension("u", 2)    # TODO what do the two uncertainty values represent?

    fillvalue = -9999
    nr = dataset.createVariable("nr", "i4", ("nr",))
    nr.units = "reach"
    nr[:] = range(1, len(Estimate.A0hat) + 1)

    nt = dataset.createVariable("nt", "i4", ("nt",))
    nt.units = "time steps"
    nt[:] = range(len(Estimate.AllQ[0]))

    reach_id = dataset.createVariable("reach_id", "i8", ("nr",))
    reach_id[:] = np.array(reachids, dtype=int)

    A0 = dataset.createVariable("A0hat", "f8", ("nr",), fill_value=fillvalue)
    A0[:] = Estimate.A0hat

    n = dataset.createVariable("nahat", "f8", ("nr",), fill_value=fillvalue)
    n[:] = Estimate.nahat
    
    q = dataset.createVariable("x1hat", "f8", ("nr",), fill_value=fillvalue)
    q[:] = Estimate.x1hat

    allq = dataset.createVariable("allq", "f8", ("nr", "nt"), fill_value=fillvalue)
    allq[:] = Estimate.AllQ

    qu = dataset.createVariable("q_u", "f8", ("nr", "u"), fill_value=fillvalue)
    qu[:] = Estimate.QhatUnc_HatAll

    dataset.close()

def main():
    inputdir = Path("/mnt/data/input")
    outputdir = Path("/mnt/data/output")
    try:
        reachtxt = inputdir.joinpath(sys.argv[1])
    except IndexError:
        reachtxt = inputdir.joinpath("reach_sets.txt") 

    reachids = get_reachids(reachtxt)
    DAll, AllObs = get_domain_obs(len(reachids))
    Qbar = retrieve_obs(DAll, inputdir, reachids, AllObs)
    C, R, Exp, P = set_up_experiment(DAll, Qbar)
    Estimate = process(DAll, AllObs, Exp, P, R, C)
    
    write_output(outputdir, reachids, Estimate)

if __name__ == "__main__":
   main()    