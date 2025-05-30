"""Program that runs MetroMan on sets of reaches and writes a NetCDF file that
contains A0, n, and Q time series.
"""

# Standard imports
import argparse
import json
import os
from pathlib import Path
import sys
import datetime

# Third party imports
from netCDF4 import Dataset
from numpy import linspace,reshape,diff,ones,array,empty,mean,zeros,putmask
import numpy as np

# Application imports
from metroman.CalcdA import CalcdA
from metroman.CalculateEstimates import CalculateEstimates
from metroman.GetCovMats import GetCovMats
from metroman.MetroManVariables import Domain,Observations,Chain,RandomSeeds,Experiment,Prior,Estimates
from metroman.MetropolisCalculations import MetropolisCalculations
from metroman.ProcessPrior import ProcessPrior
from metroman.SelObs import SelObs


def get_reachids(reachjson,index_to_run,tmp_dir,sos_bucket):
    """Extract and return a list of reach identifiers from json file.
    
    Parameters
    ----------
    reachjson : str
        Path to the file that contains the list of reaches to process
    
        
    Returns
    -------
    list
        List of reaches identifiers
    """

    #index = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
    #index=1

    if index_to_run == -235:
        index=int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
    else:
        index=index_to_run
  
    with open(reachjson) as jsonfile:
        data = json.load(jsonfile)

    if sos_bucket:
        sos_file = tmp_dir.joinpath(data[index][0]["sos"])    # just grab first in set
        if sos_bucket != 'local':
            from sos_read.sos_read import download_sos
            download_sos(sos_bucket, sos_file)
    
    return data[index]

def get_domain_obs(nr):
    """Define and return domain and observations as a tuple."""
    #this function for now is obsolete. may decide to revive it, so leaving here for now

    DAll=Domain()
    DAll.nR=nr #number of reaches
    xkm=array([25e3, 50e3, 75e3, 25e3, 50e3, 75e3])
    DAll.xkm=xkm[0:nr] #reach midpoint distance downstream [m]
    L=array([50e3, 50e3, 50e3, 50e3, 50e3, 50e3])
    DAll.L=L[0:nr]  #reach lengths, [m]
    DAll.nt=25 #number of overpasses

    AllObs=Observations(DAll)
    AllObs.sigS=1.7e-5
    AllObs.sigh=0.1
    AllObs.sigw=10

    return DAll, AllObs

def retrieve_obs(reachlist, inputdir, sosdir, Verbose):
    """ Retrieves data from SWOT and SoS files, populates observation object and
    returns : Qbar,iDelete,nDelete,SetQuality,DAll,AllObs,overlap_ts
        overlap_ts: the overlapping time indices without any bad data removed
        SetQuality : 
            0 : nominal: 3+ reaches, 4+ times (make nRmin a parameter, and link everything to ntmin) 
            1 : fewer reaches or fewer observation times
            2 : no reach with any good observation time

    -1: figure out overlapping times among reaches
     0: set up domain
     1: read observations
        1.1: check that there are enough reaches and times. if not, exit
        1.2: loop over swot input files and extract data
            1.2.1: open swot file
            1.2.2: check for nt consistency. if not exit. this should never happen
            1.2.3 read height, width and slope
            1.2.4 read Qbar from SOS
            1.2.5 read reach length and flow distance
     2: select non-fill observations

"""

    # parameters for nominal function
    ntmin=4
    ntmin_prior=4 #crashes when there are only 3 observations, weirdly... could be fixed in SelObs?
    nRmin=3
 
    # -1. figure out times by looking across all reaches and when they are measured
    #       this solution is klugey. replace once pass ids available in swot ts files

    # extract measured times
    allts=dict()
    for reach in reachlist:
        swotfile=inputdir.joinpath('swot', reach["swot"])
        swot_file_exists=os.path.exists(swotfile)
        if swot_file_exists:
            swot_dataset = Dataset(swotfile)
            nt_reach=swot_dataset.dimensions["nt"].size
            ts = list(np.round(swot_dataset["reach"]["time"][:].filled(np.nan)/3600.)) #hours
            allts[reach['reach_id']]=ts
            swot_dataset.close()
        else:
            nt_reach=0

    # determine overlapping measured times and a filter for each reach indicating which times
    #    for that reach are in the overlap
    overlap_fs=dict()
    overlap_ts = []
    for i,reach in enumerate(reachlist):
        if i==0:
            treach=[t for t in allts[reach['reach_id']]  if not np.isnan(t) ]
            overlap_ts=treach
            nt_reach=len(allts[reach['reach_id']])
            overlap_fs[reach['reach_id']]=np.full( (nt_reach,),True )
        else:
            treach=[t for t in allts[reach['reach_id']]  if not np.isnan(t) ]
            overlap_ts=list( set(overlap_ts).intersection(set(treach)))
            nt_reach=len(list(allts[reach['reach_id']]))
            overlap_fs[reach['reach_id']]=np.full( (nt_reach,),False)
            for ii,t in enumerate(list(allts[reach['reach_id']])):
                if t in overlap_ts:
                    #print('reach=',reach['reach_id'],'t=',t,'found in overlap')
                    overlap_fs[reach['reach_id']][ii]=True
        if Verbose:
            print('for reach',reach['reach_id'],'nt=',nt_reach)
            

    #loop over all reaches and fix overlap_fs if needed
    for i,reach in enumerate(reachlist):
        tarray=np.array(allts[reach['reach_id']])
        #for ii,t in enumerate(list(allts[reach['reach_id']])):
        for ii,t in enumerate(list(tarray)):
            if t in overlap_ts:
                overlap_fs[reach['reach_id']][ii]=True
            else:
                overlap_fs[reach['reach_id']][ii]=False
            # check if this is a duplicated value
            if t in list(tarray[:ii]):
                overlap_fs[reach['reach_id']][ii]=False

    overlap_ts.sort()
    nt=len(overlap_ts)

    # 0. set up domain - this could be moved to a separate function
    nr=len(reachlist)   

    DAll=Domain()
    DAll.nR=nr #number of reaches
    DAll.nt=nt

    if Verbose:
        print('Number of reaches:',nr)
        print('Total number of times (after intersecting among reach timeseries):',nt)
        print('overlapping times are: ',overlap_ts)

    # tall = [ datetime.datetime.strptime(str(t), "%Y%m%d") for t in ts ]
    epoch = datetime.datetime(2000,1,1,0,0,0)
    #tall = [ epoch + datetime.timedelta(seconds=int(t)) for t in ts ]
    tall = [ epoch + datetime.timedelta(hours=int(t)) for t in overlap_ts ]

    #print(tall)

    talli=empty(DAll.nt)
    for i in range(DAll.nt):
        dt=(tall[i]-tall[0])
        talli[i]=dt.days + dt.seconds/86400.

    if swot_file_exists:   
        AllObs=Observations(DAll)
        AllObs.sigS=1.7e-5
        AllObs.sigh=0.1
        #AllObs.sigw=10
        AllObs.sigw=20
    else:
        AllObs=0.

    # 1. read observations
    Qbar=empty(DAll.nR)
    reach_length=empty(DAll.nR)
    dist_out=empty(DAll.nR)

    if DAll.nR<nRmin:
        SetQuality=1 #initialize to deficient 
    else:
        SetQuality=0 #initialize to nominal

    # 1.1 check that there are at least some observation times
    if DAll.nt<ntmin_prior:
        if Verbose:
            print('Data issue')
            print('nR=',DAll.nR)
            print('nt=',DAll.nt)
        SetQuality=2
        iDelete=0
        nDelete=0
        return Qbar,iDelete,nDelete,SetQuality,DAll,AllObs,overlap_ts

  
    # 1.2 loop over files and extract data
    i=0
    for reach in reachlist:
        #1.2.1 open swot file
        swotfile=inputdir.joinpath('swot', reach["swot"])
        swot_file_exists=os.path.exists(swotfile)
        if swot_file_exists:
            swot_dataset = Dataset(swotfile)
            nt_reach=swot_dataset.dimensions["nt"].size
        else:
            nt_reach=0

        nt_reach_overlap=sum(overlap_fs[reach['reach_id']])

        # 1.2.2 check nt consistency
        if nt_reach_overlap != DAll.nt:
            # note this should never happen
            if Verbose:
                print('number of good observations for reach',reach['reach_id'],'does not match number of obs for set')
                print('note - this should not have happened - cause for investigation')
                print(overlap_fs[reach['reach_id']])
                print('tall=',tall)
                print('allts[reach]=',allts[reach['reach_id']])
                print('allts[reach]=',list(set(allts[reach['reach_id']])))
                print(' ... nt_reach_overlap=',nt_reach_overlap,'DAll.nt=',DAll.nt,'nt_reach=',nt_reach)
            SetQuality=2 #should test whether this is appropriate SetQuality and process prior should run
            iDelete=0
            nDelete=0
            # overlap_ts = []
            overlap_ts=list(np.delete(np.array(overlap_ts),iDelete,0))

            return Qbar,iDelete,nDelete,SetQuality,DAll,AllObs,overlap_ts 

        # 1.2.3 read height, width and slope
        h=swot_dataset["reach/wse"][0:nt_reach].filled(np.nan)
        AllObs.h[i,:]=h[overlap_fs[reach['reach_id']]]
        w=swot_dataset["reach/width"][0:nt_reach].filled(np.nan)
        AllObs.w[i,:]=w[overlap_fs[reach['reach_id']]]
        S=swot_dataset["reach/slope2"][0:nt_reach].filled(np.nan)
        AllObs.S[i,:]=S[overlap_fs[reach['reach_id']]]

        swot_dataset.close()

        # 1.2.4 read Qbar from SOS
        sosfile=sosdir.joinpath(reach["sos"])
        sos_dataset=Dataset(sosfile)
        
        sosreachids=sos_dataset["reaches/reach_id"][:]
        sosQbars=sos_dataset["model/mean_q"][:]
        k=np.argwhere(sosreachids == float(reach["reach_id"]))
        Qbar[i]=sosQbars[k].filled(np.nan) #TODO throws a warning with numpy > 1.25. fix array-scalar mapping

        if np.isnan(Qbar[i]):
             if Verbose:
                 print('Read in an invalid prior value. Stopping.')
             SetQuality=2 #cannot run process prior if Qbar is bad
        sos_dataset.close()

        #1.2.5 read reach length and flow distance
        swordfile=inputdir.joinpath('sword',reach["sword"])
        sword_dataset=Dataset(swordfile)
        swordreachids=sword_dataset["reaches/reach_id"][:]
        k=np.argwhere(swordreachids == reach["reach_id"])

        reach_lengths=sword_dataset["reaches/reach_length"][:]
        reach_length[i]=reach_lengths[k]

        dist_outs=sword_dataset["reaches/dist_out"][:]
        dist_out[i]=dist_outs[k]
        sword_dataset.close()

        i += 1

    DAll.L=reach_length
    DAll.xkm=np.max(dist_out)-dist_out + DAll.L[0]/2 #reach midpoint distance downstream [m]

    # 2. select observations that are NOT equal to the fill value

    if Verbose:
        print('before filtering bad data')
        print('nt=',DAll.nt)
        print('h=',AllObs.h)
        print('w=',AllObs.w)
        print('S=',AllObs.S)

    iDelete=np.where( np.any(np.isnan(AllObs.h),0) | np.any(np.isnan(AllObs.w),0) | np.any(np.isnan(AllObs.S),0) )

    shape_iDelete=np.shape(iDelete)
    nDelete=shape_iDelete[1]
    AllObs.h=np.delete(AllObs.h,iDelete,1)
    AllObs.w=np.delete(AllObs.w,iDelete,1)
    AllObs.S=np.delete(AllObs.S,iDelete,1)

    #overlap_ts_all=overlap_ts

    #overlap_ts=list(np.delete(np.array(overlap_ts),iDelete,0))

    DAll.nt -= nDelete
    talli=np.delete(talli,iDelete)

    if Verbose:
        print('after filtering bad data')
        print('nt=',DAll.nt)
        print('h=',AllObs.h)
        print('w=',AllObs.w)
        print('S=',AllObs.S)


    if DAll.nt<ntmin:
        
        if DAll.nt<ntmin_prior:
            SetQuality=2
        else:
            SetQuality=1

        if Verbose:
            print('After removing bad data, there are ',DAll.nt,'observations remaining. Quitting. SetQuality=',SetQuality)

        iDelete=0
        nDelete=0
        
        return Qbar,iDelete,nDelete,SetQuality,DAll,AllObs, overlap_ts	
	
    DAll.dt=empty(DAll.nt-1)
    for i in range(DAll.nt-1):
         DAll.dt[i]=(talli[i+1]-talli[i])*86400

    DAll.t=reshape(talli,[1,DAll.nt])

    # Reshape observations
    AllObs.hv=reshape(AllObs.h, (DAll.nR*DAll.nt,1))
    AllObs.Sv=reshape(AllObs.S, (DAll.nR*DAll.nt,1))
    AllObs.wv=reshape(AllObs.w, (DAll.nR*DAll.nt,1))
    return Qbar,iDelete,nDelete,SetQuality,DAll,AllObs,overlap_ts

def set_up_experiment(DAll, Qbar):
    """Define and set parameters for experiment and return a tuple of 
    Chain, Random Seed, Expriment and Prior."""

    C=Chain()
    C.N=10000
    C.Nburn=2000
    R=RandomSeeds()
    R.Seed=9

    Exp=Experiment()

    tUseMax=min(31,DAll.nt)

    #Exp.tUse=array([1,	31])
    Exp.tUse=array([1,	tUseMax-1]) #should this be 0,tUse? why is it tUseMax-1?

    Exp.nOpt=5

    P=Prior(DAll)
    P.meanQbar=mean(Qbar)
    P.covQbar=0.5
    P.eQm=0.
    P.Geomorph.Use=False
    # this is for Laterals=false
    P.AllLats.q=zeros((DAll.nR,DAll.nt))
    return C, R, Exp, P

def process(DAll, AllObs, Exp, P, R, C, Verbose,SetQuality):
    """ 
          Process observations and priors and return an estimate.

          Note - there is some sloppiness here which I think was unintended. We 
              seem to be using both "P" and "Prior" below. Pretty sure when we 
              use "Prior", we're initializing a new instance of the class, 
              rather than using the already initialized one. 

    """
    
    D,Obs,AllObs,DAll,Truth,Prior.Lats.q=SelObs(DAll,AllObs,Exp,[],Prior.AllLats)
    Prior.Lats.qv=reshape(Prior.Lats.q,(D.nR*(D.nt-1),1))
    Obs=CalcdA(D,Obs)
    AllObs=CalcdA(DAll,AllObs)
    ShowFigs=False
    DebugMode=False

    #Smin=1.7e-5
    Smin=5e-5
    Obs.S[Obs.S<Smin]=putmask(Obs.S,Obs.S<Smin,Smin) #limit slopes to a minimum value
    AllObs.S[AllObs.S<Smin]=putmask(AllObs.S,AllObs.S<Smin,Smin)

    P,jmp=ProcessPrior(P,AllObs,DAll,Obs,D,ShowFigs,Exp,R,DebugMode,Verbose)

    if SetQuality == 1:
        print('SetQuality=',SetQuality,'. returning from process function with only priors')
        Estimate,C=CalculateEstimates(None,D,Obs,P,DAll,AllObs,Exp.nOpt) 
        return Estimate,P

    Obs,P2=GetCovMats(D,Obs,Prior) #why does this use Prior instead of P? what is P2?

    C=MetropolisCalculations(P,D,Obs,jmp,C,R,DAll,AllObs,Exp.nOpt,DebugMode,Verbose)
    Estimate,C=CalculateEstimates(C,D,Obs,P,DAll,AllObs,Exp.nOpt) 
    return Estimate,P

def write_output(outputdir, reachids, Estimate, Prior,iDelete, nDelete, SetQuality,overlap_ts):
    """Write data from MetroMan run to NetCDF file in output directory. 
     
       SetQuality:  
           0 or 1: data are valid
           2: data invalid

    """

    fillvalue = -999999999999

    # 1. add back in placeholders for removed data
    iInsert=iDelete-np.arange(nDelete)
    iInsert=reshape(iInsert,[nDelete,])

    if SetQuality == 0:
        Estimate.AllQ=np.insert(Estimate.AllQ,iInsert,fillvalue,1)
        Estimate.QhatUnc_HatAllAll=np.insert(Estimate.QhatUnc_HatAllAll,iInsert,fillvalue,1)
    elif SetQuality == 1:
        Estimate.QhatAllPrior=np.insert(Estimate.QhatAllPrior,iInsert,fillvalue,1)

    # 2. write output file 
    # 2.1 file metadata
    setid = '-'.join(reachids) + "_metroman.nc"
    outfile = outputdir.joinpath(setid)

    dataset = Dataset(outfile, 'w', format="NETCDF4")
    dataset.set_id = setid    
    dataset.valid =  1 if SetQuality < 2 else 0   # TODO decide what's valid if applicable
    dataset.createDimension("nr", len(reachids))
    if SetQuality == 0:
        dataset.createDimension("nt", len(Estimate.AllQ[0]))
    elif SetQuality == 1:
        dataset.createDimension("nt", len(Estimate.QhatAllPrior[0]))
    elif SetQuality ==2:
        dataset.createDimension("nt", 0)

    # 2.2 dimensions 
    nr = dataset.createVariable("nr", "i4", ("nr",))
    nr.units = "reach"
    nr[:] = range(1, len(Estimate.A0hat) + 1)

    nt = dataset.createVariable("nt", "i4", ("nt",))
    nt.units = "time steps"
    if SetQuality == 0:
        nt[:] = range(len(Estimate.AllQ[0]))
    elif SetQuality == 1:
        nt[:] = range(len(Estimate.QhatAllPrior[0]))
    elif SetQuality == 2:
        nt[:] = range(0)

    reach_id = dataset.createVariable("reach_id", "i8", ("nr",))
    reach_id[:] = np.array(reachids, dtype=int)

    # 2.3 variables
    # 2.3.1. initialize variables
    A0 = dataset.createVariable("A0hat", "f8", ("nr",), fill_value=fillvalue)
    na = dataset.createVariable("nahat", "f8", ("nr",), fill_value=fillvalue)
    x1 = dataset.createVariable("x1hat", "f8", ("nr",), fill_value=fillvalue)
    allq = dataset.createVariable("allq", "f8", ("nr", "nt"), fill_value=fillvalue)
    qu = dataset.createVariable("q_u", "f8", ("nr", "nt"), fill_value=fillvalue)

    # 2.3.2. set variables
    if SetQuality==0: #nominal
        A0[:] = Estimate.A0hat
        na[:] = Estimate.nahat
        x1[:] = Estimate.x1hat
        allq[:] = Estimate.AllQ
        qu[:] = Estimate.QhatUnc_HatAllAll
    elif SetQuality ==1: #estimates  based on process prior only
        A0[:]=Prior.meanA0
        na[:]=Prior.meanna
        x1[:]=Prior.meanx1
        allq[:]=Estimate.QhatAllPrior
        qu[:]=fillvalue

    # 2.4 set time 
    t = dataset.createVariable("t","f8",("nt"),fill_value=fillvalue)
    t.long_name= 'swot timeseries "time" variable converted to hours and rounded to integer'
    t[:] = overlap_ts

    dataset.close()

def create_args():
    """Create and return argparsers with command line arguments."""
    
    arg_parser = argparse.ArgumentParser(description='Integrate FLPE')
    arg_parser.add_argument('-i',
                            '--index',
                            type=int,
                            default = -235,
                            help='Index to specify input data to execute on')
    arg_parser.add_argument('-r',
                            '--reachjson',
                            type=str,
                            help='Name of the reach.json',
                            default='metrosets.json')
    arg_parser.add_argument('-v',
                            '--verbose',
                            help='Indicates verbose logging',
                            action='store_true')
    arg_parser.add_argument('-s',
                            '--sosbucket',
                            type=str,
                            help='Name of the SoS bucket and key to download from, set to "local" to not download an sos',
                            default='')
    return arg_parser

def main():

    # 0 control steps
    arg_parser = create_args()
    args = arg_parser.parse_args()
    print('index: ', args.index)
    print('reach file: ', args.reachjson)
    print('verbose flag: ', args.verbose)
    print('sosbucket: ', args.sosbucket)

    # 0.1 determine the verbose flag
    if args.verbose:
        Verbose=True
    else:
        Verbose=False

    # 0.2 specify index to run. pull from command line arg or set to default = AWS
    if args.index == -235:
        index_to_run = int(os.environ.get("AWS_BATCH_JOB_ARRAY_INDEX"))
    
    else:
        index_to_run=args.index

    # 0.3 specify i/o directories
    if index_to_run == -235 or "AWS_BATCH_JOB_ID" in os.environ:
        inputdir = Path("/mnt/data/input")    
        outputdir = Path("/mnt/data/output/sets")
        if args.sosbucket != 'local':
            tmpdir = Path("/tmp")
        else:
            tmpdir = Path(os.path.join(inputdir, 'sos'))
    else:
        inputdir = Path("/home/mdurand_umass_edu/dev-confluence/mnt/input")
        outputdir = Path("/home/mdurand_umass_edu/dev-confluence/mnt/output")
        tmpdir = Path("/home/mdurand_umass_edu/dev-confluence/mnt/tmp")

    # 1 get data
    # 1.0 figure out json file. pull from command line arg or set to default
    reachjson = inputdir.joinpath(args.reachjson)

    # 1.1  get reachlist
    sos_bucket = args.sosbucket
    reachlist = get_reachids(reachjson,index_to_run,tmpdir,sos_bucket)

    if Verbose:
        print('reachlist=')
        print(reachlist)

    #1.2 retrieve obs
    if np.any(reachlist):
        if sos_bucket:
            sosdir = tmpdir
        else:
            sosdir = inputdir.joinpath("sos")
        Qbar,iDelete,nDelete,SetQuality,DAll,AllObs,overlap_ts = retrieve_obs(reachlist,inputdir,sosdir,Verbose)
    else:
        if Verbose:
            print("No reaches in list for this inversion set. ")
        print("FAIL. No good reaches. MetroMan will not run for this set. ")
        return

    # 2. run metroman
    if SetQuality == 2:
        #2.0 if inversion set is invalid, set data to fill values 
        fillvalue=-999999999999
	    #define and write fill value data
        print("FAIL. MetroMan will not run for this set. ")

        print('DAll.nt=',DAll.nt)

        DAll.nt += nDelete

        #print('iDelete=',iDelete)
        #overlap_ts=list(np.delete(np.array(overlap_ts),iDelete,0))
        
        Estimate=Estimates(DAll,DAll)
        Estimate.nahat=np.full([DAll.nR],fillvalue)
        Estimate.x1hat=np.full([DAll.nR],fillvalue)
        Estimate.A0hat=np.full([DAll.nR],fillvalue)
        Estimate.QhatUnc_HatAllAll=np.full([DAll.nR,DAll.nt],fillvalue)
        Estimate.AllQ=np.full([DAll.nR,DAll.nt],fillvalue)
        overlap_ts=np.full([DAll.nt],fillvalue)

        Estimate.QhatAllPrior=np.full([DAll.nR,DAll.nt],fillvalue)

        P=Prior(DAll)
        P.meanna=np.full([DAll.nR],fillvalue)
        P.meanx1=np.full([DAll.nR],fillvalue)
        P.meanA0=np.full([DAll.nR],fillvalue)

        DAll.nt -= nDelete
        nDelete=0
        iDelete=0
    else:
        #2.1 set up experiment
        C, R, Exp, P = set_up_experiment(DAll, Qbar)

        if Verbose:
             print('Using window',Exp.tUse)

        #2.2 process experiment
        Estimate,P = process(DAll, AllObs, Exp, P, R, C, Verbose,SetQuality)
        print("SUCCESS. MetroMan ran for this set. ")
    
    #3. write output files
    reachids = [ str(e["reach_id"]) for e in reachlist ]
    write_output(outputdir, reachids, Estimate,P,iDelete,nDelete,SetQuality,overlap_ts)

if __name__ == "__main__":
   main()    
