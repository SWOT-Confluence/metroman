"""
Module for filtering estimate
"""

from numpy import empty,nan,diag,var,linalg,mean,eye,sqrt
from metroman.CalcDelta import CalcDelta
from metroman.CalcADelta import CalcADelta
from metroman.CalcB import CalcB

import sys

def FilterEstimate(Estimate,C,DAll,AllObs):
    
    Deltax = CalcDelta(DAll.nR,DAll.nt,DAll.L)
    DeltaA = CalcADelta(DAll.nR,DAll.nt)

    B = CalcB(DAll.nR,DAll.nt)

    H=-Deltax

    y=(DeltaA@AllObs.hv) / DAll.dt * (B@AllObs.wv)

    
    Qchain=empty([DAll.nR*DAll.nt,C.N])
    Qchain[:]=nan
    for i in range(0,C.N):
        Qchain[:,i]=C.thetaAllQ[i,:,:].reshape(DAll.nR*DAll.nt)    

    
    P=diag(var(Qchain,1).reshape(DAll.nR*DAll.nt))

    R=AllObs.CA

    K=P@H.T @ linalg.inv((H@P@H.T+R))

    xminus=empty([DAll.nR*DAll.nt,C.N]); xminus[:]=nan
    xplus=empty([DAll.nR*DAll.nt,C.N]); xplus[:]=nan
    for i in range(0,C.N):
        xminus[:,i]=Qchain[:,i]    
        xplus[:,i]=xminus[:,i] + (K@(y-H@xminus[:,i].reshape([DAll.nR*DAll.nt,1]))).reshape(DAll.nR*DAll.nt)


    Qhatfv=mean(xplus[:,C.Nburn-1:C.N-1],1)
    #Estimate.QhatPostf=Qhatfv.reshape([D.nR,D.nt])    
    Estimate.QhatPostf=Qhatfv.reshape([DAll.nR,DAll.nt])    

    Ppostv=diag( (eye(DAll.nt*DAll.nR)-K@H)@P )
    Estimate.QhatPostfUnc=sqrt(Ppostv).reshape([DAll.nR,DAll.nt])
    return Estimate
