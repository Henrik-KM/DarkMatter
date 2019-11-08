%q0DistComb calculates the distribution of the combined likelihood ratio q0 for two experiments using poison
%randomized values.

%NEvents:       Wanted number of events
%sigma_etest:   The sigma used to calculate the rate
%rateSi:        The Si rate (and 1 kg-year number of expected events)
%rateGe:        The Ge rate (and 1 kg-year number of expected events)
%numRuns:       Number of q0-values to be obtained
%bSi:           Si Background
%bGe:           Ge Background

function [q0DistComb,sigmaBest]=q0DistComb(NEvents,sigma_etest,rateSi,rateGe,numRuns,bSi,bGe)

q0DistComb=zeros(numRuns,1);    %Preallocating

for j=1:numRuns
    
    sigmaBar = NEvents/sum(rateSi/sigma_etest);  %Sets sigma_bar corresponding to wanted number of Si events

    %Set number of events per energy bin for given sigma for Si and Ge resp.
    eventsPerBinSi =@ (sigma) rateSi/sigma_etest*sigma; 
    eventsPerBinGe =@ (sigma) rateGe/sigma_etest*sigma;

    sSi = eventsPerBinSi(sigmaBar);
    N_bar_sb_Si = poissrnd(sSi+bSi);    %Random number of Si events from poison dist of mean s+b.

    sGe = eventsPerBinGe(sigmaBar);
    N_bar_sb_Ge = poissrnd(sGe+bGe);    %Random number of Ge events from poison dist of mean s+b.
    
    
    %Log of likelihood functions for Si and Ge separately
    LSi =@(sigma) sum(N_bar_sb_Si.*log(eventsPerBinSi(sigma)+bSi)-(eventsPerBinSi(sigma)+bSi));
    LGe =@(sigma) sum(N_bar_sb_Ge.*log(eventsPerBinGe(sigma)+bGe)-(eventsPerBinGe(sigma)+bGe));
    
    %Log of product of likelihood functions
    L=@(sigma) LSi(sigma)+LGe(sigma);

    maxL=@(sigma) -L(sigma);    %Made negative for fminbnd below

    options = optimset('tolX',10^-100,'MaxFunEvals',100000000,'MaxIter',100000000);
    sigmaBest(j)=fminbnd(maxL,0,1,options);     %Find best sigma

    if L(0)>L(sigmaBest(j))||sigmaBest(j)<0    %In case there is no local min in [0,1] set sigmabest=0
        sigmaBest(j)=0;
    end
   
    q0DistComb(j) = -2*(L(0)-L(sigmaBest(j)));     %Calculate likelihood ratio
end
