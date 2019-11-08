%q0Dist calculates the distribution of the likelihood ratio q0 using poison
%randomized values.

%NEvents:       Wanted number of events
%sigma_etest:   The sigma used to calculate the rate
%rate:          The rate (and 1 kg-year number of expected events)
%numRuns:       Number of q0-values to be obtained
%b:             Background

function [q0Dist,sigmaBest,N]=q0Dist(NEvents,sigma_etest,rate,numRuns,b)

q0Dist=zeros(numRuns,1);    %Preallocating
sigmaBest=zeros(numRuns,1);
N=zeros(numRuns,length(rate));

sigmaBar = NEvents/sum(rate/sigma_etest);    %Sets sigma_bar corresponding to wanted number of events

%eventsPerBinSi =@ (sigma) rate/sigma_etest.*sigma;   %Sets the number of events per energy bin for given sigma

s = rate/sigma_etest.*sigmaBar;

for j=1:numRuns
    
    N_bar_sb = poissrnd(s+b);       %Random number of events from poison dist of mean s+b.
    
    
    L =@(sigma) sum(N_bar_sb.*log(rate./sigma_etest.*sigma+b)-(rate./sigma_etest.*sigma+b));      %Log of likelihood function

    maxL=@(sigma) -L(sigma);    %Made negative for fminbnd below

    options = optimset('tolX',10^-80,'MaxFunEvals',100000000,'MaxIter',100000000);
    sigmaBest(j)=10^(fminbnd(maxL,0,1e-20,options));     %Find best sigma
    
    
    if L(0)>L(sigmaBest(j))||sigmaBest(j)<0       %In case there is no local min in [0,1] set sigmabest=0
        sigmaBest(j)=0;
    end

    q0Dist(j) = -2*(L(0)-L(sigmaBest(j)));     %Calculate likelihood ratio 
    N(j,:)=N_bar_sb;
end