%% This routine uses the crystal form factor column vectors of DM-FF.mat to 
%calculate the event rate, sigma detection range and required number of
%events for detection. This process is also iterated over a range of
%masses. Also small changes to the code allows it to be used to test the
%dependency of other parameters. 

load('DM-FF.mat')       %Load crystal form factor column vectors
load('Background.mat')  %Load background values
bgdc(1,2) = 1.6*10^7;

F_si=zeros(900,500);        %Preallocating form factor matrices
F_ge=zeros(900,500);

%Rearranging the column vectors into matrices with E_e on x axis and q on y
%axis. deltaE_e=0.1 eV and deltaq=0.02*alpha*m_e 
for i=1:500
    for j=1:900
        F_si(j,i)=C2((i-1)*900+j);
        F_ge(j,i)=C3((i-1)*900+j);
    end
end

format long
q=linspace(1,900,900);                              %unitless row vector
E_e=linspace(1,500,500);                            %unitless row vector

%declare constants
sec2year = 60*60*24*365.25;                         %seconds to years
amu2kg = 1.660538782*10^(-27);                      %conversion atomic units to kg
alpha = 1/137;                                      %fine structure const
m_e = 0.511*1e6;                                    %electron mass
rho_chi = 0.4*1e9;                                  %local dm density [eV/cm^3]
c = 2.99792458*1e10;                                %speed of light [cm/s]
v0 = 230*1e5;                                       %typical velocity [cm/s]
vE = 240*1e5;                                       %velocity earth [cm/s]
vesc = 600*1e5;                                     %escape velocity earth


F_DM =@(q,n) (alpha*m_e/q)^n;                       %F_DM function


sigma_etest = 1e-37;                                %testparam for crossection

MTarget = 1;                                        %Mass of target detector, kg

qunit = 0.02*alpha*m_e;                             %unit for q vector
Eunit = 0.1;                                        %unit for E_e vector
Emin_Si = 1.1;                                      %bandgap energy for Si
Emin_Ge = 0.67;
iterations = 50;
temp1 = ones(1,10)*1500;
temp2 = ones(1,10)*500;
temp3 = ones(1,60)*350;
temp4 = ones(1,20)*350;
%temp5 = ones(1,10)*100;
expectedMaxNReq = [5e5 3e5 8e4 2e4 8e3 4e3 2e3 1.7e3 1.3e3];
expectedMaxNReq = [expectedMaxNReq  temp1 temp2 temp3 temp4];
k = expectedMaxNReq./iterations;

for n=2:2:2
for exposureIteration = 1:1
%Choose F_DM n=(0 -> F_DM=1, 1 -> F_DM=1/q, 2 -> F_DM=1/q^2)
%n=2;
l=[linspace(1,9,9),linspace(10,1000,100)];  %masses to be iterated over (in MeV)
%l=2;
% Choose High Voltage or CCD
HV=false;

DCHigh=true;        %Choose DC value 5e-3(true) or 1e-7(false)

%Set exposure time

if HV
    exposureSi=10;      %in kg*years
    exposureGe=50;     %in kg*years
else
    exposureSi=[0.1 1];     %in kg*years
    exposureGe=0.1;       %in kg*years
end

dRdESi = zeros(1,length(E_e));                     %preallocate rate array
dRdEGe = zeros(1,length(E_e));

%Semiconductur properties
epsSi=3.6;      %Linearization constant for secondary processes
epsGe=2.9;
E_gapSi=1.11;   %Bandgap
E_gapGe=0.67;

%Calc. signal Q eq. 3.14, i.e. number of excited electrons per event of
%energy E_e
QSi=1+floor((E_e*Eunit-E_gapSi)/epsSi);
QGe=1+floor((E_e*Eunit-E_gapGe)/epsGe);

Background_HV = dat(:,1);

bSi = 0.01*ones(1,max(QSi));
bGe = 0.01*ones(1,max(QGe));

if HV
    bSi(1:3)=exposureSi(exposureIteration)*Background_HV*365.24;
    bGe(1:3)=exposureGe*Background_HV*365.24;
else if DCHigh
        bSi=0.01*ones(1,max(QSi));
        bGe=0.01*ones(1,max(QGe));
        bSi(1:4)=exposureSi(exposureIteration)*bgdc(:,2)*365.24;      %Si background
        bGe(1:4)=exposureGe*bgdc(:,2)*365.24;      %Ge background
    else
        bSi=0.01*ones(1,max(QSi));
        bGe=0.01*ones(1,max(QGe));
        bSi(1:4)=exposureSi*bgdc(:,3)*365.24;      %Si background
        bGe(1:4)=exposureGe*bgdc(:,3)*365.24;      %Ge background
    end
end


count=0;    %basic counter to avoid lots of zeros in matrices due to masses not being linear.

%Preallocating resluts vectors
N5reqSi=zeros(1,length(l));
N3reqSi=zeros(1,length(l));
NreqGe=zeros(1,length(l));
NreqComb=zeros(1,length(l));
sig5Si=zeros(1,length(l));
sig3Si=zeros(1,length(l));
sigGe=zeros(1,length(l));
sigComb=zeros(1,length(l));

%Start of iterations
for iter=100:length(l)
    j = l(iter);
    startClock=clock;   %Clock for estimating time left
    count=count+1;      %Count upwards, start on 1
    
    m_chi=(1*j)*10^6;   %Set mass to MeV
    mu = (m_e*m_chi)/(m_e+m_chi);                       %reduced mass
    v_min =@(q,E_e,m_chi) c.*(q./(2.*m_chi)+E_e./q);    %V_min function
    
    %Computes the diff. event rate dR/dE by summing over all q-values in each
    %energy bin
    for i=round(Emin_Si/Eunit):length(E_e)
        
        dRdESi(i) = exposureSi(exposureIteration)*dRdE_Si(m_chi,i,n,rho_chi,MTarget,sigma_etest,qunit,Eunit,q,F_si);
    end

    Events_Si = round(sum(dRdESi));     %Just a check of total number of 


    RSi=zeros(1,500);       %Preallocating rates for each energy bin


    %Sum interval and set all points in interval to value of sum
    %(for plotting the rates in different intervals)

    if HV
        for i=1:max(QGe)
            RSi=RSi+0.9*sum(dRdESi.*(QSi==i)).*(QSi==i);
            RGe=RGe+0.9*sum(dRdEGe.*(QGe==i)).*(QGe==i);
        end
    else
        eff15=[1 0.62 0.48 0.41 0.36];
        for i=1:5
            RSi=RSi+eff15(i)*sum(dRdESi.*(QSi==i)).*(QSi==i);
        end
        
        for i=6:max(QGe)
            RSi=RSi+0.36*sum(dRdESi.*(QSi==i)).*(QSi==i);
        end
    end


    
    %Create vectors with only one value for each energy interval
    for i=1:max(QSi)
        tmp=find(QSi==i);
        TOTrE_Si(i)=RSi(tmp(1));
    end
    
    
    %Start values to remove zeros in beginning of plot
    %(the log will take care of it anyway but...)
    startSi=length(QSi)-length(nonzeros(QSi))+1;
    RSi=RSi(startSi:length(RSi));
    
   % Start of the statistical calculations
    numberOfSimulations=1000;       %Set number of q0 simulations

    rateSi=TOTrE_Si;        %Si rate input
 
    sigmaMax=10^-35;

    
    wantedNbrOfEvents=zeros(1,iterations);              %Preallocating
        %Step size
    %Calculating significance for different numbers of events.
    for i=1:iterations

        wantedNbrOfEvents(i) =k(iter)*i;         %Number of measured events

        [q0Si,sigmaBestSi,NSi]=q0Dist(wantedNbrOfEvents(i),sigma_etest,rateSi,numberOfSimulations,bSi);
        %[q0Ge,sigmaBestGe,NGe]=q0Dist(wantedNbrOfEvents(i),sigma_etest,rateGe,numberOfSimulations,bGe);
        %q0Comb=q0DistComb(wantedNbrOfEvents(i),sigma_etest,rateSi,rateGe,numberOfSimulations,bSi,bGe);

        ZSi(count,i)=sqrt(median(q0Si));    %Calculate Z-value
        %ZGe(count,i)=sqrt(median(q0Ge));
        %ZComb(count,i)=sqrt(median(q0Comb));
    end
    
    
    %ZrelComb=ZComb-5;
    %NreqComb=NreqComb+ k*(find(ZrelComb>0,1)-1);
        
    %ZrelSi=ZSi-5;
    %NreqSi=NreqSi+ k*(find(ZrelSi>0,1)-1);
        
    %ZrelGe=ZGe-5;
    %NreqGe=NreqGe+ k*(find(ZrelGe>0,1)-1);
        
        

%     aSi=polyfit(log(wantedNbrOfEvents(2:length(wantedNbrOfEvents))),log(ZSi(2:length(wantedNbrOfEvents))),1);
%     
%     NreqSi(count)=(5/exp(aSi(2)))^(1/aSi(1));
    
    Z5relSi=ZSi(count,:)-5;
    Z3relSi=ZSi(count,:)-3;
    N5reqSi(count) = k(iter)*find(Z5relSi>0,1);
    N3reqSi(count) = k(iter)*find(Z3relSi>0,1);
%     ktmpSi(count)=(ZrelSi(NtmpSi(count))-ZrelSi(NtmpSi(count)-1))/k;
%     NreqSi(count)=k*NtmpSi(count)-ZrelSi(NtmpSi(count))/ktmpSi(count);
    sig5Si(count)=N5reqSi(count)/sum(rateSi/sigma_etest);
    sig3Si(count)=N3reqSi(count)/sum(rateSi/sigma_etest);
    
    %Pretty stuff for knowing how much time is left
    disp([num2str(count/length(l)*100),'%']);

    estTimeLeft=(length(l)-count)*etime(clock,startClock);
    
    tmp=mod(estTimeLeft,3600);
    hours=round((estTimeLeft-tmp)/3600);
    seconds=mod(estTimeLeft-3600*hours,60);
    minutes=round((tmp-seconds)/60);
   
    
    disp(['Estimated time left: ',num2str(hours),'h ',num2str(minutes),'m ',num2str(seconds),'s']);
    
 end
N5reqSiSave(exposureIteration+n,:) = N5reqSi;
N3reqSiSave(exposureIteration+n,:) = N3reqSi;
sig5SiSave(exposureIteration+n,:) = sig5Si;
sig3SiSave(exposureIteration+n,:) = sig3Si;
end 
end
%% Plotting
figure(1)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

plot(l,sig5Si,'blue','LineWidth',1);
hold on;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('$m_{\chi}$ (MeV)');
ylabel('$\overline{\sigma}_e$ (cm$^2$)');
set(gca,'fontsize', 16);

% plot(l,sigGe,':','Color','red','LineWidth',1.7);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');

legend('Si', 'Ge');

grid on
hold off;
%% Plotting real
figure(2)

subplot(1,2,2)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
plot(l,N5reqSi,'Color',[0,0,1],'LineWidth',1);
hold on
set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
set(gca,'fontsize', 16);

%axis([1 1000 1e-42 1e-38]);

plot(l,NreqGe,':','Color',[1,0,0],'LineWidth',1.5);
plot(l,NreqComb,'--','Color',[0,0.6,0],'LineWidth',1.5);
legend('Si','Ge','Combined');
grid on
xlabel('$m_{\chi}$ (MeV)');
ylabel('$\overline{\sigma}_e$ (cm$^2$)');

t=text(2,2e-39,'F$_\mathrm{DM}=1$');
%t=text(1.5,2e-39,'F$_\mathrm{DM}=\left(\frac{\alpha m_e}{q}\right)^2$');
t.FontSize=16;



%% Plot Nreq from sigreq

figure(2)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
plot(l,sig5Si*sum(rateSi/sigma_etest),'Color',[0,0,1],'LineWidth',1);
hold on
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
set(gca,'fontsize', 16);

plot(l,sigGe*sum(rateGe/sigma_etest),':','Color',[1,0,0],'LineWidth',1.5);
plot(l,sigComb*sum(rateGe/sigma_etest),'--','Color',[0,0.6,0],'LineWidth',1.5);
legend('Si','Ge','Combined');
grid on
xlabel('$m_{\chi}$ (MeV)');
ylabel('$\overline{\sigma}_e$ (cm$^2$)');

t=text(2,2e-39,'F$_\mathrm{DM}=\left(\frac{\alpha m_e}{q}\right)^2$');
t.FontSize=16;
hold off


%% Plot Z vs N
figure(3)

plot((0:19),ZSi,'o','color','b');
hold on

plot(linspace(0,20),exp(aSi(2)).*linspace(0,20).^(aSi(1)),'b');

%plot((0:19),ZComb,'x','color','b');

plot(linspace(0,20),5*ones(100,1),'--','color','r');
grid on
plot(N5reqSi*ones(100,1),linspace(0,8),':','color','r');

text(14,4.5,[num2str(N5reqSi),' events']);
legend('Simulated data','Fitted approximation')
hold off;

%% Comparing diff between Si and Ge with improvement of combining

figure(2)

subplot(1,2,1)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

for o=1:109
    if sig5Si(o)>sigGe(o)
        diffvec(o)=sigGe(o)-sigComb(o);
    else
        diffvec(o)=sig5Si(o)-sigComb(o);
    end
    
end



plot(abs(sig5Si-sigGe),sig5Si-sigComb,'Color',[0,0,1],'LineWidth',1);
hold on
%set(gca, 'XScale', 'log');
%set(gca, 'YScale', 'log');
set(gca,'fontsize', 16);

%plot(l,sigGe,':','Color',[1,0,0],'LineWidth',1.5);
%plot(l,sigComb,'--','Color',[0,0.6,0],'LineWidth',1.5);
legend('Si','Ge','Combined');
grid on
xlabel('$m_{\chi}$ (MeV)');
ylabel('$\overline{\sigma}_e$ (cm$^2$)');

t=text(2,2e-39,'F$_\mathrm{DM}=1$');
%t=text(1.5,2e-39,'F$_\mathrm{DM}=\left(\frac{\alpha m_e}{q}\right)^2$');
t.FontSize=16;