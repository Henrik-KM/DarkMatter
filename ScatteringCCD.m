%% This routine uses the crystal form factor column vectors of DM-FF.mat to
%calculate the event rate, sigma detection range and required number of
%events for detection. This process is also iterated over a range of
%masses. Also small changes to the code allows it to be used to test the
%dependency of other parameters.

clear
clc

load('DM-FF.mat')       %Load crystal form factor column vectors
load('Background.mat')  %Load background values

bgdc(1,2)=1.6e7;


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

v0init =  230;
vEinit = 240;
vescinit = 600;

v0Array = linspace(200,260,50)*1e5;                                       %typical velocity [cm/s]
vEArray = linspace(200,280,50)*1e5;                      %velocity earth [cm/s]
vescArray = linspace(500,700,50)*1e5;                                     %escape velocity earth


F_DM =@(q,n) (alpha*m_e/q)^n;                       %F_DM function


sigma_etest = 1e-37;                                %testparam for crossection

MTarget = 1;                                        %Mass of target detector, kg

qunit = 0.02*alpha*m_e;                             %unit for q vector
Eunit = 0.1;                                        %unit for E_e vector
Emin_Si = 1.1;                                      %bandgap energy for Si
Emin_Ge = 0.67;

for n=0:2:2
for relevantv = 1:3
for iterVar = 1:size(v0,2)
    if relevantv == 1
        v0 = v0Array(iterVar);
        vesc = vescinit;
        vE = vEinit;
    end
        if relevantv == 2
        v0 = v0init;
        vesc = vescArray(iterVar);
        vE = vEinit;
        end
            if relevantv == 3
        v0 = v0init;
        vesc = vescinit;
        vE = vEArray(iterVar);
    end

for CCD_low = 1
    %Choose F_DM n=(0 -> F_DM=1, 1 -> F_DM=1/q, 2 -> F_DM=1/q^2)   
    % Choose High Voltage or CCD
    HV=false;
    
    DCHigh=false;        %Choose DC value 5e-3(true) or 1e-7(false)
    
    %Set exposure time
    
    if CCD_low == 0
        exposureSi=0.1;     %in kg*years
    end
    if CCD_low == 1
        exposureSi = 1;
    end
    
    
    dRdESi = zeros(1,length(E_e));                     %preallocate rate array
    
    %Semiconductur properties
    epsSi=3.6;      %Linearization constant for secondary processes
    epsGe=2.9;
    E_gapSi=1.11;   %Bandgap
    E_gapGe=0.67;
    
    %Calc. signal Q eq. 3.14, i.e. number of excited electrons per event of
    %energy E_e
    QSi=1+floor((E_e*Eunit-E_gapSi)/epsSi);
    
    bSi=exposureSi*background(E_e,Eunit,QSi,true);

    if DCHigh
        bSi(1:4)=exposureSi*bgdc(:,2)*365.24;      %Si background
        %bGe(1:4)=exposureGe*bgdc(:,2)*365.24;      %Ge background
    end
    if DCHigh == false
        bSi(1:4)=bSi(1:4)+exposureSi*bgdc(:,3)'*365.24;     %Si background
        %bGe(1:4)=exposureGe*bgdc(:,3)*365.24;      %Ge background
        
    end
    
    l=[linspace(1,9,9),linspace(10,1000,100)];  %masses to be iterated over (in MeV)
    %l=1;
    
    Z = 5;
  
    count=0;    %basic counter to avoid lots of zeros in matrices due to masses not being linear.
    
    
           %Step size
        if CCD_low==1 
        k_values = [10,10,10,10,10,10,10,0*linspace(1,1,102)+1];
        startbiases = [900,0*linspace(1,1,108)];
        iter = [100,100,100,100,60,60,60*linspace(1,1,103)];
        end
        if CCD_low==0 
        k_values = [10,linspace(1,1,108)];
        startbiases = [200,0*linspace(1,1,108)];
        iter = [100,300,100,100,80,60,50*linspace(1,1,103)];
        end
        
    %Preallocating resluts vectors
    NreqSi_5=zeros(1,length(l));
    sigSi_5=zeros(1,length(l));
    
    %Start of iterations
    for r=19
        j = l(r);
        startClock=clock;   %Clock for estimating time left
        count=count+1;      %Count upwards, start on 1
        
        m_chi=(1*j)*10^6;   %Set mass to MeV
        mu = (m_e*m_chi)/(m_e+m_chi);                       %reduced mass
        v_min =@(q,E_e,m_chi) c.*(q./(2.*m_chi)+E_e./q);    %V_min function
        
        %Computes the diff. event rate dR/dE by summing over all q-values in each
        %energy bin
        for i=round(Emin_Si/Eunit):length(E_e)
            
            dRdESi(i) = exposureSi*dRdE_Si(m_chi,i,n,rho_chi,MTarget,sigma_etest,qunit,Eunit,q,F_si,v0,vE,vesc);
        end
        
        %for i =round(Emin_Ge/Eunit):length(E_e)
        
        %    dRdEGe(i)= exposureGe*dRdE_Ge(m_chi,i,n,rho_chi,MTarget,sigma_etest,qunit,Eunit,q,F_ge);
        %end
        
        Events_Si = round(sum(dRdESi));     %Just a check of total number of
        %Events_Ge = round(sum(dRdEGe));     % events for 1 kg year exposure
        
        
        RSi=zeros(1,500);       %Preallocating rates for each energy bin
        %RGe=zeros(1,500);
        
        
        %Sum interval and set all points in interval to value of sum
        %(for plotting the rates in different intervals)
        

        if HV == false
            eff15=[1 0.62 0.48 0.41 0.36];
            for i=1:5
                RSi=RSi+eff15(i)*sum(dRdESi.*(QSi==i)).*(QSi==i);
                %RGe=RGe+eff15(i)*sum(dRdEGe.*(QGe==i)).*(QGe==i);
            end
            
            for i=6:max(QSi)
                RSi=RSi+0.36*sum(dRdESi.*(QSi==i)).*(QSi==i);
                %RGe=RGe+0.36*sum(dRdEGe.*(QGe==i)).*(QGe==i);
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
        %rateGe=TOTrE_Ge;        %Ge rate input
        
        sigmaMax=10^-35;
        
        
        wantedNbrOfEvents=zeros(1,15);              %Preallocating
        
        
       
        
        k = k_values(r);
        startbias=startbiases(r);
        
        %Calculating significance for different numbers of events.
        for i=1:iter(r)
            
            wantedNbrOfEvents(i) =k*i+startbias;         %Number of measured events
            
            [q0Si,sigmaBestSi,NSi]=q0Dist(wantedNbrOfEvents(i),sigma_etest,rateSi,numberOfSimulations,bSi);
            %[q0Ge,sigmaBestGe,NGe]=q0Dist(wantedNbrOfEvents(i),sigma_etest,rateGe,numberOfSimulations,bGe);
            %q0Comb=q0DistComb(wantedNbrOfEvents(i),sigma_etest,rateSi,rateGe,numberOfSimulations,bSi,bGe);
            
            ZSi(count,i)=sqrt(median(q0Si));    %Calculate Z-value
            %ZGe(count,i)=sqrt(median(q0Ge));
            %ZComb(count,i)=sqrt(median(q0Comb));
        end
        
        
        
        ZrelSi_5=ZSi(iterVar,count,:)-5;

        NreqSi_5(iterVar,count) = k*find(ZrelSi_5>0,1)+startbias;
 
        sigSi_5(iterVar,count)=NreqSi_5(iterVar,count)/sum(rateSi/sigma_etest);
       
        
        
        %Pretty stuff for knowing how much time is left
        disp([num2str((count+length(l)*CCD_low)/(2*length(l))*100),'%']);
        
        estTimeLeft=(2*length(l)-count-CCD_low*length(l))*etime(clock,startClock);
        
        tmp=mod(estTimeLeft,3600);
        hours=round((estTimeLeft-tmp)/3600);
        seconds=mod(estTimeLeft-3600*hours,60);
        minutes=round((tmp-seconds)/60);
        
        disp(['Estimated time left: ',num2str(hours),'h ',num2str(minutes),'m ',num2str(seconds),'s']);
    end
    
%     if CCD_low == 0
%         save('FDM_q2_wanted_CCD_01_Z5xxx','NreqSi_5')
%         save('FDM_q2_sig_CCD_01_Z5xxx','sigSi_5')
%         save('FDM_q2_wanted_CCD_01_Z3xxx','NreqSi_3')
%         save('FDM_q2_sig_CCD_01_Z3xxx','sigSi_3')
%     end
%     if CCD_low == 1
%         save('FDM_q2_wanted_CCD_1_Z5xxx','NreqSi_5')
%         save('FDM_q2_sig_CCD_1_Z5xxx','sigSi_5')
%         save('FDM_q2_wanted_CCD_1_Z3xxx','NreqSi_3')
%         save('FDM_q2_sig_CCD_1_Z3xxx','sigSi_3')
%     end
%     save('FDM_q2_ZSi_runxxx','ZSi')
end
end
if n==0
if relevantv == 1
    v0sigF1 = sigSi_5;
    v0NreqF1 = NreqSi_5;
end
if relevantv == 2
    vescsigF1 = sigSi_5;
    vescNreqF1 = NreqSi_5;
end
if relevantv == 3
    vEsigF1 = sigSi_5;
    vENreqF1 = NreqSi_5;
end
else
    if relevantv == 1
    v0sigF1q2 = sigSi_5;
    v0NreqF1q2 = NreqSi_5;
end
if relevantv == 2
    vescsigF1q2 = sigSi_5;
    vescNreqF1q2 = NreqSi_5;
end
if relevantv == 3
    vEsigF1q2 = sigSi_5;
    vENreqF1q2 = NreqSi_5;
end
end
end

end

%% Plotting

import = load('sig_CCD_01_Z5.mat','sigSi_5');
sigSi = import.sigSi_5;

figure(1)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

plot(l,sigSi,'blue','LineWidth',1);
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
plot(l,NreqSi,'Color',[0,0,1],'LineWidth',1);
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
plot(l,sigSi*sum(rateSi/sigma_etest),'Color',[0,0,1],'LineWidth',1);
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
plot(NreqSi*ones(100,1),linspace(0,8),':','color','r');

text(14,4.5,[num2str(NreqSi),' events']);
legend('Simulated data','Fitted approximation')
hold off;

%% Comparing diff between Si and Ge with improvement of combining

figure(2)

subplot(1,2,1)
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

for o=1:109
    if sigSi(o)>sigGe(o)
        diffvec(o)=sigGe(o)-sigComb(o);
    else
        diffvec(o)=sigSi(o)-sigComb(o);
    end
    
end



plot(abs(sigSi-sigGe),sigSi-sigComb,'Color',[0,0,1],'LineWidth',1);
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

%%

mass_test = [linspace(1,9,9),linspace(10,1000,100)];

figure(1)
clf

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

sig_test_Z5_01 = load('sig_CCD_01_Z5.mat','sigSi_5');
sig_FDM_1_exp_01 = sig_test_Z5_01.sigSi_5;
sig_test_Z5_01 = load('sig_CCD_1_Z5.mat','sigSi_5');
sig_FDM_1_exp_1 = sig_test_Z5_01.sigSi_5;
sig_test_Z5_01 = load('FDM_q2_sig_CCD_01_Z5.mat','sigSi_5');
sig_FDM_q2_exp_01 = sig_test_Z5_01.sigSi_5;
sig_test_Z5_01 = load('FDM_q2_sig_CCD_1_Z5.mat','sigSi_5');
sig_FDM_q2_exp_1 = sig_test_Z5_01.sigSi_5;

plot((mass_test),(sig_FDM_1_exp_1),'b')
hold on
plot((mass_test),(sig_FDM_1_exp_01),'r')
hold on
plot((mass_test),(sig_FDM_q2_exp_1),'b--')
hold on
plot((mass_test),(sig_FDM_q2_exp_01),'r--')
hold on


set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xlabel('$m_{\chi}$ (MeV)');
ylabel('$\overline{\sigma}_e$ (cm$^2$)');
legend('$F_{DM}=1$, exposure = 1', '$F_{DM}=1$, exposure = 0.1', '$F_{DM}\propto1/q^2$, exposure = 1', '$F_{DM}\propto1/q^2$, exposure = 0.1');

grid on
hold off;