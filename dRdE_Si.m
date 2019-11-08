function dRdE_Si = dRdE_Si(m_chi,E_e,n,rho_c,MTarget,sigmatest,qunit,Eunit,qarray,ffarray,v0,vE,vesc)

% dRdE_Si Computes the number of scattering events in an Si detector.
%   m_chi:      Mass of dark matter particle.
%   E_e:        Electron energy of the interaction.
%   n:          Determines dark matter form factor. F_DM = (const/q)^n.
%   rho_c:      Dark matter density.
%   MTarget:    Weight of Si target detector (kg).
%   sigmatest:  Sigma_e in equation 3.13. Usually 1e-37.
%   qunit:      Unit of which the energies q are binned.
%   Eunit:      Unit of which the electron energies are binned.
%   qarray:     Array of q-values in interest.
%   ffarray:    Array of the crystal form factor for Si. 500 E-bins across,
%   900 q values in each.
%   When in doubt, use the standard ffarray in this package, SiArray.dat.
%   This requires qarray = linspace(1,900,900).
%   Standard qunit is 0.02/137*0.511*10^6.
%   Standard Eunit is 0.1. 

m_e = 0.511*1e6;                    %mass of electron [eV]
amu2kg = 1.660538782*10^(-27);      %conversion amu to kg
c = 2.99792458*1e10;                %speed of light
alpha = 1/137;                      %fine structure constant
wk = 2/137;                         %TBA
sec2year = 60*60*24*365.25;         %conversion sec to year
MCell = 2*28.0855*amu2kg;           %Mass of unit cell in Si
my = (m_e*m_chi)/(m_e+m_chi);       %Reduced mass
Eprefactor_Si = 2;                     %Prefactor for calculation

F_DM =@(q,n) (alpha*m_e/q)^n;       %Dark matter form factor definition
FFarraySi = ffarray;                %Just assigns the formfactor array (Si)
q=qarray;                           %Just assigns the q values, unitless.

prefac = c^2*sec2year*(rho_c/m_chi)*(MTarget/MCell)*sigmatest*(alpha*m_e^2)/(my)^2;     %Constant prefactor 

%Calculate all q-values in the energy bin.
for i=1:length(q)
    
    drde(i)=(Eunit*qunit*Eta(i*qunit,E_e*Eunit,m_chi,v0,vE,vesc)*(F_DM(i*qunit,n))^2*Eprefactor_Si*wk*wk*FFarraySi(i,E_e))/(i*qunit*qunit*Eunit*2*2*wk);
end

dRdE_Si = prefac*sum(drde); %Sum over q-values


