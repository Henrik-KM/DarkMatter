function dRdE_Ge = dRdE_Ge(m_chi,E_e,n,rho_c,MTarget,sigmatest,qunit,Eunit,qarray,ffarray)

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
wk = 2*alpha;                       %TBA
mu = (m_e*m_chi)/(m_e+m_chi);       %Reduced mass
sec2year = 60*60*24*365.25;         %conversion sec to year

MCell = 2*72.63*amu2kg;                                     %Mass of unit cell in Ge
bohr = 2.6817268*10^(-4);                                  %1 bohr in eV
VCell = 316.4269*bohr^3;                                   %Volume of cell in Bohr^3 converted to natural units
%Eprefactor_Ge = 2*pi^2/(alpha*m_e^2*VCell);                %Prefactor for calculation
Eprefactor_Ge = 1.8;

F_DM =@(q,n) (alpha*m_e/q)^n;       %Dark matter form factor definition
FFarrayGe = ffarray;                %Just assigns the formfactor array (Si)
q=qarray;                           %Just assigns the q values, unitless.


prefac =c^2*sec2year*rho_c*MTarget*sigmatest*alpha*m_e^2/(mu^2*MCell*m_chi);     %Prefactor 
%j=round(E_e/Eunit);
%Calculate all q-values in the energy bin.
for i=1:length(q)
    drde(i)=(Eunit*qunit*Eta(i*qunit,E_e*Eunit,m_chi)*(F_DM(i*qunit,n))^2*Eprefactor_Ge*wk*wk*FFarrayGe(i,E_e))/(i*qunit*qunit*Eunit*2*2*wk);
end

dRdE_Ge = prefac*sum(drde);


