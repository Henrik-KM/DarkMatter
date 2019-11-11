function getBackgroundGe = getBackgroundGe
%This function builds the background signal for Germanium. It renormalizes
%the obtained PDF to an energy range 0.1-50 eV. Output is a vector ranging
%from 0.1:50 eV in steps of 0.1 eVin terms of #/eV/kg/year.

f_bg=zeros(1,20*10^4); %Steps of 0.1 eV

E_analysis = linspace(1,20*10^4,20*10^4); %Analysis range. 

Eunit = 0.1; % Unit is 0.1 eV.

%Compute PDF as given by paper.
for i = 1:length(E_analysis)
    f_bg(i) = backgroundGe_asInPaper(E_analysis(i).*Eunit).*Eunit;
end

%normalized PDF?
norm = trapz(f_bg); %Sanity check, 0.9962 is good enough.

%Renormalization process. f_distr is the new PDF. 
A=trapz(f_bg)./trapz(f_bg(1:500)); 
f_distr = A*f_bg;

%Obtain final eventrate in /eV/kg/year
dRbdEGe = (0.1.*365.*f_distr(1:500))./trapz(f_distr(1:500));

getBackgroundGe = dRbdEGe;
