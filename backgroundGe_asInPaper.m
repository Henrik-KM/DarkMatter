function backgroundGe_asInPaper = backgroundGe_asInPaper(E)
%This function obtains the PDF as proposed in arXiv:1402.7137. Parameters
%used are specific to Run 3 of CDMSlite. Uncertainties has been omitted.
%Output is in keV^(-1).

%Constants associated for Run 3 of CDMSlite.
A_divby_N   = [5.7, 15.2, 9.43, 18.7].*10^(-3); %[keV^(-1)]
mu          = [11.1, 1.26, 0.13, 0.04]*10^3;    %eV
sigma       = [108, 36.3, 13.9, 13];            %eV

%Analytic expression for the PDF. 
background = zeros(1,4);
for i = 1:4
    background(i) = A_divby_N(i)*0.5*(1+erf((E-mu(i))/(sqrt(2)*sigma(i))));
end

%Output. Convert keV^(-1) to eV^(-1).
backgroundGe_asInPaper = (0.005 + sum(background)).*10^(-3);