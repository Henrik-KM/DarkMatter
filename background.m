%%This function gives the background of silicon or germanium detectors in
%the form of a vector with the expected number of events for 1 kg year
%exposure summed over the energy intervals specified by the input vector Q.


function bg=background(E_e,Eunit,Q,runSi)


if runSi
    
    btmp = 0.073*ones(1,500); %Background signal, Si.
    
else 
    btmp = getBackgroundGe; %Background signal, Ge.
end

bg=zeros(1,max(Q));     %Preallocating background vector
%Sum over energy intervals corresponding to different values in Q
for i=1:max(Q)
    bg(i)=sum(btmp.*(Q==i));
end