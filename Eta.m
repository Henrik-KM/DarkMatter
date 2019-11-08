function Eta = Eta(q, E_e, m_chi,v0,vE,vesc)
% Computes the velocity function used in computing dRdE.
% Modify v0, vE, vesc to modify the behaviour of this function.

c = 2.99792458*1e10;                                                    %speed of light [cm/s]
% v0 = 230*1e5;                                                           %typical velocity [cm/s]
% vE = 240*1e5;                                                           %velocity earth [cm/s]
% vesc = 600*1e5;                                                         %escape velocity earth
KK = v0^3*(-2*exp(-vesc^2/v0^2)*pi*(vesc/v0)+pi^(3/2)*erf(vesc/v0));    %normalization constant [cm/s]^(-3)

v_min =@(q,E_e,m_chi) c*(q/(2*m_chi)+(E_e/q));                          %Vmin velocity function

prefact = (v0^2*pi)/(2*vE*KK);                                          %Prefactor

if v_min(q,E_e,m_chi) < (vesc - vE)
    Eta = prefact*(-4*exp(-vesc^2/v0^2)*vE+sqrt(pi)*v0*(erf((v_min(q,E_e,m_chi)+vE)/v0)-erf((v_min(q,E_e,m_chi)-vE)/v0)));
else
    if v_min(q,E_e,m_chi) < (vesc + vE) && (v_min(q,E_e,m_chi)) > (vesc - vE)
        Eta = prefact*(-2*exp(-vesc^2/v0^2)*(vE+vesc-v_min(q,E_e,m_chi))+sqrt(pi)*v0*(erf(vesc/v0)-erf((v_min(q,E_e,m_chi)-vE)/v0)));
    else
        if v_min(q,E_e,m_chi) > (vesc + vE)
        Eta = 0;
        end
    end
end