% Function to simulate greenhouse lettuce
% x(1)    = x.dw
% x(2)    = x.co2In
% x(3)    = x.tempIn
% x(4)    = x.vapIn

% d(1)    = d.rad;
% d(2)    = d.co2Out
% d(3)    = d.tempOut
% d(4)    = d.vapOut

% u(1)    = u.co2
% u(2)    = u.vent
% u(3)    = u.heat

function dxdt = GHLettuce_ODE(t,x,u,d,p)
p       = DefineParameters;
dxdt = [p.alfaBeta*(...
    (1-exp(-p.laiW*x(1)))*p.photI0.*d(1).* ...
    (                  -p.photCO2_1*x(3).^2+p.photCO2_2*x(3)-p.photCO2_3).*(x(2)-p.photGamma) ...
    ./(p.photI0*d(1)+(-p.photCO2_1*x(3).^2+p.photCO2_2*x(3)-p.photCO2_3).*(x(2)-p.photGamma)))...
    - p.Wc_a*x(1).*2.^(0.1*x(3)-2.5);
    
    1/p.CO2cap*(-((1-exp(-p.laiW*x(1)))*p.photI0.*d(1).* ...
    (                  -p.photCO2_1*x(3).^2+p.photCO2_2*x(3)-p.photCO2_3).*(x(2)-p.photGamma) ...
    ./(p.photI0*d(1)+(-p.photCO2_1*x(3).^2+p.photCO2_2*x(3)-p.photCO2_3).*(x(2)-p.photGamma)))...
    + p.CO2c_a.*x(1).*2.^(0.1*x(3)-2.5) + u(1)/1e6 - (u(2)/1e3+p.leak).*(x(2)-d(2)));
    
    1/p.aCap*(u(3) - (p.ventCap*u(2)/1e3+p.trans_g_o).*(x(3)-d(3)) + p.rad_o_g*d(1));
    
    1/p.H2Ocap*((1-exp(-p.laiW*x(1)))*p.evap_c_a.*(p.satH2O1./(p.R*(x(3)+p.T)).*...
    exp(p.satH2O2*x(3)./(x(3)+p.satH2O3))-x(4)) - (u(2)/1e3+p.leak).*(x(4)-d(4)))];
end

function y = GHLettuce_ODE_conv(x,u,d,p,h)   %units problem
     y = [1*x(1);
        1e-3*co2dens2ppm(x(3),x(2));
        x(3);
        vaporDens2rh(x(3), x(4))];
end
