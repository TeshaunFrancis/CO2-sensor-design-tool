%Original Solution by Rechnitz et al [Response Time Characteristics of the pCO2 Electrode]
% y = [H^+] mols/L
% x = pCO2 (mmHg)
%function [varout,y3,yquad,ylin,y4] = H_CO2(varargin)
function [varout,y3,y4] = H_CO2(varargin)
    varin = varargin{1};
    dye = varargin{2};
    nahco3 = varargin{3};
    if length(varargin) == 4
        pka = varargin{4};
    else
        pka = 7.7;
    end
    k0 = 4.47E-5;
    k1 = 2.6E-3;
    k2 = 1.72E-4;
    kT = k0*k1*k2;
    k3 = 5.59E-11;
    kw = 1E-14;
    ka = 10^(-pka);
    D1 = 0.2E-3;
    D2 = 0.1E-3;
    D = dye;
    %
    N1 = 1.4E-3;
    N2 = 4.4E-3;
    N3 = 8.4E-3;
    N = nahco3;
    %
    %syms x y
    %original = (y^3)-y*x*kT-2*x*kT*k3-y*kw+(y^2)*N == 0;
    %optochemical = (y^3)-y*x*kT-2*x*kT*k3-y*kw+(y^2)*N-(y^2)*(ka*D/(ka+y)) == 0;
    %opto_1 = -y*x*kT+(y^2)*N-(y^2)*(ka*D/(ka+y)) == 0;
    %opto_2 = (y^2)*N+y*N*ka-y*D*ka-y*x*kT-x*kT*ka == 0 ;
    % x = 1:0.1:40;
    x = varin;
    y = (sqrt((-D*ka + ka*N-kT*x).^2+4*ka*N*kT*x)+D*ka-ka*N+kT*x)/(2*N);%???? why squared?? oh cause quadractic formula
    y2 = sqrt((x*kT+D*ka)/N);
    y3 = 100*(y2-y)./y;
    varout = y;
    %% I need H2CO3, HCO3, and CO3 as well
    H2CO3 = k1*k0*x;
    HCO3 = H2CO3*k2./y;
    CO3 = HCO3*k3./y;
    y4 = [y,H2CO3,HCO3,CO3];
    
    % f= [f 200*k3*N/kT];
    % disp(f);
    % figure;
    % plot(x,y);
    % hold on
    % plot(x,y2);
    % yyaxis right
    % plot(x,y3);
    % hold off;
end