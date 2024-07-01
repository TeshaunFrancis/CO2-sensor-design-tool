function varout = CO2_H(varargin)
    %Inverted function so that we can solve the equation easily
    %Original Solution by Rechnitz et al [Response Time Characteristics of the pCO2 Electrode]
    % y = [H^+] mols/L
    % x = pCO2 (mmHg)
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
    
    %( (y.^3)-y*(1E-14)+(y.^2)*N-(y.^2).*(10^(-7.7)*2.0000E-04./(10^(-7.7)+y)) )./(y.*1.99898400000000E-11+2*1.99898400000000E-11*5.59E-11);
    
    %opto_1 = -y*x*kT+(y^2)*N-(y^2)*(ka*D/(ka+y)) == 0;
    %opto_2 = (y^2)*N+y*N*ka-y*D*ka-y*x*kT-x*kT*ka == 0 ;    
    y = varin;
    x_num = (y.^3)-y*kw+(y.^2)*N-(y.^2).*(ka*D./(ka+y));
    x_den = y*kT+2*kT*k3;
    x = x_num ./ x_den;
    y_quad = (sqrt((-D*ka + ka*N-kT*x).^2+4*ka*N*kT*x)+D*ka-ka*N+kT*x)/(2*N);
    y_lin = (x*kT+D*ka)/N;  
    percentout = ((y.^3))./ x_den;
    percentout2 = ((y.^2)*N)./ x_den;
    percentout3 = (-(y.^2).*(ka*D./(ka+y)))./ x_den;
    percentout4 = ((y.^2)*N)./x_num;
    varout = [x;y_quad;y_lin;percentout;percentout2;percentout3;percentout4];
        
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

