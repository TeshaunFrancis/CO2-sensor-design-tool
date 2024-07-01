function [A,HA] = A_H(varargin)
%function [y1,y2,y3,y4,y5,y6,y7,y8] = A_H(varin,D,L,pka)
if length(varargin) > 4
    eA = varargin{5};
    eHA = varargin{6};
else
    eA = 24000;
    eHA = 20000;
end
H = varargin{1};
D = varargin{2};
L = varargin{3};
pka = varargin{4};

ka = 10.^-pka;
%This is Absolute A-
A = L * D * eA * (ka ./ (ka + H)); %just pH sensitive portion 
%This is at 405 
HA = L * D * eHA * (H ./ (ka + H)) + L * D * 10000 * (ka ./ (ka + H)); %just pH sensitive portion 

%y8 = ((-L*D*e454*ka)-(varin-L*D*e415)*ka) ./ (varin - L * D * e415);
end

