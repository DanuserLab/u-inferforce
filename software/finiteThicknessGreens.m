function [G]=finiteThicknessGreens(i,j,x,y,E,h)
% finiteThicknessGreens calculates Green's function only in tangential direction
% This assumes the thickness (h) in gel
% Adapted from Merkel et al. Cell Force Microscopy on Elastic Layers of Finite
% Thickness, Biophys J, 2007
% This Greensfunction is only valid for v=0.5;
% G = | A1-(x^2-y^2)/r^2*A2     -2*x*y/r^2*A2               |
%        | -2*x*y/r^2*A2                A1+(x^2-y^2)/r^2*A2   |
% where A1 = A1B * (.12*exp(-0.43*h/r)+.88*exp(-0.83*h/r))
%             A2 = A2B * (1+1.22*h/r+1.31*(h/r)^2.23)*exp(-1.25*h/r)
%             A1B = (1+v)*(2-v)/(2*pi*E*r)
%             A2B = - (1+v)*v/(2*pi*E*r)
% Sangyoon Han April 2013
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
%
% This file is part of TFM_Package.
% 
% TFM_Package is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TFM_Package is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with TFM_Package.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
v=0.5;

r=sqrt(x.^2+y.^2);
s=r/h;
% A1B=(1+v).*(2-v)./(2.*pi.*E.*r);
% A1 = A1B .* (.12.*exp(-0.43.*s)+.88.*exp(-0.83.*s));
A1 = (1+v).*(2-v)./(2.*pi.*E.*r) .* (.12.*exp(-0.43.*s)+.88.*exp(-0.83.*s));
% A2B = - (1+v)*v./(2*pi*E*r);
% A2 = A2B .* (1+1.22.*s+1.31.*(s).^2.23)*exp(-1.25.*s);
A2 = - (1+v).*v./(2.*pi.*E.*r) .* (1+1.22.*s+1.31.*(s).^2.23).*exp(-1.25.*s);
if i==1 && j==1 
    G = A1-(x.^2-y.^2)./(r.^2).*A2;
elseif (i==1 && j==2) || (i==2 && j==1) 
    G = -2.*x.*y./(r.^2).*A2;
elseif i==2 && j==2 
    G = A1+(x.^2-y.^2)./(r.^2).*A2;
else
    display('something went wrong')
end

%remove the NaN if the Greensfunction has been evaluated at zero.
nanMat=isnan(G);
if sum(nanMat(:))>0
    %display('The Boussinesq Greensfunction has been evaluated at zero')
    %display('To resolve this, the maximal value has been set')
    %G(nanMat)=max(max(G));
    
    %display('To resolve this, this value has been set to zero')
    G(nanMat)=0;
end
%remove the inf if the Greensfunction has been evaluated at zero.
infMat=isinf(G);
if sum(infMat(:))>0
    G(infMat)=max(max(G(~infMat)))*100;
    %display('To resolve this, this value has been set to zero')
%     G(infMat)=0;
end