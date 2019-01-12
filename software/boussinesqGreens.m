%This Greensfunction is only valid for v=0.5 if not specified;
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
function [G]=boussinesqGreens(i,j,x,y,E,v)
if nargin <6
    v=0.5;
end

r=sqrt(x.^2+y.^2);
preFactor=(1+v)./(pi.*E.*r.^3);

if i==1 && j==1 
    G=preFactor.*((1-v).*r.^2+v.*x.^2);
elseif (i==1 && j==2) || (i==2 && j==1) 
    G=v*preFactor.*x.*y;
elseif i==2 && j==2 
    G=preFactor.*((1-v).*r.^2+v.*y.^2);
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