function [force]=assumedForce(j,x,y)
xshift=7;
yshift=7;
if j==1
    force=x*0;
else
    force=(heaviside(1-((x-xshift).^2+(y-yshift).^2)).*(exp(-((x-xshift).^2+(y-yshift).^2))-1/exp(1)));
end

% if j==1
%     force=x*0;
% else
%     force=(heaviside(1-(x.^2+y.^2)).*(exp(-(x.^2+y.^2))-1/exp(1)));
% end
%
% Copyright (C) 2026, Danuser Lab - UTSouthwestern 
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
