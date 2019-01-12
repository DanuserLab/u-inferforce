function points = spleval(f)
%SPLEVAL Evaluation of a spline or spline curve.
%
% points = spleval(f)
%
% Computes points on the given spline or spline curve f between
% its extreme breaks.
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

% Original routine fnplt by C. de Boor / latest change: Oct. 25, 1997
% Simplified by Per Christian Hansen, IMM, 04/16/98.

% Set default number of points.
npoints = 300;

if (f.form(1)=='B'), f = sp2pp(f); end

[breaks,coefs,l,k,d] = ppbrk(f);
x = breaks(1) + (0:npoints)*((breaks(l+1)-breaks(1))/npoints);
v=ppual(f,x);

if (d==1), points=[x;v]; else points = v; end