function [insideIdx] = maskVectors(dispMatX,dispMatY,bwstackImg)
Npoints = length(dispMatX);
insideIdx = false(Npoints,1);

for ii=1:Npoints
    if round(dispMatY(ii)) <= size(bwstackImg,1) && ...
            round(dispMatY(ii)) > 0 && ...
            round(dispMatX(ii)) <= size(bwstackImg,2) && ...
            round(dispMatX(ii)) > 0 && ...
            bwstackImg(round(dispMatY(ii)),round(dispMatX(ii)))
        insideIdx(ii) = true;
    end
end
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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
