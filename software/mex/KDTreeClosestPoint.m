%KDTREECLOSESTPOINT for every query point in queryPts, find the closest point belonging to inPts
% 
% [idx, dist] = KDTreeClosestPoint(inPts,queryPts)
% 
% This function returns the index of the input point closest to each inPts.
% Supports 1D, 2D or 3D point sets.
%
% Input:
% 
%     inPts - an MxK matrix specifying the input points to test for distance
%     from the query points, where M is the number of points and K is the
%     dimensionality of the points.
% 
%     queryPts - an NxK matrix specifying the query points.
% 
% Output:
% 
%   idx - Nx1 array, the n-th element of which gives the index of
%   the input point closest to the the n-th query point.
% 
%   dist - Nx1 array, the n-th element of which gives the corresponding 
%   distance between the closest input point and the n-th query point.
%
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
