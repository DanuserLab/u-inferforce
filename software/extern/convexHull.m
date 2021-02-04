function [hull, inds] = convexHull(points)
%CONVEXHULL Convex hull of a set of points
%
%   POLY = convexHull(POINTS)
%   Computes the convex hull of the set of points POINTS. This function is
%   mainly a wrapper to the convhull function, that format the result to a
%   polygon.
%
%   [POLY INDS] = convexHull(POINTS)
%   Also returns the indices of convex hull vertices within the original
%   array of points.
%
%   Example
%     % Draws the convex hull of a set of random points
%     pts = rand(30,2);
%     drawPoint(pts, '.');
%     hull = convexHull(pts);
%     hold on; 
%     drawPolygon(hull);
%
%     % Draws the convex hull of a paper hen
%     x = [0 10 20  0 -10 -20 -10 -10  0];
%     y = [0  0 10 10  20  10  10  0 -10];
%     poly = [x' y'];
%     hull = convexHull(poly);
%     figure; drawPolygon(poly);
%     hold on; axis equal;
%     drawPolygon(hull, 'm');
%
%   See also
%   polygons2d, convhull
%
%
% Copyright (C) 2021, Danuser Lab - UTSouthwestern 
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

% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-04-08,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.

% checkup on array size
if size(points, 1) < 3
    hull = points;
    inds = 1:size(points, 1);
    return;
end

% compute convex hull by calling the 'convhull' function
inds = convhull(points(:,1), points(:,2));
hull = points(inds, :);
