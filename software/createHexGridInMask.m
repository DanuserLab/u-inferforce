function [xvec,yvec]=createHexGridInMask(spacing,bgdMask)
% createHexGridInMask creates grid points in bgdMask with spacing (in
% pixel).
% Sangyoon Han June 2013
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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

% First get the corners (here we take field no.1 but that is arbitrary):
x1D = find(any(bgdMask,1));
y1D = find(any(bgdMask,2));


LeftUpperCorner(1:2) = [x1D(1), y1D(1)];
RightLowerCorner(1:2) = [x1D(end), y1D(end)];

oldSize=RightLowerCorner-LeftUpperCorner;

% Use given spacing for edgeLength
edgeLength=spacing+mod(spacing,2); % make it a multiple of 2

% then the hight of the triangles is:
triHight=ceil(edgeLength*sqrt(3)/2);

% The hexagonal grid is the superposition of two regular grids with
% yspacing=edgeLength and xspacing=2*triHight with a relative
% yshift=edgeLength/2 and xshift=triHight.

% Grid1:
xval1=LeftUpperCorner(1):(2*triHight):RightLowerCorner(1);
yval1=LeftUpperCorner(2):edgeLength:RightLowerCorner(2);

% the second starts at the shifted origin:
xval2=(LeftUpperCorner(1)+triHight):(2*triHight):RightLowerCorner(1);
yval2=(LeftUpperCorner(2)+edgeLength/2):edgeLength:RightLowerCorner(2);

% Now place the new grids in the center of the old bounds:
% first get the new size:
xmax=max([xval1 xval2]);
ymax=max([yval1 yval2]);

% The left upper corner is still the same by definition:
newSize =[xmax-LeftUpperCorner(1)    ymax-LeftUpperCorner(2)];

% calculate the size difference and the needed shift:
sizeDiff=oldSize-newSize;
sizeDiff=sizeDiff-mod(sizeDiff,2); % make it a multiple of 2

if sizeDiff(1)<0 || sizeDiff(2)<0
    display('something went wrong! Nothing has been done');
    return;
end

% shift all grid points:
xval1=xval1+sizeDiff(1)/2;
xval2=xval2+sizeDiff(1)/2;

yval1=yval1+sizeDiff(2)/2;
yval2=yval2+sizeDiff(2)/2;

% Now we have everything to create the two regular grids:
[xgrid1,ygrid1] = meshgrid(xval1, yval1);
xvec1=reshape(xgrid1,[],1);
yvec1=reshape(ygrid1,[],1);

[xgrid2,ygrid2] = meshgrid(xval2, yval2);
xvec2=reshape(xgrid2,[],1);
yvec2=reshape(ygrid2,[],1);

% Now fuse the two vectors:
xvec=vertcat(xvec1,xvec2);
yvec=vertcat(yvec1,yvec2);

bothVec=horzcat(xvec,yvec);
% These two vectors should be sorted for later purposis:
bothVec=sortrows(bothVec,[1 2]);
xvec=bothVec(:,1);
yvec=bothVec(:,2);

% Mask points with bgdMask
Npoints = length(bothVec);
outsideIdx = false(Npoints,1);

for ii=1:Npoints
    if bgdMask(round(bothVec(ii,2)),round(bothVec(ii,1)))
        outsideIdx(ii) = true;
    end
end
xvec=bothVec(outsideIdx,1);
yvec=bothVec(outsideIdx,2);
