function [xvec,yvec]=createHexGridFromDisplField(displField,forceInteger)
% For Benedikt's software to work, the displacement field has to be
% interpolated on a rectangular grid, with an even number of grid points
% along each edge. Furthermore, one has to make sure that the noisy data 
% has not to be extrapolated. This may happen along the edges. To prevent
% this, extract the corner of the displacement grid, calculate how often 
% (even number) the optimal gridspacing fits into each dimension, then 
% place the regular grid centered to the orignal bounds. Thereby make sure 
% that the edges have been eroded to a certain extend.
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

if nargin<2 || isempty(forceInteger)
    forceInteger=0;
end

% First get the corners (here we take field no.1 but that is arbitrary):
LeftUpperCorner(1:2) = [min(displField(1).pos(:,1)), min(displField(1).pos(:,2))];
RightLowerCorner(1:2) = [max(displField(1).pos(:,1)), max(displField(1).pos(:,2))];
oldSize=RightLowerCorner-LeftUpperCorner;

% Next get the average bead density (here we take field no.1 but that is arbitrary):
numBeads=length(displField(1).pos);
areaImg=prod(RightLowerCorner-LeftUpperCorner);
interBeadDensity=numBeads/areaImg;

% from that calculate the edge length within the hexagonal grid (the edge
% length of the regular triangles.):
if forceInteger==1
    edgeLength=ceil(sqrt(2/(sqrt(3)*interBeadDensity)));
    edgeLength=edgeLength+mod(edgeLength,2); % make it a multiple of 2
else
    edgeLength=sqrt(2/(sqrt(3)*interBeadDensity));
end

% then the hight of the triangles is:
if forceInteger==1
    triHight=ceil(edgeLength*sqrt(3)/2);
else
    triHight=edgeLength*sqrt(3)/2;
end

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
if forceInteger==1
    sizeDiff=oldSize-newSize;
    sizeDiff=sizeDiff-mod(sizeDiff,2); % make it a multiple of 2
else
    sizeDiff=oldSize-newSize;
end

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
%it might be needed to exchange x and y at the end;