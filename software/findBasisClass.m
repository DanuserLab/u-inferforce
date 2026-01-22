function [foundClass]=findBasisClass(basisClass,currNeighPos,currCtrNodePos)
if nargin<3
    currCtrNodePos=[0 0];
end

foundClass=[];

% translate all positions to the center [0.0]:
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
numCurrNeigh=size(currNeighPos,1);
pts2=currNeighPos-repmat(currCtrNodePos,numCurrNeigh,1);
pts2=sortrows(pts2);

numClass=length(basisClass);
% run through all classes and see if you find a match:
for classNo=1:numClass
    % first check if the number of neighbors is consistent:
    check1=(basisClass(classNo).numNeigh==numCurrNeigh);
    
    % then check if the neighbors are consistent:
    if check1
        % this can be improved! By refining not all identical
        % classes are found since the ordering might be not exactly
        % be the same!
        % check2=sum(sum(abs(basisClass(classNo).neighPos-(myMesh.neigh(j).pos-repmat(myMesh.p(j,:),length(myMesh.neigh(j).cand),1)))))<10^4*eps;
        
        pts1=basisClass(classNo).neighPos;        
        
        % sort these points (pts2 have been sorted above):
        pts1=sortrows(pts1);        
        
        % calculate the difference between the points
        check2=sum(sum(abs(pts1-pts2)))<10^4*eps;
    end
    
    % if both checks are positive:
    if check1 && check2
        foundClass=vertcat(foundClass,classNo);
    end
end
