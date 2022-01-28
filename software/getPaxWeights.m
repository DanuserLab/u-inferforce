function paxWeights = getPaxWeights(forceMesh,paxImage,x_vec,y_vec,ux_vec,uy_vec)
% this function calculates weight for each force nodes according to
% paxImage. The purpose is to get appropriate weight at a position where
% paxillin signal is so that position with no cell (or paxillin) will have
% more penalty ( minimum should be 1 and maximum is subjective. I'll say
% it's 100 for now. In case the paxillin signal is in between the two
% nodes, the distance to the signal will proportionally reduce the effect.
% Again, 1 for highest paxillin signal, 100 for no signal
%
% Copyright (C) 2022, Danuser Lab - UTSouthwestern 
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

% each node postition should be examined at its deformed position
paxWeights = zeros(2*forceMesh.numBasis); % both for x-direction and y-
% [grid_mat,iu_mat, ~,~] = interp_vec2grid([x_vec y_vec], [ux_vec uy_vec],[],forceMesh);
xH = zeros(forceMesh.numBasis);
wtBar = waitbar(0,'Please wait, building weight matrix from paxillin channel');
paxImage = double(paxImage);
pId = paxImage/max(paxImage(:));
alpha = graythresh(pId);
bwPI = im2bw(pId,alpha);
bwPI2 = bwmorph(bwPI,'clean');

for jj=1:forceMesh.numBasis
    plotClass=forceMesh.basis(jj).class;
    integrandx_jj = forceMesh.basisClass(plotClass).basisFunc(1).f_intp_x;% this is basically the same as integrandx_ii
    jjx = forceMesh.basis(jj).node(1);
    jjy = forceMesh.basis(jj).node(2);
    allNeighPos = forceMesh.neigh(forceMesh.basis(jj).nodeID).pos;
    xmin = min(allNeighPos(:,1));
    xmax = max(allNeighPos(:,1));
    ymin = min(allNeighPos(:,2));
    ymax = max(allNeighPos(:,2));
    % boundary area of the basis node in the form of mask;
    [xmat,ymat] = meshgrid(xmin:xmax,ymin:ymax);
    zmat_jj = integrandx_jj(xmat-jjx,ymat-jjy);

    % crop paxImage with basisClass area considering displacement at the
    % location closest to the force node
    distToNode = sqrt(sum(([x_vec y_vec]- ...
                            ones(length(x_vec),1)*[jjx jjy]).^2,2));
                        
    [minD,closest_ind] = min(distToNode);
    dispVec = [ux_vec(closest_ind) uy_vec(closest_ind)];
    
    meshLx = size(zmat_jj,2);
    meshLy = size(zmat_jj,1);
    
    paxCropped = imcrop(bwPI2,[jjx+dispVec(1)-floor(meshLx/2)+0.5 jjy+dispVec(2)-floor(meshLy/2)+0.5 meshLx-1 meshLy-1]); %it accounts for deformation
    
    paxDistProd = zmat_jj.*paxCropped;
    xH(jj,jj) = sum(paxDistProd(:));
    if ishandle(wtBar) && mod(jj,10)==0
        waitbar(jj/forceMesh.numBasis,wtBar,'Please wait, building weight matrix with paxillin channel');
    end
end
xH=xH/max(xH(:)); % normalize it with max value.
xH = diag(diag(eye(size(xH)))./diag(xH));%xH\eye(size(xH)); % now node with maximum paxillin signal is 1 while node with small signal will become very large or not very much.
maxFactor=100;
xH(isinf(xH))=maxFactor;
% xH = xH*(maxFactor-1)/(max(xH(:))-1)+(max(xH(:))-maxFactor)/(max(xH(:))-1);
paxWeights(1:forceMesh.numBasis,1:forceMesh.numBasis) = xH;
paxWeights(forceMesh.numBasis+1:2*forceMesh.numBasis,forceMesh.numBasis+1:2*forceMesh.numBasis) = xH;
