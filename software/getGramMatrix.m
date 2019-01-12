function [R,H]=getGramMatrix(forceMesh)
% this function obtains Gram matrix of the function hi(x) (basis function).
% Another output is Cholesky factorized matrix (R) of H.
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

H = zeros(2*forceMesh.numBasis); % both for x-direction and y-
% We make H for x-direction and copy later for y
xH = zeros(forceMesh.numBasis);
bItself = false;
bLateral = false;
bDiagonal = false;
gramItself=[];
gramLateral = [];
gramDiagonal = [];

% wtBar = waitbar(0,'Please wait, building Gram Matrix...');
for jj=1:forceMesh.numBasis
    plotClass=forceMesh.basis(jj).class;
    integrandx_jj = forceMesh.basisClass(plotClass).basisFunc(1).f_intp_x;% this is basically the same as integrandx_ii
    jjx = forceMesh.basis(jj).node(1);
    jjy = forceMesh.basis(jj).node(2);
    
    if ~bItself
        gramItself = discreteConvolution(forceMesh,jj,jj,0.1,integrandx_jj,jjx,jjy);
        bItself = true;
        xH(jj,jj) = gramItself;
    else
        xH(jj,jj) = gramItself;
    end
    neighBaseNoCand = forceMesh.nodeIDtoBaseNo(forceMesh.neigh(forceMesh.basis(jj).nodeID).cand);
    neighBaseNo = neighBaseNoCand(isfinite(neighBaseNoCand));
    for ii = neighBaseNo
        if forceMesh.basis(ii).node(1)==forceMesh.basis(jj).node(1) || forceMesh.basis(ii).node(2)==forceMesh.basis(jj).node(2)
            if ~bLateral
                gramLateral = discreteConvolution(forceMesh,jj,ii,0.1,integrandx_jj,jjx,jjy);
                bLateral = true;
                xH(jj,ii) = gramLateral;
            else
                xH(jj,ii) = gramLateral;
            end
        else
            if ~bDiagonal
                gramDiagonal = discreteConvolution(forceMesh,jj,ii,0.1,integrandx_jj,jjx,jjy);
                bDiagonal = true;
                xH(jj,ii) = gramDiagonal;
            else
                xH(jj,ii) = gramDiagonal;
            end
        end
    end
    %Below is very slow!
%     for ii=1:forceMesh.numBasis
%         % definition of field boudary (this can be just an overlapped region)
%         % center point of jj-th function
%         % if ii-th node is among the neighbors (IDs) of the jj-th node,
%         % calculate the integral between the two basis functions.
%         if any(forceMesh.basis(ii).nodeID == [forceMesh.neigh(forceMesh.basis(jj).nodeID).cand forceMesh.basis(jj).nodeID])
%             if forceMesh.basis(ii).node(1)==forceMesh.basis(jj).node(1) && forceMesh.basis(ii).node(2)==forceMesh.basis(jj).node(2)
%                 if ~bItself
%                     gramItself = discreteConvolution(forceMesh,jj,ii,1,integrandx_jj,jjx,jjy);
%                     bItself = true;
%                     xH(jj,ii) = gramItself;
%                 else
%                     xH(jj,ii) = gramItself;
%                 end
%             elseif forceMesh.basis(ii).node(1)==forceMesh.basis(jj).node(1) || forceMesh.basis(ii).node(2)==forceMesh.basis(jj).node(2)
%                 if ~bLateral
%                     gramLateral = discreteConvolution(forceMesh,jj,ii,1,integrandx_jj,jjx,jjy);
%                     bLateral = true;
%                     xH(jj,ii) = gramLateral;
%                 else
%                     xH(jj,ii) = gramLateral;
%                 end
%             else
%                 if ~bDiagonal
%                     gramDiagonal = discreteConvolution(forceMesh,jj,ii,1,integrandx_jj,jjx,jjy);
%                     bDiagonal = true;
%                     xH(jj,ii) = gramDiagonal;
%                 else
%                     xH(jj,ii) = gramDiagonal;
%                 end
%             end
% %            % get it back to the original coordinates
% %            figure, surf(xmat,ymat,zmat_jj)
% %            figure, surf(xmat,ymat,zmat_ii)
% %            integral2(integrandx_jj,xmin-jjx,xmax-jjx,ymin-jjy,ymax-jjy,'MaxFunEvals',10^10,'AbsTol',5e-10);
% %            % we can use integral2, but it takes much more time. Thus we
% %            stick with discrete method.
%         end
%     end
%     % Update the waitbar
%     if ishandle(wtBar) && mod(jj,10)==0
%         waitbar(jj/forceMesh.numBasis,wtBar,['Please wait, building Gram Matrix...' num2str(jj)]);
%     end
end
% close(wtBar)

xH=xH/max(xH(:)); % normalize it with max value.
H(1:forceMesh.numBasis,1:forceMesh.numBasis) = xH;
H(forceMesh.numBasis+1:2*forceMesh.numBasis,forceMesh.numBasis+1:2*forceMesh.numBasis) = xH;

%Cholesky factorization
xR = chol(xH);
R(forceMesh.numBasis+1:2*forceMesh.numBasis,forceMesh.numBasis+1:2*forceMesh.numBasis) = xR;
R(1:forceMesh.numBasis,1:forceMesh.numBasis) = xR;

function xH = discreteConvolution(forceMesh,jj,ii,refineFactor,integrandx_jj,jjx,jjy)
% get all 12 neighbor points from two nodes and get min and max
allNeighPos = [forceMesh.neigh(forceMesh.basis(jj).nodeID).pos; forceMesh.neigh(forceMesh.basis(ii).nodeID).pos];
xmin = min(allNeighPos(:,1));
xmax = max(allNeighPos(:,1));
ymin = min(allNeighPos(:,2));
ymax = max(allNeighPos(:,2));
% do the integral of jj-th function for debug - volume should
% be 66.667 since the shape is pyramid now.
iix = forceMesh.basis(ii).node(1);
iiy = forceMesh.basis(ii).node(2);
% refineFactor = 1;%1e-1; %1e-1 was too slow
[xmat,ymat] = meshgrid(xmin:refineFactor:xmax,ymin:refineFactor:ymax);
zmat_jj = integrandx_jj(xmat-jjx,ymat-jjy);
zmat_ii = integrandx_jj(xmat-iix,ymat-iiy);
% convolution in zmat_jj .* zmat_ii
xH = sum(sum(zmat_jj.*zmat_ii))/(refineFactor^2);


