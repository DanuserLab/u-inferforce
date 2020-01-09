function [M]=calcFwdMap(x_vec_u, y_vec_u, forceMesh, E,span,meshPtsFwdSol,method)
if nargin<5
    baseSpan=1:forceMesh.numNodes;
    method = 'fast';
else
    baseSpan=span;
end

if nargin<6
    meshPtsFwdSol=[];
    method = 'fast';
end

Bounds = forceMesh.bounds;
numNodes = forceMesh.numNodes;
ux=zeros(length(x_vec_u),2*numNodes);
uy=zeros(length(y_vec_u),2*numNodes);
ux2ndhalf=zeros(length(x_vec_u),numNodes);
uy2ndhalf=zeros(length(y_vec_u),numNodes);

if isempty(gcp('nocreate'))
    try
        parpool local;
    catch
        matlabpool;
    end
end

parfor j=baseSpan
    display([num2str(j),' of: ',num2str(length(baseSpan))]);
    %limits for integration: integrate base function only over their
    %respective support:
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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
    
    xmin=Bounds(j).x(1);
    xmax=Bounds(j).x(2);
    ymin=Bounds(j).y(1);
    ymax=Bounds(j).y(2);

    if strcmp(method,'fast')
        [ux(:,j),uy(:,j)]=fwdSolution(x_vec_u,y_vec_u,E,xmin,xmax,ymin,ymax,forceMesh.base(j).f_intp_x,forceMesh.base(j).f_intp_y,'fft',[],meshPtsFwdSol);
        [ux2ndhalf(:,j),uy2ndhalf(:,j)]=fwdSolution(x_vec_u,y_vec_u,E,xmin,xmax,ymin,ymax,forceMesh.base(j+forceMesh.numNodes).f_intp_x,forceMesh.base(j+forceMesh.numNodes).f_intp_y,'fft',[],meshPtsFwdSol);
    elseif strcmp(method,'conv_free')
        [ux(:,j),uy(:,j)]=fwdSolution(x_vec_u,y_vec_u,E,xmin,xmax,ymin,ymax,forceMesh.base(j).f_intp_x,forceMesh.base(j).f_intp_y,'conv_free',[],meshPtsFwdSol);
        [ux2ndhalf(:,j),uy2ndhalf(:,j)]=fwdSolution(x_vec_u,y_vec_u,E,xmin,xmax,ymin,ymax,forceMesh.base(j+forceMesh.numNodes).f_intp_x,forceMesh.base(j+forceMesh.numNodes).f_intp_y,'conv_free',[],meshPtsFwdSol);
    end
end
ux(:,numNodes+1:2*numNodes) = ux2ndhalf;
uy(:,numNodes+1:2*numNodes) = uy2ndhalf;
M=vertcat(ux,uy);

% plot an example to see if it works correctly
ind=10;
if forceMesh.numNodes>ind-1
    xmin=min(x_vec_u);
    ymin=min(y_vec_u);
    xmax=max(x_vec_u);
    ymax=max(y_vec_u);
    figure(11)
    quiver(x_vec_u,y_vec_u,ux(:,ind),uy(:,ind))
    hold on
    quiver(x_vec_u,y_vec_u,ux(:,ind+forceMesh.numNodes),uy(:,ind+forceMesh.numNodes))
    xlim([xmin xmax])
    ylim([ymin ymax])
    hold off
end