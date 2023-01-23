function [myMesh]=createMeshAndBasisFastBEM(x_vec,y_vec,keepBDPts,basisClassInit,doPlot)
% to make it very efficient, x_vec and y_vec should be such that they yield
% a hexagonal lattice.
%
% Copyright (C) 2023, Danuser Lab - UTSouthwestern 
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

if nargin < 3 || isempty(keepBDPts)
    keepBDPts=true;
end

if nargin < 4 || isempty(basisClassInit)
    basisClass=[];
else
    % this is needed only for refining the mesh!
    basisClass=basisClassInit;
end

if nargin < 5 || isempty(doPlot)
    doPlot=0;
end

% first check if the points represent a square lattice:
[p,t,dt,x_vec,y_vec,isSq]=delaunaySqLattice(x_vec,y_vec);

if ~isSq
    %delaunay_mesh=delaunay(x_vec,y_vec);
    dt=DelaunayTri(x_vec,y_vec);
    % dt has two fields:
    % dt.X             : contains the position of the nodes
    % dt.Triangulation : contains triangles identified by three node numbers.
    %                    The number, is the position/row of points in dt.X

    % we read this out because the access of the Delaunay-structure takes so
    % long:
    % bdPtsID = convexHull(dt);
    p=dt.X;
    t=dt.Triangulation;
end
% find the neighbors to each node, the bounds of each basis function and
% the convex hull of the grid:
[neigh,bounds,bdPtsID]=findNeighAndBds(p,t);


myMesh.dt=dt;  % DelaunayTri structure
myMesh.p=p;    % this is the same as dt.X
myMesh.neigh=neigh;
myMesh.bounds=bounds;
myMesh.numNodes=length(myMesh.p(:,1));

%create the basis functions and interpolate them using the Delaunay Triangulation:
zeroNodes=[];
baseNo=0;
for j=1:myMesh.numNodes
    % classify the basis function.
    % check if the current basis function falls into an existing basis class:
    noClass=false;
    foundClass=[];
    if isempty(basisClass)
        noClass=true;
        newClNo=1;
    else
        currNeighPos   = myMesh.neigh(j).pos;
        currCtrNodePos = myMesh.p(j,:);
        [foundClass]=findBasisClass(basisClass,currNeighPos,currCtrNodePos);        
        if isempty(foundClass)
            noClass=true;
            numClass=length(basisClass);
            newClNo=numClass+1;
        elseif length(foundClass)>1
            error('Basis Classes should be unique, something went wrong!')
        end
    end
    
    if isempty(foundClass)
        matchClass=newClNo;
    else
        matchClass=foundClass;
    end
    
    bdCriteriaOK=true;
    % if bdPts can be neglected
    if ~keepBDPts
        % test if the current pt is part of the boundary. If yes, the
        % criteria that the pts is not on the bd is not fulfilled.
        bdCriteriaOK = ~ismember(j, bdPtsID);         
    end
    if noClass && bdCriteriaOK
        % sort in the position of the neighboring nodes. The center node is
        % positioned at the origin.
        
        % Each neighbor around the center node yields a triangle:
        basisClass(newClNo).centerPos=[0, 0];
        basisClass(newClNo).numNeigh =length(myMesh.neigh(j).cand);
        basisClass(newClNo).neighPos =myMesh.neigh(j).pos-repmat(myMesh.p(j,:),basisClass(newClNo).numNeigh,1);
        
        % Make its own Delaunay-Triangulation, that is super fast:
        % Put all the points together:
        allPts   = vertcat(vertcat(basisClass(newClNo).centerPos,basisClass(newClNo).neighPos));
        dtBaseSup= DelaunayTri(allPts(:,1),allPts(:,2));
        
        % The two basis functions are one at the center, zero otherwise.
        % They are a combination of the two column vectors:
        % all zero:
        f_disc_0  =zeros(length(allPts),1);
        
        % first component is 1:
        f_disc_1   =f_disc_0;
        f_disc_1(1)=1;              
        
        % There might be nans if the points where the function is evaluated
        % is not within the support.        
        basisClass(newClNo).basisFunc(1).f_intp_x= @(x,y) nan2zeroTriScatteredInterp(x,y,dtBaseSup,f_disc_1,'linear');
        basisClass(newClNo).basisFunc(1).f_intp_y= @(x,y) nan2zeroTriScatteredInterp(x,y,dtBaseSup,f_disc_0,'linear'); % only zeros
        
        basisClass(newClNo).basisFunc(2).f_intp_x= @(x,y) nan2zeroTriScatteredInterp(x,y,dtBaseSup,f_disc_0,'linear'); % only zeros
        basisClass(newClNo).basisFunc(2).f_intp_y= @(x,y) nan2zeroTriScatteredInterp(x,y,dtBaseSup,f_disc_1,'linear'); 
        
        basisClass(newClNo).dtBaseSup=dtBaseSup;
        
        % calculate the area of the support. allPts, not only includes
        % neighbors but also the center which might be important for basis
        % functions at the boundary, i.e. here, the center may be part of
        % the convex hull!!!
        allPts=vertcat(basisClass(newClNo).centerPos,basisClass(newClNo).neighPos);
        [~, areaSupport] = convhull(allPts(:,1),allPts(:,2));
        basisClass(newClNo).area=areaSupport;
        
        % calculate the volume of the basis-function. This is needed later
        % on to calculate the forces from stresses in an efficient way.
        % Create a pyramid of height 1 over the support of the base
        % function:
        % First fill in the polygone-points of the support, they have
        % height zero.
        % The x,y-components are the points of the basis function:
        X=allPts;
        % The z-components:
        X(:,3)=zeros( size(allPts,1),1);        
        % now fill in the unit peak with height 1 above the center node (0,0) 
        X(end+1,:)=[0 0 1];
        
        % calculate the convex hull:
        [~, unitVolume] = convhulln(X);
        basisClass(newClNo).unitVolume=unitVolume;
    end
    
    % take it only as a base function if it has enough connectivity, else
    % sort it into the nodes where forces are set to zero. These nodes are
    % not considered in the force reconstruction:
    if bdCriteriaOK
        baseNo=baseNo+1;
        % we have to tread each as two basis function later on!
        basis(baseNo).class = matchClass;
        basis(baseNo).node  = myMesh.p(j,:);
        basis(baseNo).nodeID= j;
        basis(baseNo).numNeigh   = length(myMesh.neigh(j).cand);
        basis(baseNo).area       =  basisClass(basis(baseNo).class).area;
        basis(baseNo).unitVolume =  basisClass(basis(baseNo).class).unitVolume;  
        myMesh.nodeIDtoBaseNo(j) = baseNo;
    else
        myMesh.nodeIDtoBaseNo(j) = NaN;
        zeroNodes = vertcat(zeroNodes,myMesh.p(j,:));
    end
end

myMesh.basisClass = basisClass;
myMesh.basis      = basis;
myMesh.numBasis   = length(myMesh.basis);
myMesh.zeroNodes  = zeroNodes;
myMesh.keepBDPts  = keepBDPts;

if doPlot==1
    marker=['or','ob','og','ok','oy','oc','om','sr','sb','sg','sk','sy','sc','sm','dr','db','dg','dk','dy','dc','dm','*r','*b','*g','*k','*y','*c','*m'];
    % plot the basisfunction classes:
    figure(1)
    triplot(myMesh.dt.Triangulation,myMesh.dt.X(:,1),myMesh.dt.X(:,2))
    hold on;
    for j=1:length(myMesh.basis)
        plot(myMesh.basis(j).node(1),myMesh.basis(j).node(2),marker(2*myMesh.basis(j).class-1:2*myMesh.basis(j).class))
    end
    hold off;
    
    % plot an example to see if it works correctly
    ind=350;
    if length(x_vec)>ind-1
        xmin=min(x_vec);
        ymin=min(y_vec);
        xmax=max(x_vec);
        ymax=max(y_vec);
        
        dx=xmax-xmin;
        dy=ymax-ymin;
        
        
        pointsPerEdge=round(sqrt(length(x_vec)));
        [x_fine y_fine]=meshgrid(linspace(-dx,dx,20*pointsPerEdge) , linspace(-dy,dy,20*pointsPerEdge));
        
        figure(50)
        % triplot(myMesh.dt
        triplot(myMesh.dt.Triangulation,myMesh.dt.X(:,1),myMesh.dt.X(:,2))
        hold on
        ID=myMesh.basis(ind).nodeID;
        plot(myMesh.p(myMesh.neigh(ID).cand,1),myMesh.p(myMesh.neigh(ID).cand,2),'or')
        plot(myMesh.p(ID,1),myMesh.p(ID,2),'ob')
        plot([myMesh.bounds(ID).x(1) myMesh.bounds(ID).x(1) myMesh.bounds(ID).x(2) myMesh.bounds(ID).x(2) myMesh.bounds(ID).x(1)],[myMesh.bounds(ID).y(1) myMesh.bounds(ID).y(2) myMesh.bounds(ID).y(2) myMesh.bounds(ID).y(1) myMesh.bounds(ID).y(1)],'k')
        
        % the class to be plotted is:
        plotClass=myMesh.basis(ind).class;
        % the shift is:
        plotShift=myMesh.basis(ind).node;
        
        quiver(x_fine+plotShift(1),y_fine+plotShift(2),myMesh.basisClass(plotClass).basisFunc(1).f_intp_x(x_fine,y_fine),myMesh.basisClass(plotClass).basisFunc(1).f_intp_y(x_fine,y_fine),'b')
        quiver(x_fine+plotShift(1),y_fine+plotShift(2),myMesh.basisClass(plotClass).basisFunc(2).f_intp_x(x_fine,y_fine),myMesh.basisClass(plotClass).basisFunc(2).f_intp_y(x_fine,y_fine),'g')
        xlim([xmin xmax])
        ylim([ymin ymax])
        hold off
    end
end

function vOut=nan2zeroTriScatteredInterp(x,y,dtIn,vIn,method)
    F=TriScatteredInterp(dtIn,vIn,method);
    vOut=F(x,y);
    checkVec=isnan(vOut);
    vOut(checkVec)=0;
end

end