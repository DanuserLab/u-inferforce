function [fx fy x_out y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,method)

if nargin < 5 || isempty(method) || strcmp(method,'new')
    % we collect the nodeID for the basis functions. The position of these
    % nodes is then given by p(nodeID,:). We go back to the initial p in
    % order to make sure that the interpolation is performed on exactly the
    % same triangulation as we started with. Otherwise there might arise
    % problems at the boundary, see below.
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
    nodeIDvec=vertcat(forceMesh.basis(:).nodeID);
    % The corresponding values are:
    nonZeroVec=horzcat(sol_coef(1:forceMesh.numBasis),sol_coef(forceMesh.numBasis+1:end));    
    
    % Create the force field from this information:
    pos=forceMesh.p;
    
    % some Basis function might have been set to zero due to insufficient
    % support:
    vec=zeros(size(pos));
    vec(nodeIDvec,:)=nonZeroVec;
    
    % interpolate only if necessary
    if ~isempty(x_out) && ~isempty(y_out)
        % display('This might yield wrong results for square lattice since the Matlab triangulation doesnt agree with Achim''s triangulations!!!')
        % then interpolate the solution:
        fx_intp=TriScatteredInterp(pos(:,1),pos(:,2),vec(:,1),'linear');
        fy_intp=TriScatteredInterp(pos(:,1),pos(:,2),vec(:,2),'linear');
        
        fx=fx_intp(x_out,y_out);
        fy=fy_intp(x_out,y_out);
        
        checkVec=isnan(fx);
        fx(checkVec)=0;
        
        checkVec=isnan(fy);
        fy(checkVec)=0;
    else
        x_out=pos(:,1);
        y_out=pos(:,2);
        
        fx=vec(:,1);
        fy=vec(:,2);
    end    

elseif strcmp(method,'old')
    %if no points of interests are specified, i.e x_out, y_out, then forces are 
    %calculated on the nodes of the force mesh:                                                      
    if nargin<5 || isempty(x_out)
        x_out=forceMesh.p(:,1);
        y_out=forceMesh.p(:,2);
    end
    fx=0*x_out;
    fy=0*y_out;


    %reconstructed forces are obtained by multiplying the coef with the
    %resepctive basis.
    %This can be done much quicker, since most of the functions evaluate to zero 
    %e.g. fx(x0,y0) is non-zero only for 1:forceMesh.numNodes and base functions
    %base(j).f_intp_x which have support that comprises (x0,y0)
    for j=1:forceMesh.numBasis
        class  = forceMesh.basis(j).class;
        xShift = forceMesh.basis(j).node(1);
        yShift = forceMesh.basis(j).node(2);

        % These are the contributions from baseFunc(1):
        fx = fx+sol_coef(j)*forceMesh.basisClass(class).basisFunc(1).f_intp_x(x_out-xShift,y_out-yShift);
        fy = fy+sol_coef(j)*forceMesh.basisClass(class).basisFunc(1).f_intp_y(x_out-xShift,y_out-yShift);

        % These are the contributions from baseFunc(2):
        fx = fx+sol_coef(j+forceMesh.numBasis)*forceMesh.basisClass(class).basisFunc(2).f_intp_x(x_out-xShift,y_out-yShift);
        fy = fy+sol_coef(j+forceMesh.numBasis)*forceMesh.basisClass(class).basisFunc(2).f_intp_y(x_out-xShift,y_out-yShift);
    end
else
    display('Please specify the method')
end

% The following implementation suffers problems at the field boundary when basis
% functions (e.g. with insufficient support) are disregarded. The
% interpolation of the solution in these boundary regions is not always correct
% since the triangulation might be different from the initial triangulation
% saved in forceMesh.dt. This results from simply appending the zero nodes at
% the end of the non-zero nodes. However, this problem occurs only
% at the boundary but not inside the domain. This "minor" problem has been
% solved by both implementations above:


% pos=vertcat(forceMesh.basis(:).node);
% vec=horzcat(sol_coef(1:forceMesh.numBasis),sol_coef(forceMesh.numBasis+1:end));
% 
% % append the zero nodes:
% if isfield(forceMesh,'zeroNodes')
%     pos=vertcat(pos,forceMesh.zeroNodes);
%     vec=vertcat(vec,zeros(size(forceMesh.zeroNodes)));
% end
% 
% if ~isempty(x_out) && ~isempty(y_out)
%     % then interpolate the solution:
%     fx_intp=TriScatteredInterp(pos(:,1),pos(:,2),vec(:,1),'linear');
%     fy_intp=TriScatteredInterp(pos(:,1),pos(:,2),vec(:,2),'linear');
%     
%     fx=fx_intp(x_out,y_out);
%     fy=fy_intp(x_out,y_out);
%     
%     checkVec=isnan(fx);
%     fx(checkVec)=0;
%     
%     checkVec=isnan(fy);
%     fy(checkVec)=0;
% else
%     x_out=pos(:,1);
%     y_out=pos(:,2);
%     
%     fx=vec(:,1);
%     fy=vec(:,2);
% end