function [pos_f,force,forceMesh,M,pos_u,u,sol_coef,sol_mats]=...
    reg_FastBEM_TFM(grid_mat, displField, frame, yModu_Pa, pRatio, regParam, varargin)
% Synopsis [pos_f force forceMesh M pos_u u sol_coef sol_mats]=reg_FastBEM_TFM(grid_mat, displField, frame, yModu_Pa, pRatio, regParam, meshPtsFwdSol, solMethodBEM,varargin)
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

% Input check
ip =inputParser;
ip.addRequired('grid_mat');
ip.addRequired('displField',@isstruct);
ip.addRequired('frame',@isscalar);
ip.addRequired('yModu_Pa',@isscalar);
ip.addRequired('pRatio',@(x) isequal(pRatio,.5));
ip.addRequired('regParam',@isscalar);
ip.addOptional('meshPtsFwdSol',[],@(x)isscalar(x) ||isempty(x));
ip.addOptional('solMethodBEM','QR',@ischar);
ip.addParamValue('basisClassTblPath','',@ischar);
ip.addParamValue('LcurveDataPath','',@ischar);
ip.addParamValue('LcurveFigPath','',@ischar);
ip.addParamValue('wtBar',-1,@isscalar);
ip.addParamValue('imgRows',@isscalar);
ip.addParamValue('imgCols',@isscalar);
ip.addParamValue('thickness',472,@isscalar); % default assuming 34 um with 72 nm/pix resolution
ip.addParamValue('useLcurve',false,@islogical);
ip.addParamValue('lcornerOptimal','optimal',@ischar);
ip.addParamValue('LcurveFactor',@isscalar);
ip.addParamValue('paxImg',[],@ismatrix);
ip.addParamValue('forceMesh',[],@isstruct);
ip.addParamValue('pixelSize',@isscalar);
ip.addParamValue('strictBEM',false,@islogical);
ip.addParamValue('fwdMap',[],@ismatrix);
ip.parse(grid_mat, displField, frame, yModu_Pa, pRatio, regParam, varargin{:});
meshPtsFwdSol=ip.Results.meshPtsFwdSol;
solMethodBEM=ip.Results.solMethodBEM;
basisClassTblPath=ip.Results.basisClassTblPath;
LcurveDataPath=ip.Results.LcurveDataPath;
LcurveFigPath=ip.Results.LcurveFigPath;
lcornerOptimal=ip.Results.lcornerOptimal;
wtBar=ip.Results.wtBar;
imgRows = ip.Results.imgRows;
imgCols = ip.Results.imgCols;
thickness = ip.Results.thickness;    
paxImage = ip.Results.paxImg;
pixelSize = ip.Results.pixelSize;
forceMesh = ip.Results.forceMesh;
useLcurve = ip.Results.useLcurve;
LcurveFactor = ip.Results.LcurveFactor;
strictBEM = ip.Results.strictBEM;
M = ip.Results.fwdMap;

if isempty(grid_mat)
    % If no mesh is specified for the forces, we create a hexagonal mesh
    % that will be centered in the field of view. Here, only the first
    % frame will be used:
    [xvec,yvec]=createHexGridFromDisplField(displField,1);

    % plot the grid points:
    % figure(1)
    % plot(xvec,yvec,'o');
else
    xvec=grid_mat(:,:,1);
    yvec=grid_mat(:,:,2);
    
    xvec=xvec(:);
    yvec=yvec(:);
end

display('1.) Creating mesh & basis [~5sec]:...');
tic;
% keepBDPts=true; %this might lead to unmatching forward map that lead to
% diagonalized traction map
keepBDPts=false; % change made SH 2016.07.19
doPlot=0;
% strictBEM = false;
% if strcmp(solMethodBEM,'1NormReg') || strcmp(solMethodBEM,'1NormRegLaplacian')
%     strictBEM = true;
% end
if strictBEM
    xvec = displField(frame).pos(:,1);
    yvec = displField(frame).pos(:,2);
    idxNonan = ~isnan(displField(frame).vec(:,1));
    xvec = xvec(idxNonan);
    yvec = yvec(idxNonan);

    forceMesh=createMeshAndBasis(xvec,yvec,doPlot);
elseif isempty(paxImage)
    forceMesh=createMeshAndBasisFastBEM(xvec,yvec,keepBDPts,[],doPlot);
elseif isempty(forceMesh)
    forceMesh=createMeshAndBasisFromAdhesions(xvec,yvec,paxImage,displField(frame),pixelSize);
end
toc;
display('Done: mesh & basis!');

if strictBEM
    [fx,fy,x_out,y_out,M,pos_u,u,sol_coef,sol_mats] = ...
        BEM_force_reconstruction(displField(frame).pos(:,1),displField(frame).pos(:,2),...
        displField(frame).vec(:,1),displField(frame).vec(:,2),forceMesh,yModu_Pa,regParam,...
        [],[],'fast',meshPtsFwdSol,solMethodBEM,'wtBar',wtBar,'thickness',thickness,'useLcurve',useLcurve,...
        'LcurveFactor',LcurveFactor,'LcurveDataPath',LcurveDataPath, 'LcurveFigPath',LcurveFigPath,...
        'strictBEM',strictBEM);    
elseif isempty(paxImage)
    [fx,fy,x_out,y_out,M,pos_u,u,sol_coef,sol_mats] = ...
        BEM_force_reconstruction(displField(frame).pos(:,1),displField(frame).pos(:,2),...
        displField(frame).vec(:,1),displField(frame).vec(:,2),forceMesh,yModu_Pa,regParam,...
        [],[],'fast',meshPtsFwdSol,solMethodBEM,'basisClassTblPath',basisClassTblPath,...
        'wtBar',wtBar,'imgRows',imgRows,'imgCols',imgCols,'thickness',thickness,'useLcurve',useLcurve,...
        'LcurveFactor',LcurveFactor,'LcurveDataPath',LcurveDataPath, 'LcurveFigPath',LcurveFigPath,'fwdMap',M,...
        'lcornerOptimal',lcornerOptimal);
    % The units of fx and fy are the same as the input E, that is ususally Pa!
else
    xmin = min(forceMesh.p(:,1));
    xmax = max(forceMesh.p(:,1));
    ymin = min(forceMesh.p(:,2));
    ymax = max(forceMesh.p(:,2));
    [x_out,y_out] = meshgrid(xmin:xmax,ymin,ymax);
    x_out = reshape(x_out,[],1);
    y_out = reshape(y_out,[],1);
    
    [fx,fy,x_out,y_out,M,pos_u,u,sol_coef,sol_mats] = ...
        BEM_force_reconstruction(displField(frame).pos(:,1),displField(frame).pos(:,2),...
        displField(frame).vec(:,1),displField(frame).vec(:,2),forceMesh,yModu_Pa,regParam,...
        x_out,y_out,'slow',meshPtsFwdSol,solMethodBEM,'wtBar',wtBar,'imgRows',imgRows,...
        'imgCols',imgCols,'thickness',thickness,'paxImg',paxImage,'fwdMap',M);
    % The units of fx and fy are the same as the input E, that is ususally Pa!
end

pos_f=horzcat(x_out,y_out);
force=horzcat(   fx,   fy);

% 
% figure(100)
% quiver(x_out,y_out,fx,fy,'b')
% hold on
% quiver(displField(1).pos(:,1),displField(1).pos(:,2),displField(1).vec(:,1),displField(1).vec(:,2),'r')
% title('red: displacement field, blue: force field')
% hold off