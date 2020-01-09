function [] = setROIfromForcemap(MD)
% This function shows you (FTTC-based) TFM map and lets you draw ROI for
% your BEM (L1 or L2) force reconstruction. FTTC should've been processed
% before selecting this function:
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

% Load TFM map
tfmPackage=MD.getPackage(MD.getPackageIndex('TFMPackage'));
% Load force process
forceProc = tfmPackage.getProcess(4);
p = forceProc.funParams_;
% Load displ process
try
    displProc = tfmPackage.getProcess(3);
catch
    displProc = tfmPackage.getProcess(2);
end
displField=displProc.loadChannelOutput;
% SDC proc
try
    SDCProc = tfmPackage.getProcess(1);
catch
    SDCProc = tfmPackage.getProcess(2);
end
% Get force map
% try
%     tMap = load(forceProc.outFilePaths_{2});
%     tMap = tMap.tMap;
% catch
% If there is nothing, run the package with FTTC setting
funParams = forceProc.funParams_;
funParams.method='FTTC';
funParams.useLcurve=false;
forceProc.setPara(funParams);
%     forceProc.run
[reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displField,1,0);
endInd = numel(displField);
[grid_mat,iu_mat, i_max,j_max] = interp_vec2grid(displField(endInd).pos, displField(endInd).vec,[],reg_grid);
[pos_f,~,force,~,~,~] = reg_fourier_TFM(grid_mat, iu_mat, p.YoungModulus,...
    p.PoissonRatio, MD.pixelSize_/1000, gridSpacing, i_max, j_max, p.regParam);
curForceField.pos=pos_f;
curForceField.vec=force;
curDispField=displField(endInd);
[tMapIn, tmax, tmin, cropInfo] = generateHeatmapShifted(curForceField,curDispField,0);    
%     tMap = load(forceProc.outFilePaths_{2});
%     tMap = tMap.tMap;
% end
refFrame = double(imread(SDCProc.outFilePaths_{2,1}));

% Use mask of first frame to filter bead detection
firstMask = refFrame>0; %false(size(refFrame));
cur_tMap = zeros(size(firstMask));
cur_tMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapIn{1};

% Show the map
h2=figure; imshow(cur_tMap,[tmin tmax]), colormap jet;

% Draw the ROI
disp(['Draw rectangle for ROI for ' MD.movieDataPath_ '.'])
if isempty(MD.roiMaskPath_)
    h=imrect;
else
    roiMask=imread(MD.roiMaskPath_);
    boundROI=bwboundaries(roiMask);
    hold on, plot(boundROI{1}(:,2),boundROI{1}(:,1),'w')
    h=imrect;
end
ROI_rect = wait(h);
roiMask=createMask(h);

% Save it as ROI mask associated with MD
roiPath=[MD.outputDirectory_ filesep 'roiMask.tif'];
imwrite(roiMask,roiPath);
MD.setROIMaskPath(roiPath);
% maskArray = imread(MD.roiMaskPath_);
MD.roiMask=roiMask;
% maskArray = MD.getROIMask;
close(h2)
disp('ROI created!')
















