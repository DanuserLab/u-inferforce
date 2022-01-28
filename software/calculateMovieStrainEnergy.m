function []=calculateMovieStrainEnergy(movieData, varargin)
%function [SE_FOV,SE_FAs,SE_CellSeg]=quantifyMovieStrainEnergy(curMD)
%quantifies from traction map the strain energy for entire field of view
%(FOV), segmented FAs, and cell segmentation when the cell segmentation
%information is there.
% Unit is in femto Joule (1e-15 J)
% Sangyoon Han March, 2016
% calculateMovieStrainEnergy calculate the strain energy and total force
% out of the force field and the displacement field
%
% calculateMovieStrainEnergy 
%
% SYNOPSIS calculateMovieStrainEnergy(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
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

% Sangyoon J. Han, May 2017

disp(['Working on ' movieData.getFullPath '...'])
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('StrainEnergyCalculationProcess',1,0);
%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(StrainEnergyCalculationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
strainEnergyCalcProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(strainEnergyCalcProc,paramsIn);

%% Output setup
%% Backup the original vectors to backup folder
if exist(p.OutputDirectory,'dir')
    contents=dir(p.OutputDirectory);
    if any(arrayfun(@(x) ~x.isdir,contents)) % There should be something inside. Otherwise, backing up
        display('Backing up the original data')
        ii = 1;
        backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
        while exist(backupFolder,'dir')
            backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
            ii=ii+1;
        end
        mkdir(backupFolder);
        copyfile(p.OutputDirectory, backupFolder,'f')
    end
end
mkClrDir(p.OutputDirectory);
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',strainEnergyCalcProc.getName());
end
nFrames = movieData.nFrames_;

%% Load the forcefield
TFMPackage = movieData.getPackage(movieData.getPackageIndex('TFMPackage'));
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};

% tractionImgFolder=p.OutputDirectory;
% imgPath = [tractionImgFolder filesep 'imgs'];
% dataPath = [tractionImgFolder filesep 'data'];
% tifForcePath = [imgPath filesep 'tifForce'];
% tifBSensorPath = [imgPath filesep 'tifBSensor'];
% epsPath = [imgPath filesep 'eps'];
% figPath = [imgPath filesep 'figs'];
% if ~exist(tifForcePath,'dir')
%     system(['mkdir -p ' imgPath]);
%     system(['mkdir -p ' dataPath]);
%     system(['mkdir -p ' epsPath]);
%     system(['mkdir -p ' figPath]);
%     system(['mkdir -p ' tifForcePath]);
%     system(['mkdir -p ' tifBSensorPath]);
%     mkdir(imgPath);
%     mkdir(dataPath);
%     mkdir(epsPath);
%     mkdir(figPath);
%     mkdir(tifForcePath);
%     mkdir(tifBSensorPath);
% end
% Set up the output directories
outputFile{1,1} = [p.OutputDirectory filesep 'strainEnergyInFOV.mat'];
outputFile{1,2} = [p.OutputDirectory filesep 'strainEnergyInCell.mat'];
outputFile{1,3} = [p.OutputDirectory filesep 'forceBlobs.mat'];
if p.exportCSV
    outputFile{1,4} = [p.OutputDirectory filesep 'strainEnergyInFOV.csv'];
    outputFile{1,5} = [p.OutputDirectory filesep 'totalForceInFOV.csv'];
    outputFile{1,6} = [p.OutputDirectory filesep 'strainEnergyInCell.csv'];
    outputFile{1,7} = [p.OutputDirectory filesep 'totalForceInCell.csv'];
    outputFile{1,8} = [p.OutputDirectory filesep 'strainEnergyInForceBlobs.csv'];
    outputFile{1,9} = [p.OutputDirectory filesep 'indivForceInForceBlobs.csv'];
end    
outputFile{1,10} = [p.OutputDirectory filesep 'MovieStrainEnergyData_all.mat'];
strainEnergyCalcProc.setOutFilePaths(outputFile);

logMsg='Loading traction map...';
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end

tractionMaps=load(forceFieldProc.outFilePaths_{2});
try
    tMap = tractionMaps.tMap; % this is currently in Pa per pixel (1pix x 1pix)
catch
    % Set up the output directories
    outputDir = fullfile(forceFieldProc.funParams_.OutputDirectory,'tractionMaps');
    fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
    numStr = @(frame) num2str(frame,fString);
    outFileTMap=@(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
    tMap=cell(1,nFrames);
end
yModulus = forceFieldProc.funParams_.YoungModulus;
%% Load the displfield
iCorrectedDisplFieldProc = 3;
CorrectedDisplFieldProc=TFMPackage.processes_{iCorrectedDisplFieldProc};
logMsg='Loading displacement map...';
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
if ~isempty(CorrectedDisplFieldProc)
    try
        displMaps=load(CorrectedDisplFieldProc.outFilePaths_{2});
        dMap=displMaps.dMap; % this is currently in pix
    catch
        disp('Assuming displacement proportional to tMap because there was no iCorrectedDisplFieldProc found')
        dMap=cellfun(@(x) x/(yModulus),tMap,'UniformOutput',false);
    end
else
    disp('Assuming displacement proportional to tMap because there was no iCorrectedDisplFieldProc found')
    dMap=cellfun(@(x) x/(yModulus),tMap,'UniformOutput',false);
end
%% Calculate strain energy for FOV
% gridSpacing=1;
pixSize_mu=movieData.pixelSize_*1e-3; % in um/pixel
% factorConvert=gridSpacing^2*pixSize_mu^3/10^6;
% factorConvert=gridSpacing^2*(pixSize_mu*1e-6)^3;
areaConvert=pixSize_mu^2; % in um2/pixel
SE_FOV=struct('SE',zeros(nFrames,1),'area',zeros(nFrames,1),'SEDensity',zeros(nFrames,1));
totalForceFOV = zeros(nFrames,1);
SE_Cell=struct('SE',zeros(nFrames,1),'area',zeros(nFrames,1),'SEDensity',zeros(nFrames,1));
totalForceCell = zeros(nFrames,1);
SE_Blobs=struct('SE',zeros(nFrames,1),'nFA',zeros(nFrames,1),'areaFA',zeros(nFrames,1),...
    'SEDensity',zeros(nFrames,1),'avgFAarea',zeros(nFrames,1),'avgSEperFA',zeros(nFrames,1));
totalForceBlobs = struct('force',zeros(nFrames,1),'avgTraction',[],...
    'maxTraction',[],'forceBlobPixelIdxList',[]);
%% Get SDC
% Cell Boundary Mask 
iChan=2;
iBeadChan=1;
iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1); 
existSDC=false;
if ~isempty(iSDCProc)
    SDCProc=movieData.processes_{iSDCProc};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
    existSDC=true;
end
%% Get the cell boundary from the best mask
iMaskProcess=movieData.getProcessIndex('MaskIntersectionProcess');
existMask=false;
if ~isempty(iMaskProcess)
    maskProc = movieData.getProcess(iMaskProcess);
    existMask = true;
else
    iMaskProcess=movieData.getProcessIndex('MaskRefinementProcess');
    if ~isempty(iMaskProcess)
        maskProc = movieData.getProcess(iMaskProcess);
        existMask = true;
    end
end
%% Calculate strain energy and total force
tic
minSize = 20; % in pixel
minTraction = 100; % in Pa
forceField=load(forceFieldProc.outFilePaths_{1});
forceField = forceField.forceField;
gridSapcing=forceField(1).pos(2,2)-forceField(1).pos(1,2);
borderWidth=2*gridSapcing;
nTopBlobs=150; % the number of top maxForceBlobs in single frames
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];

logMsg='Quantifying strain energy and total force';
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
for ii=1:nFrames
    % Make sure if each tmap has its contents
    curTMap=tMap{ii};
    curDMap=dMap{ii};
    if isempty(curTMap)
        try
            curTMap=load(outFileTMap(ii),'cur_tMap');
            curTMap = curTMap.cur_tMap;
        catch
            curDisplField = CorrectedDisplFieldProc.loadChannelOutput(ii);
            [curTMapInsert,~,~,cropInfo]=generateHeatmapShifted(forceField(ii),curDisplField,0);
            curTMap = zeros(size(curDMap));
            curTMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3))=curTMapInsert{1};
        end
    end
    if p.useFOV || 1 % I will calculate this anyway
        % segment
        maskFOV = curTMap>0;%~isnan(curTMap); 
        % To crop, you can calculate the top, bottom, left, and right columns of the mask's "true" area by using any(). 
        anyCol=any(maskFOV);
        anyRow=any(maskFOV,2);
        row1 = find(anyRow,1);
        row2 = find(anyRow,1,'last');
        col1 = find(anyCol,1);
        col2 = find(anyCol,1,'last');

        tMapFOV = curTMap(row1:row2, col1:col2);
        dMapFOV = curDMap(row1:row2, col1:col2);

        % Erode FOV mask to remove the edge effect 
        maskShrunkenBorder = true(size(tMapFOV));
        maskShrunkenBorder = bwmorph(maskShrunkenBorder,'erode',borderWidth);

        SE_FOV.SE(ii)=1/2*sum(sum(dMapFOV.*tMapFOV.*maskShrunkenBorder))*(pixSize_mu*1e-6)^3*1e15; % in femto-Joule=1e15*(N*m)
        SE_FOV.area(ii)=sum(maskShrunkenBorder(:))*areaConvert; % this is in um2
        SE_FOV.SEDensity(ii)=SE_FOV.SE(ii)/SE_FOV.area(ii)*1e3; % J/m2
        
        totalForceFOV(ii) = sum(sum(tMapFOV.*maskShrunkenBorder))*areaConvert*1e-3; % in nN
    end
    
    if existMask && p.useCellMask
        maskCell = maskProc.loadChannelOutput(iChan,ii);
        if existSDC
            if isa(SDCProc,'EfficientSubpixelRegistrationProcess')
                maxX = 0;
                maxY = 0;
            else        
                maxX = ceil(max(abs(T(:, 2))));
                maxY = ceil(max(abs(T(:, 1))));
            end
            Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
            % Apply subpixel-wise registration to original masks

            Ibw = padarray(maskCell, [maxY, maxX]);
            maskCell = imtransform(Ibw, Tr, 'XData',[1 size(Ibw, 2)],'YData', [1 size(Ibw, 1)]);
        end
        tMapCell = curTMap(maskCell);
        dMapCell = curDMap(maskCell);

        SE_Cell.SE(ii)=1/2*sum(sum(dMapCell.*tMapCell))*(pixSize_mu*1e-6)^3*1e15; % in femto-Joule=1e15*(N*m)
        SE_Cell.area(ii)=sum(maskCell(:))*areaConvert; % this is in um2
        SE_Cell.SEDensity(ii)=SE_Cell.SE(ii)/SE_Cell.area(ii)*1e3; % J/m2
        totalForceCell(ii) = sum(sum(tMapCell))*areaConvert*1e-3; % in nN
    end
        
    if p.performForceBlobAnalysis
        maskForceBlob = blobSegmentThresholdTFM(tMapFOV,minSize,0,maskShrunkenBorder);
        maskForceBlob = bwmorph(maskForceBlob,'dilate',1);
%         maskForceBlob = padarray(maskForceBlob,[borderWidth borderWidth]);
        maskHighTraction=tMapFOV>minTraction;
        maskForceBlob = maskForceBlob & maskHighTraction;

        SE_Blobs.SE(ii)=1/2*sum(dMapFOV(maskForceBlob).*tMapFOV(maskForceBlob))*(pixSize_mu*1e-6)^3; % this is in Newton*m.
        SE_Blobs.SE(ii)=SE_Blobs.SE(ii)*1e15; % this is now in femto-Joule
        
        stats=regionprops(maskForceBlob,tMapFOV,'Area','PixelIdxList','Centroid','MinIntensity','MaxIntensity','MeanIntensity','WeightedCentroid');
        SE_Blobs.nFA(ii)=numel(stats);
        SE_Blobs.areaFA(ii) = sum(maskForceBlob(:))*areaConvert; % in um2
        SE_Blobs.avgFAarea(ii) = SE_Blobs.areaFA(ii)/SE_Blobs.nFA(ii);
        SE_Blobs.avgSEperFA(ii) = SE_Blobs.SE(ii)/SE_Blobs.nFA(ii); % still in femto-J
        SE_Blobs.SEDensity(ii)=SE_Blobs.SE(ii)/SE_Blobs.areaFA(ii)*1e3; % J/m2
        individualForceBlobs = arrayfun(@(x) x.MeanIntensity,stats);
        individualForceBlobMax = arrayfun(@(x) x.MaxIntensity,stats);
    %     individualForceBlobMin = arrayfun(@(x) x.MinIntensity,stats);
%         individualForceBlobCenters = arrayfun(@(x) x.WeightedCentroid,stats,'UniformOutput',false);
%         individualForceBlobAreas = arrayfun(@(x) x.Area,stats);

        totalForceBlobs.force(ii,1)=sum(tMapFOV(maskForceBlob))*areaConvert*1e-3; % in nN %*(pixSize_mu*1e-6)^2*1e9; % in nN
        totalForceBlobs.avgTraction{ii,1}=(individualForceBlobs); % in Pa
        totalForceBlobs.maxTraction{ii,1}=(individualForceBlobMax); % in Pa
        totalForceBlobs.forceBlobPixelIdxList{ii,1}=arrayfun(@(x) x.PixelIdxList,stats,'UniformOutput',false);
        % find an adhesion that contains top three max traction
    %     [individualForceBlobMaxSorted, topIDs]=sort(individualForceBlobMax,'descend');
    end
    
    % See if there is overlap between interiorMask and maskAdhesion - for
    % this I have to use AND operator
    
%     for k=1:SE_Blobs.nFA(ii)
% %         SE_Blobs.individualSE(ii).SE(k)=1/2*sum(dMapFOV(stats(k).PixelIdxList).*tMapFOV(stats(k).PixelIdxList))*(pixSize_mu*1e-6)^3*1e15;
% %         SE_Blobs.individualSE(ii).area(k)=stats(k).Area*areaConvert;
% %         SE_Blobs.individualSE(ii).SED(k)=SE_Blobs.individualSE(ii).SE(k)/SE_Blobs.individualSE(ii).area(k)*1e3;% J/m2.
%     end
%     [maxT]=max(tMapFOV(:));
%     SE_Blobs.maxT(ii) = maxT;
%     [SE_Blobs.maxSEFA(ii),~]=max(SE_Blobs.individualSE(ii).SE);
%     [SE_Blobs.maxSED_FA(ii),~]=max(SE_Blobs.individualSE(ii).SED);
    
    % Update the waitbar
    if feature('ShowFigureWindows')
        tj=toc;
        waitbar(ii/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-ii)/ii)]));
    end
end


%% Save
logMsg='Saving...';
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
if p.useFOV || 1
    save(outputFile{1},'SE_FOV','totalForceFOV','-v7.3');
    if p.exportCSV
        tableSE_FOV=struct2table(SE_FOV);
        writetable(tableSE_FOV,outputFile{4})
        tableForceFOV=table(totalForceFOV,'VariableNames',{'totalForceFOV'});
        writetable(tableForceFOV,outputFile{5})
    end
end
if p.useCellMask
    save(outputFile{2},'SE_Cell','totalForceCell','-v7.3'); % need to be updated for faster loading. SH 20141106
    if p.exportCSV
        tableSE_Cell=struct2table(SE_Cell);
        writetable(tableSE_Cell,outputFile{6})
        tableForceCell=table(totalForceCell,'VariableNames',{'totalForceCell'});
        writetable(tableForceCell,outputFile{7})
    end
end
if p.performForceBlobAnalysis
    save(outputFile{3},'SE_Blobs','totalForceBlobs','-v7.3');
    if p.exportCSV
        tableSE_Blobs=struct2table(SE_Blobs);
        writetable(tableSE_Blobs,outputFile{8})
        tableForceBlobs=table(totalForceBlobs.force,'VariableNames',{'totalForceBlobs'});
        writetable(tableForceBlobs,outputFile{9})
    end
end
save(outputFile{1,10},'SE_Blobs','totalForceBlobs', 'SE_Cell','totalForceCell','SE_FOV','totalForceFOV','-v7.3')
%% Close waitbar
if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished calculating strain energy and total force!')
