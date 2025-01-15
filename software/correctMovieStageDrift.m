function correctMovieStageDrift(movieData,varargin)
% correctMovieStageDrift corrects stage drift using reference frame and channel
%
% correctMovieStageDrift 
%
% SYNOPSIS correctMovieStageDrift(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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

% Sebastien Besson, Sep 2011
% Andrew R. Jamieson Feb. 2017 - clean up code, and update MATLAB functions.

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('BeadTrackingCorrectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(StageDriftCorrectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
displFieldProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(displFieldProc,paramsIn);

%% Backup the original vectors to backup folder
if exist(p.OutputDirectory,'dir')
    display('Backing up the original data')
    ii = 1;
    backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
    while exist(backupFolder,'dir')
        backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
        ii=ii+1;
    end
    try
        mkdir(backupFolder);
    catch
        system(['mkdir -p ' backupFolder]);
    end
    copyfile(p.OutputDirectory, backupFolder,'f')
end
mkClrDir(p.OutputDirectory);
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',displFieldProc.getName());
else
    wtBar=-1;
end

% Reading various constants
imDirs  = movieData.getChannelPaths;
imageFileNames = movieData.getImageFileNames;
nFrames = movieData.nFrames_;

% Set up the input directories (input images)
inFilePaths = cell(3,numel(movieData.channels_));
for j = p.ChannelIndex
    inFilePaths{1,j} = imDirs{j};
end
displFieldProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(3,numel(movieData.channels_));
mkClrDir(p.OutputDirectory);
for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i)];
    mkClrDir(outFilePaths{1,i});
end
iTFM = movieData.getPackageIndex('TFMPackage');
TFMpackage = movieData.getPackage(iTFM);
SDCProcess = TFMpackage.getProcess(1);
[~,refName,refExt]=fileparts(SDCProcess.funParams_.referenceFramePath);
outFilePaths{2,p.ChannelIndex(1)} = [p.OutputDirectory filesep refName refExt];
outFilePaths{3,p.ChannelIndex(1)} = [p.OutputDirectory filesep 'transformationParameters.mat'];
displFieldProc.setOutFilePaths(outFilePaths);

%% --------------- Stage drift correction ---------------%%% 

disp('Starting correcting stage drift...')
% Anonymous functions for reading input/output
outFile=@(chan,frame) [outFilePaths{1,chan} filesep imageFileNames{chan}{frame}];

% Loading reference channel images and cropping reference frame
refFrame = double(imread(p.referenceFramePath));
croppedRefFrame = imcrop(refFrame,p.cropROI);

% Detect beads in reference frame
p.doSubPixReg = true;
beadsChannel = movieData.channels_(p.ChannelIndex(1));

stack = zeros([movieData.imSize_ nFrames]);
disp('Loading stack...');
for j = 1:nFrames, stack(:,:,j) = double(beadsChannel.loadImage(j)); end


if p.doPreReg % Perform pixel-wise registration by auto-correlation
    % Initialize subpixel transformation array
    preT=zeros(nFrames,2);

    disp('Calculating pixel-wise pre-registration...')
    logMsg = 'Please wait, performing pixel-wise pre-registration';
    timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
    tic;
   
    % Caclulate autocorrelation of the reference frame
    selfCorr = normxcorr2(croppedRefFrame,refFrame);
    [maxAutoScore , imax] = max(abs(selfCorr(:)));
    [rowPosInRef, colPosInRef] = ind2sub(size(selfCorr),imax(1));

    if ishandle(wtBar), waitbar(0,wtBar,sprintf(logMsg)); end
    for j= 1:nFrames       
        % Find the maximum of the cross correlation between the template and
        % the current image:
        xCorr  = normxcorr2(croppedRefFrame,stack(:,:,j)); 
        [maxXScore , imax] = max(abs(xCorr(:)));
        if maxXScore/maxAutoScore<0.3
            % this means the position roi-ed is wrong. skipping ..
            preT(j,:)=[0 0];
        else
            [rowPosInCurr, colPosInCurr] = ind2sub(size(xCorr),imax(1));

            % The shift is thus:
            rowShift = rowPosInCurr-rowPosInRef;
            colShift = colPosInCurr-colPosInRef;

            % and the Transformation is:
            preT(j,:)=-[rowShift colShift];
        end        
        % Update the waitbar
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*nFrames/j-tj)]));
        end
    end
    
    % Apply pixel-wise registration to thre reference channel
    disp('Applying pixel-wise pre-registration...')
    logMsg = 'Please wait, applying pixel-wise pre-registration';
    timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
    tic;

      
    if ishandle(wtBar), waitbar(0,wtBar,sprintf(logMsg)); end
    for j = 1:nFrames
        
        Tr = affine2d([1 0 0; 0 1 0; fliplr(preT(j, :)) 1]);
        stack(:,:,j) = imwarp(stack(:,:,j), Tr);
        
        % Update the waitbar
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*nFrames/j-tj)]));
        end
    end  
else
    % Initialize transformation array
    preT=zeros(nFrames,2);
end

if p.doSubPixReg
    disp('Determining PSF sigma from reference frame...')
    % Adaptation of psfSigma from bead channel image data
    psfSigma = getGaussianPSFsigmaFromData(refFrame,'Display',false);
    if isnan(psfSigma) || psfSigma>movieData.channels_(1).psfSigma_*3  
        if strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'Widefield') || movieData.pixelSize_>130
            psfSigma = movieData.channels_(1).psfSigma_*2; %*2 scale up for widefield
        elseif strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'Confocal')
            psfSigma = movieData.channels_(1).psfSigma_*0.79; %*4/7 scale down for  Confocal finer detection SH012913
        elseif strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'TIRF')
            psfSigma = movieData.channels_(1).psfSigma_*3/7; %*3/7 scale down for TIRF finer detection SH012913
        else
            error('image type should be chosen among Widefield, confocal and TIRF!');
        end
        disp(['PSF sigma could not be determined by data due to abnormal distribution. Determined sigma using microscope setting: ' num2str(psfSigma)])
    else
        disp(['Determined sigma: ' num2str(psfSigma)])
    end

    assert(~isempty(psfSigma), ['Channel ' num2str(p.ChannelIndex(1)) ' have no '...
        'estimated PSF standard deviation. Pleae fill in the emission wavelength '...
        'as well as the pixel size and numerical aperture of the movie']);

    disp('Detecting beads in the reference frame...')
    pstruct = pointSourceDetection(croppedRefFrame, psfSigma, 'alpha', p.alpha);
    beads = [pstruct.x' pstruct.y'];

    % Select only beads  which are minCorLength away from the border of the
    % cropped reference frame
    beadsMask = true(size(croppedRefFrame));
    erosionDist=ceil((p.minCorLength+1+floor(p.maxFlowSpeed))/2);
    beadsMask(erosionDist:end-erosionDist,erosionDist:end-erosionDist)=false;
    indx=beadsMask(sub2ind(size(beadsMask),ceil(beads(:,2)),ceil(beads(:,1))));
    beads(indx,:)=[];
    assert(size(beads,1)>=10, ['Insufficient number of detected beads (less than 10): current number: ' num2str(length(beads)) '.']);
end

% Initialize transformation array
T=zeros(nFrames,2);
minNumVectors = 10;
if p.doSubPixReg
    disp('Calculating subpixel-wise registration...')
    logMsg = 'Please wait, performing sub-pixel registration';
    timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
    tic;
    
    % Perform sub-pixel registration
    if ishandle(wtBar), waitbar(0,wtBar,sprintf(logMsg)); end
    flow=cell(nFrames,1);
    for j= 1:nFrames

        % Stack reference frame and current frame and track beads displacement
        corrStack =cat(3,imcrop(refFrame,p.cropROI),imcrop(stack(:,:,j),p.cropROI));

        delta = trackStackFlow(corrStack,beads,...
            p.minCorLength,p.minCorLength,'maxSpd',p.maxFlowSpeed,'mode','accurate');
        
        %The transformation has the same form as the registration method from
        %Sylvain. Here we take simply the median of the determined flow
        %vectors. We take the median since it is less distorted by outliers.
        finiteFlow  = ~isinf(delta(:,1));
        nonNanFlow = ~isnan(delta(:,1));

        % If there is only few flow tracked and the distribution is
        % non-normal, we have to discard this tracking
        numTrackedVectors = sum(nonNanFlow);
        numAllVectors = length(nonNanFlow);
        if numTrackedVectors>minNumVectors && numTrackedVectors/numAllVectors>0.5 
            hNormal= kstest(delta(nonNanFlow,2));
            if hNormal
                T(j,:)=-nanmedian(delta(finiteFlow,2:-1:1),1);
            else            
                disp(['The ' num2str(numTrackedVectors) ' tracked vectors are non-normal. Assigning zeros in ' num2str(j) 'th frame...'])
            end
        else            
            disp(['There are only ' num2str(numTrackedVectors) ' tracked vectors among ' ...
                num2str(numAllVectors) ' total vectors. Assigning zeros in ' num2str(j) 'th frame...'])
        end

        % Remove infinite flow and save raw displacement under [pos1 pos2] format
        % into image coordinate system
        flow{j} = [beads(finiteFlow,2:-1:1) beads(finiteFlow,2:-1:1)+delta(finiteFlow,2:-1:1)];
       
        % Update the waitbar
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-j)/j)]));
        end
    end
end
T=T+preT;
disp('Applying stage drift correction...')
logMsg = @(chan) ['Please wait, correcting stage drift for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;

nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;
for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    refFrame = double(imread(p.referenceFramePath));
    
    for j= 1:nFrames
%          if i==1
%              %  Apply pixel-wise registration to the reference channel
%              Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(round(T(j, :))) 1]);
%          else
        try
             Tr = affine2d([1 0 0; 0 1 0; fliplr(T(j, :)) 1]);
        catch
            % Apply subpixel-wise registration to other channels
             Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(j, :)) 1]);
        end
        I = double(movieData.channels_(iChan).loadImage(j));
        try
            I2 = imwarp(I, Tr);%, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
        catch
            I2 = imtransform(I, Tr, 'XData',[1 size(I, 2)],'YData', [1 size(I, 1)]);
        end
          
        % Statistically test the local maxima to extract (significant) speckles
        imwrite(uint16(I2), outFile(iChan,j));
        
        % Update the waitbar
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
end

imwrite(uint16(refFrame), outFilePaths{2,p.ChannelIndex(1)});
if p.doSubPixReg
    save(outFilePaths{3,p.ChannelIndex(1)},'preT','T','flow');
else
    save(outFilePaths{3,p.ChannelIndex(1)},'preT','T');
end
% Close waitbar
if ishandle(wtBar), close(wtBar); end
disp('Finished correcting stage drift!')
