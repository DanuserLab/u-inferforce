function correctMovieDisplacementField(movieData,varargin)
% correctMovieDisplacementField calculate the displacement field
%
% correctMovieDisplacementField 
%
% SYNOPSIS correctMovieDisplacementField(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
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
% Sangyoon Han, from Oct 2014
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(DisplacementFieldCorrectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
displFieldCorrProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(displFieldCorrProc,paramsIn);

%% Backup the original vectors to backup folder
if exist(p.OutputDirectory,'dir')
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
mkClrDir(p.OutputDirectory);
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',displFieldCorrProc.getName());
end

% Reading various constants
nFrames = movieData.nFrames_;

% Check displacement field process
iDisplFieldCalcProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,1);     
if isempty(iDisplFieldCalcProc)
    error(['Displacement field calculation has not been run! '...
        'Please run displacement field calculation prior to force field calculation!'])   
end

displFieldCalcProc=movieData.processes_{iDisplFieldCalcProc};
if ~displFieldCalcProc.checkChannelOutput
    error(['The channel must have a displacement field ! ' ...
        'Please calculate displacement field to all needed channels before '...
        'running force field calculation!'])
end
displParams = displFieldCalcProc.funParams_;
inFilePaths{1} = displFieldCalcProc.outFilePaths_{1};
displFieldCorrProc.setInFilePaths(inFilePaths);

% Set up the output directories
outputFile{1,1} = [p.OutputDirectory filesep 'displField.mat'];
outputFile{2,1} = [p.OutputDirectory filesep 'dispMaps.mat'];
mkClrDir(p.OutputDirectory);
displFieldCorrProc.setOutFilePaths(outputFile);

% get firstMask
iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
pDistProc = displFieldCalcProc.funParams_;
if ~isempty(iSDCProc)
    SDCProc=movieData.processes_{iSDCProc};
    if ~SDCProc.checkChannelOutput(pDistProc.ChannelIndex)
        error(['The channel must have been corrected ! ' ...
            'Please apply stage drift correction to all needed channels before '...
            'running displacement field calclation tracking!'])
    end
    refFrame = double(imread(SDCProc.outFilePaths_{2,pDistProc.ChannelIndex}));
else
    refFrame = double(imread(pDistProc.referenceFramePath));
end
firstMask=refFrame>0;
%% --------------- Displacement field correction ---------------%%% 

disp('Starting correcting displacement field...')
% Anonymous functions for reading input/output
displField=displFieldCalcProc.loadChannelOutput;


disp('Detecting and filtering vector field outliers...')
logMsg = 'Please wait, detecting and filtering vector field outliers';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;
useGrid=displParams.useGrid;

% %Parse input, store in parameter structure
% pd = parseProcessParams(displFieldCalcProc,paramsIn);

% Perform vector field outlier detection
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
if feature('ShowFigureWindows'),parfor_progress(nFrames); end

outlierThreshold = p.outlierThreshold;
parfor j= 1:nFrames
% for j= 1:nFrames
    % Outlier detection
    dispMat = [displField(j).pos displField(j).vec];
    % Take out duplicate points (Sangyoon)
    [dispMat,~,~] = unique(dispMat,'rows'); %dispMat2 = dispMat(idata,:),dispMat = dispMat2(iudata,:)
    displField(j).pos=dispMat(:,1:2);
    displField(j).vec=dispMat(:,3:4);

    if ~isempty(outlierThreshold)
        if useGrid
            if j==1
                disp('In previous step, PIV was used, which does not require the current filtering step. skipping...')
            end
        else
            [outlierIndex,sparselyLocatedIdx,~,neighborhood_distance(j)] = detectVectorFieldOutliersTFM(dispMat,outlierThreshold,1);
            %displField(j).pos(outlierIndex,:)=[];
            %displField(j).vec(outlierIndex,:)=[];
            dispMat(outlierIndex,3:4)=NaN;
            dispMat(sparselyLocatedIdx,3:4)=NaN;
        end
        % I deleted this part for later gap-closing
        % Filter out NaN from the initial data (but keep the index for the
        % outliers)
%         ind= ~isnan(dispMat(:,3));
%         dispMat=dispMat(ind,:);
 
        displField(j).pos=dispMat(:,1:2);
        displField(j).vec=dispMat(:,3:4);

        % I deleted this part because artificially interpolated vector can
        % cause more error or false force. - Sangyoon June 2013
%         % Filling all NaNs with interpolated displacement vectors -
%         % We also calculate the interpolated displacements with a bigger correlation length.
%         % They are considered smoothed displacements at the data points. Sangyoon
%         dispMat = [dispMat(:,2:-1:1) dispMat(:,2:-1:1)+dispMat(:,4:-1:3)];
%         intDisp = vectorFieldSparseInterp(dispMat,...
%             displField(j).pos(:,2:-1:1),...
%             pd.minCorLength,pd.minCorLength,[],true);
%         displField(j).vec = intDisp(:,4:-1:3) - intDisp(:,2:-1:1);
    end
    
    % Update the waitbar
%     if mod(j,5)==1 && feature('ShowFigureWindows')
%         tj=toc;
%         waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-j)/j)]));
%     end
    if feature('ShowFigureWindows'), parfor_progress; end
end
if feature('ShowFigureWindows'), parfor_progress(0); end

if p.fillVectors
    % Now this is the real cool step, to run trackStackFlow with known
    % information of existing displacement in neighbors
    % Check optional process Flow Tracking
    pStep2 = displParams;
    if ~isempty(iSDCProc)
        s = load(SDCProc.outFilePaths_{3,pStep2.ChannelIndex},'T');
        residualT = s.T-round(s.T);
        refFrame = double(imread(SDCProc.outFilePaths_{2,pStep2.ChannelIndex}));
    else
        refFrame = double(imread(pStep2.referenceFramePath));
        residualT = zeros(nFrames,2);
    end
    logMsg = 'Please wait, retracking untracked points ...';
    timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
    tic
    nFillingTries=5;
    for j= 1:nFrames
        % Read image and perform correlation
        if ~isempty(iSDCProc)
            currImage = double(SDCProc.loadChannelOutput(pStep2.ChannelIndex(1),j));
        else
            currImage = double(movieData.channels_(pStep2.ChannelIndex(1)).loadImage(j));
        end
        nTracked=1000;
        nFailed=0;
        for k=1:nFillingTries
            % only un-tracked vectors
            unTrackedBeads=isnan(displField(j).vec(:,1));
            ratioUntracked = sum(unTrackedBeads)/length(unTrackedBeads);
            if ratioUntracked<0.0001 || (nTracked==0 && nFailed>30)
                break
            end
            currentBeads = displField(j).pos(unTrackedBeads,:);
            neighborBeads = displField(j).pos(~unTrackedBeads,:);
            neighborVecs = displField(j).vec(~unTrackedBeads,:);
            % Get neighboring vectors from these vectors (meanNeiVecs)
            [idx] = KDTreeBallQuery(neighborBeads, currentBeads, (1+5*k/nFillingTries)*neighborhood_distance(j)); % Increasing search radius with further iteration
%             [idx] = KDTreeBallQuery(neighborBeads, currentBeads, (2-1.5*k/nFillingTries)*neighborhood_distance(j)); % Increasing search radius with further iteration
            % In case of empty idx, search with larger radius.
            emptyCases = cellfun(@isempty,idx);
            mulFactor=1;
            while any(emptyCases)
                mulFactor=mulFactor+0.5;
                idxEmpty = KDTreeBallQuery(neighborBeads, currentBeads(emptyCases,:), mulFactor*(1+5*k/nFillingTries)*neighborhood_distance(j));
                idx(emptyCases)=idxEmpty;
                emptyCases = cellfun(@isempty,idx);
            end
            % Subsample idx to reduce computing time
            % Calculate the subsampling rate
            leap = cellfun(@(x) max(1,round(length(x)/100)),idx,'Unif',false);
            idx = cellfun(@(x,y) x(1:y:end,1),idx,leap,'Unif',false);
            closeNeiVecs = cellfun(@(x) neighborVecs(x,:),idx,'Unif',false);
        %     meanNeiVecs = cellfun(@mean,closeNeiVecs,'Unif',false);

            [v,nTracked] = trackStackFlowWithHardCandidate(cat(3,refFrame,currImage),currentBeads,...
                pStep2.minCorLength,pStep2.minCorLength,'maxSpd',pStep2.maxFlowSpeed,...
                'mode',pStep2.mode,'hardCandidates',closeNeiVecs);%,'usePIVSuite', pStep2.usePIVSuite);
            if nTracked==0
                nFailed=nFailed+1;
            else
                nFailed=0;
            end

        %     displField(j).pos(unTrackedBeads,:)=currentBeads; % validV is removed to include NaN location - SH 030417
            displField(j).vec(unTrackedBeads,:)=[v(:,1)+residualT(j,2) v(:,2)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
        end
        disp(['Done for frame ' num2str(j) '/' num2str(nFrames) '.'])
        % Update the waitbar
        if feature('ShowFigureWindows')
            tj=toc;
            waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-j)/j)]));
        end
    end
    %Filtering again
    parfor j= 1:nFrames
        % Outlier detection
        dispMat = [displField(j).pos displField(j).vec];
        % Take out duplicate points (Sangyoon)
        [dispMat,~,~] = unique(dispMat,'rows'); %dispMat2 = dispMat(idata,:),dispMat = dispMat2(iudata,:)
        displField(j).pos=dispMat(:,1:2);
        displField(j).vec=dispMat(:,3:4);

        [outlierIndex,sparselyLocatedIdx] = detectVectorFieldOutliersTFM(dispMat,outlierThreshold*3,1);
        %displField(j).pos(outlierIndex,:)=[];
        %displField(j).vec(outlierIndex,:)=[];
        dispMat(outlierIndex,3:4)=NaN;
        dispMat(sparselyLocatedIdx,3:4)=NaN;

        displField(j).pos=dispMat(:,1:2);
        displField(j).vec=dispMat(:,3:4);
        if feature('ShowFigureWindows'), parfor_progress; end
    end
    if feature('ShowFigureWindows'), parfor_progress(0); end
end
% Here, if nFrame>1, we do inter- and extrapolation of displacement vectors
% to prevent sudden, wrong force field change.
if nFrames>1 && ~displParams.useGrid
    disp('Performing displacement vector gap closing ...')
    % Depending on stage drift correction, some beads can be missed in certain
    % frames. Now it's time to make the same positions for all frames
    % go through each frame and filter points to the common ones in
    % iMinPointFrame - this needs to be improved by checking intersection
    % of all frames to find truly common beads, once there is error here.
    mostCommonPos = displField(1).pos;
    for ii= 2:nFrames
        commonPos=intersect(displField(ii).pos,mostCommonPos,'rows');
        mostCommonPos = commonPos;
    end
    for ii= 1:nFrames
        [commonPos,ia,~]=intersect(displField(ii).pos,mostCommonPos,'rows');
        displField(ii).pos = commonPos;
        displField(ii).vec = displField(ii).vec(ia,:);
    end
    % going through each point, see if there is NaN at each displacment
    % history and fill the gap
    logMsg = 'Performing displacement vector gap closing ...';

    nPoints = length(displField(1).pos(:,1));
    for k=1:nPoints
        % build each disp vector history
        curVecX = arrayfun(@(x) x.vec(k,1),displField);
        curVecY = arrayfun(@(x) x.vec(k,2),displField);
        if any(isnan(curVecX)) && sum(~isnan(curVecX))/nFrames>0.6
            t = 1:length(curVecX);
            t_nn = t(~isnan(curVecX));
            curVecX2 = interp1(t_nn,curVecX(~isnan(curVecX)),t,'linear');
            curVecY2 = interp1(t_nn,curVecY(~isnan(curVecX)),t,'linear');
            for ii=find(isnan(curVecX))
                displField(ii).vec(k,:) = [curVecX2(ii) curVecY2(ii)];
            end
        else
            continue
        end
        if mod(k,5)==1 && feature('ShowFigureWindows')
            tj=toc;
            waitbar(k/nPoints,wtBar,sprintf([logMsg timeMsg(tj*(nPoints-k)/k)]));
        end
    end
else
    if displParams.useGrid
        disp('In previous step, PIV was used, which does not require the current gap closing step. skipping...')
    end
end

% Find rotational registration
if p.doRotReg, displField=perfRotReg(displField); end 

%% Displacement map creation - this is shifted version
[dMapIn, dmax, dmin, cropInfo,dMapXin,dMapYin,reg_grid] = generateHeatmapShifted(displField,displField,0);
% Insert traction map in forceField.pos 
disp('Generating displacement maps ...')
dMap = cell(1,nFrames);
dMapX = cell(1,nFrames);
dMapY = cell(1,nFrames);
displFieldShifted(nFrames)=struct('pos','','vec','');
for ii=1:nFrames
    % starts with original size of beads
    cur_dMap = zeros(size(firstMask));
    cur_dMapX = zeros(size(firstMask));
    cur_dMapY = zeros(size(firstMask));
    cur_dMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapIn{ii};
    cur_dMapX(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapXin{ii};
    cur_dMapY(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapYin{ii};
    dMap{ii} = cur_dMap;
    dMapX{ii} = cur_dMapX;
    dMapY{ii} = cur_dMapY;
    % Shifted displField vector field
    [grid_mat,iu_mat, ~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
   
    [displFieldShiftedpos,displFieldShiftedvec, ~, ~] = interp_vec2grid(grid_mat+iu_mat, iu_mat,[],grid_mat); %1:cluster size
    pos = [reshape(displFieldShiftedpos(:,:,1),[],1) reshape(displFieldShiftedpos(:,:,2),[],1)]; %dense
    disp_vec = [reshape(displFieldShiftedvec(:,:,1),[],1) reshape(displFieldShiftedvec(:,:,2),[],1)]; 

    displFieldShifted(ii).pos = pos;
    displFieldShifted(ii).vec = disp_vec;
end
disp('Saving ...')
save(outputFile{1},'displField','displFieldShifted','-v7.3');
save(outputFile{2},'dMap','dMapX','dMapY','-v7.3'); % need to be updated for faster loading. SH 20141106
displFieldCorrProc.setTractionMapLimits([dmin dmax])

%% Close waitbar
if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished correcting displacement field!')