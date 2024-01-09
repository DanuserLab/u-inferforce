function calculateMovieForceField(movieData,varargin)
% calculateMovieForceField calculate the displacement field
%
% calculateMovieForceField 
%
% SYNOPSIS calculateMovieForceField(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
% Last updated by Sangyoon Han, Oct 2014
% Updates by Andrew R. Jamieson, Feb 2017

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('ForceFieldCalculationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(ForceFieldCalculationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
forceFieldProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(forceFieldProc,paramsIn);
p.usePaxImg = false;
p.saveBEMparams = true;
% p.lastToFirst = false;
p.LcurveFactor = 10;
p.divideConquer = 1; % If this is 9, grid is divided by 9 sub-grids where force field will be calculated to reduce memory usage. It's under refined construction.
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows'),
    wtBar = waitbar(0,'Initializing...','Name',forceFieldProc.getName());
    wtBarArgs={'wtBar',wtBar};
else
    wtBar=-1;
    wtBarArgs={};
end

% Reading various constants
nFrames = movieData.nFrames_;

% Check optional process Displacement field correction first
iDisplFieldProc =movieData.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);     
if isempty(iDisplFieldProc)
    iDisplFieldProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);     
end

if isempty(iDisplFieldProc)
    error(['Displacement field calculation has not been run! '...
        'Please run displacement field calculation prior to force field calculation!'])
end

displFieldProc=movieData.processes_{iDisplFieldProc};

if ~displFieldProc.checkChannelOutput()
    error(['Missing displacement field ! Please apply displacement field '...
        'calculation/correction  before running force field calculation!'])
end

% define resolution depending on the grid information in displacementField
% step
iDisplFieldCalProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);     
displFieldCalProc=movieData.processes_{iDisplFieldCalProc};
pDisp = parseProcessParams(displFieldCalProc);
try
    pDisp.useGrid;
catch
    pDisp.useGrid = false;
end
if pDisp.useGrid || pDisp.addNonLocMaxBeads
    p.highRes = false;
else
    p.highRes = true;
end
% Set up the input file
inFilePaths{1} = displFieldProc.outFilePaths_{1};
forceFieldProc.setInFilePaths(inFilePaths);

% Set up the output file
outputFile{1,1} = [p.OutputDirectory filesep 'forceField.mat'];
outputFile{2,1} = [p.OutputDirectory filesep 'tractionMaps.mat'];
outputFile{3,1} = [p.OutputDirectory filesep 'BEMParams.mat'];
% if  ~strcmp(p.solMethodBEM,'QR')
if p.useLcurve
    outputFile{4,1} = [p.OutputDirectory filesep 'Lcurve.fig'];
    outputFile{5,1} = [p.OutputDirectory filesep 'LcurveData.mat'];
end

% Add a recovery mechanism if process has been stopped in the middle of the
% computation to re-use previous results
firstFrame =1; % Set the starting frame to 1 by default
lastFrame=nFrames;
if exist(outputFile{1},'file')
    % Check analyzed frames
    s=load(outputFile{1},'forceField');
    frameForceField=~arrayfun(@(x)isempty(x.pos),s.forceField);
    
    if ~all(frameForceField) && ~all(~frameForceField)
        % Look at the first non-analyzed frame
        if p.lastToFirst 
            firstFrame = find(frameForceField);
            numFramesAnalyzed=length(firstFrame);
            lastFrame = firstFrame(1)-1;
        else
            firstFrame = find(~frameForceField,1);
            numFramesAnalyzed=firstFrame-1;
        end
        % Ask the user if display mode is active
        if ishandle(wtBar)
            recoverRun = questdlg(...
                ['A force field output has been dectected with ' ...
                num2str(numFramesAnalyzed) ' analyzed frames. Do you' ...
                ' want to use these results and continue the analysis'],...
                'Recover previous run','Yes','No','Yes');
            if ~strcmpi(recoverRun,'Yes'), firstFrame=1; lastFrame=nFrames; end
        end
%         if ~usejava('desktop')
%             firstFrame=1;
%         end
    end
end
% asking if you want to reuse the fwdMap again (you have to make sure that
% you are solving the same problem with only different reg. param.) -SH
reuseFwdMap = 'No';
if strcmpi(p.method,'FastBEM') && exist(outputFile{3,1},'file') 
    if usejava('desktop')
        reuseFwdMap = questdlg(...
            ['BEM parameters were dectected. Do you' ...
            ' want to use these parameter and overwrite the results?'],...
            'Reuse Fwdmap','Yes','No','No');
    end
end

% Backup the original vectors to backup folder
if firstFrame==1 && (strcmpi(reuseFwdMap,'No') || strcmpi(p.method,'FTTC')) && exist(outputFile{1,1},'file')
    display('Backing up the original data')
    backupFolder = [p.OutputDirectory ' Backup']; % name]);
    if exist(p.OutputDirectory,'dir')
        ii = 1;
        while exist(backupFolder,'dir')
            backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
            ii=ii+1;
        end
        mkdir(backupFolder);
        copyfile(p.OutputDirectory, backupFolder,'f')
    end
    forceField(nFrames)=struct('pos','','vec','','par','');
    forceFieldShifted(nFrames)=struct('pos','','vec','');
    displErrField(nFrames)=struct('pos','','vec','');
    distBeadField(nFrames)=struct('pos','','vec','');
    mkClrDir(p.OutputDirectory);
    M = [];
elseif strcmpi(reuseFwdMap,'Yes') && strcmpi(p.method,'FastBEM') && exist(outputFile{3,1},'file')
    fwdMapFile = load(outputFile{3,1},'M');
    M = fwdMapFile.M;
elseif strcmpi(p.method,'FastBEM') && strcmpi(reuseFwdMap,'No')  && exist(outputFile{3,1},'file') 
    % Load old force field structure 
    forceField=s.forceField;
    M = [];
else
    mkClrDir(p.OutputDirectory);
    M = [];
end
forceFieldProc.setOutFilePaths(outputFile);

%% --------------- Force field calculation ---------------%%% 

disp('Starting calculating force  field...')
if ~isempty(movieData.roiMaskPath_)
    maskArray = imread(movieData.roiMaskPath_);
else
    maskArray = movieData.getROIMask;
end
if min(min(maskArray(:,:,1))) == 0
    iStep2Proc = movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);
    step2Proc = movieData.processes_{iStep2Proc};
    pDisp = parseProcessParams(step2Proc,paramsIn);

    % Use mask of first frame to filter displacementfield
    iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    displFieldOriginal=displFieldProc.loadChannelOutput;
    
    if ~isempty(iSDCProc)
        SDCProc=movieData.processes_{iSDCProc};
        if ~SDCProc.checkChannelOutput(pDisp.ChannelIndex)
            error(['The channel must have been corrected ! ' ...
                'Please apply stage drift correction to all needed channels before '...
                'running displacement field calclation tracking!'])
        end
        %Parse input, store in parameter structure
        refFrame = double(imread(SDCProc.outFilePaths_{2,pDisp.ChannelIndex}));

        % Use mask of first frame to filter bead detection
        firstMask = refFrame>0; %false(size(refFrame));
        tempMask = maskArray(:,:,1);
        if isa(SDCProc,'EfficientSubpixelRegistrationProcess')
            s = load(SDCProc.outFilePaths_{3,pDisp.ChannelIndex},'T');
            T = s.T;
            meanYShift = round(T(1,1));
            meanXShift = round(T(1,2));
            firstMask = circshift(tempMask,[meanYShift meanXShift]);
            % Now I blacked out erroneously circularaly shifted bead image
            % portion - SH 20171008
            if meanYShift>=0 %shifted downward
                firstMask(1:meanYShift,:)=0;
            else %shifted upward
                firstMask(end+meanYShift:end,:)=0;
            end
            if meanXShift>=0 %shifted right hand side
                firstMask(:,1:meanXShift)=0;
            else %shifted left
                firstMask(:,end+meanXShift:end)=0;
            end
            
        else
            % firstMask(1:size(tempMask,1),1:size(tempMask,2)) = tempMask;
            tempMask2 = false(size(refFrame));
            y_shift = find(any(firstMask,2),1);
            y_lastNonZero = find(any(firstMask,2),1,'last');
            x_shift = find(any(firstMask,1),1);
            
            x_lastNonZero = find(any(firstMask,1),1,'last');
            % It is possible that I created roiMask based on tMap which is
            % based on reg_grid. In that case, I'll have to re-size firstMask accordingly
            % check if maskArray is made based on channel
            if (y_lastNonZero-y_shift+1)==size(tempMask,1) ...
                    && (x_lastNonZero-x_shift+1)==size(tempMask,2)
                tempMask2(y_shift:y_shift+size(tempMask,1)-1,x_shift:x_shift+size(tempMask,2)-1) = tempMask;
                if (y_shift+size(tempMask,1))>size(firstMask,1) || x_shift+size(tempMask,2)>size(firstMask,2)
                    firstMask=padarray(firstMask,[y_shift-1 x_shift-1],'replicate','post');
                end
            elseif size(firstMask,1)==size(tempMask,1) ... 
                    && size(firstMask,2)==size(tempMask,2) % In this case, maskArray (or roiMask) is based on reg_grid
                disp('Found that maskArray (or roiMask) is based on reg_grid')
                tempMask2 = tempMask;
            else
                error('Something is wrong! Please check your roiMask!')
            end
            firstMask = tempMask2 & firstMask;
        end
               
%         firstMask = false(size(refFrame));
%         tempMask = maskArray(:,:,1);
%         firstMask(1:size(tempMask,1),1:size(tempMask,2)) = tempMask; % This was wrong
        displField = filterDisplacementField(displFieldOriginal,firstMask);
    else
        firstMask = maskArray(:,:,1);
        displField = filterDisplacementField(displFieldOriginal,firstMask);
    end        
else
    displField=displFieldProc.loadChannelOutput;
    iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    if ~isempty(iSDCProc)
        SDCProc=movieData.processes_{iSDCProc};
        if ~SDCProc.checkChannelOutput(pDisp.ChannelIndex)
            error(['The channel must have been corrected ! ' ...
                'Please apply stage drift correction to all needed channels before '...
                'running displacement field calclation tracking!'])
        end
        %Parse input, store in parameter structure
        refFrame = double(imread(SDCProc.outFilePaths_{2,pDisp.ChannelIndex}));
        firstMask = false(size(refFrame));
    else
        firstMask = maskArray(:,:,1);
    end
end
   
% % Prepare displacement field for BEM
% if strcmpi(p.method,'fastBEM')
%     displField(end).par=0; % for compatibility with Achim parameter saving
%     displField=prepDisplForBEM(displField,'linear');
% end
% I don't think this is necessary anymore. - SH

% For Benedikt's software to work, the displacement field has to be
% interpolated on a rectangular grid, with an even number of grid points
% along each edge. Furthermore, one has to make sure that the noisy data 
% has not to be extrapolated. This may happen along the edges. To prevent
% this, extract the corner of the displacement grid, calculate how often 
% (even number) the optimal gridspacing fits into each dimension, then 
% place the regular grid centered to the orignal bounds. Thereby make sure 
% that the edges have been eroded to a certain extend. This is performed by
% the following function.
if ~p.highRes
    [reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displField,0.9,0); % we have to lower grid spacing because there are some redundant or aggregated displ vectors when additional non-loc-max beads were used for tracking SH170311
else
    [reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displField,1.0,0); %no dense mesh in any case. It causes aliasing issue!
end
distToBead = zeros(size(reg_grid(:,:,1)));
distToBead = distToBead(:);

disp('Calculating force field...')
logMsg = 'Please wait, calculating force field';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

if p.lastToFirst 
    frameSequence = lastFrame:-1:1;
else
    frameSequence = firstFrame:nFrames;
end

if strcmpi(p.method,'FastBEM')
    % if FastBEM, we calculate forward map and mesh only in the first frame
    % and then use parfor for the rest of the frames to calculate forces -
    % SH
    i=frameSequence(1); % For the first frame
    [grid_mat,~, ~,~] = interp_vec2grid(displField(i).pos, displField(i).vec,[],reg_grid);
    
    % basis function table name adjustment
    expectedName = ['basisClass' num2str(p.YoungModulus/1000) 'kPa' num2str(gridSpacing) 'pix.mat'];
    
    % Check if path exists - this step might not be needed because
    % sometimes we need to build a new table with a new name.
    if isempty(p.basisClassTblPath)
        disp(['Note, no path given for basis tables, outputing to movieData.mat path: ' ...
            movieData.movieDataPath_]);
        p.basisClassTblPath = fullfile(movieData.movieDataPath_, expectedName);
    else
        if exist(p.basisClassTblPath,'file')==2 
            disp('BasisFunctionFolderPath is valid.');
            if numel(whos('basisClassTbl', '-file', p.basisClassTblPath)) ~= 1
                disp(['basisFunction.mat not valid!' p.basisClassTblPath '. Will build a new basisFunction to this name.']);
            end
        else
            disp('New basisFunctionFolderPath is entered. Will build a new table and save in this path.');
        end
%         assert(exist(p.basisClassTblPath,'file')==2, 'basisFunctionFolderPath not valid!');
%         assert(numel(whos('basisClassTbl', '-file', p.basisClassTblPath)) == 1, ['basisFunction.mat not valid!' p.basisClassTblPath]);
    end

    % Sanity check the paths.
    basisFunctionFolderPath = fileparts(p.basisClassTblPath);
    assert(exist(basisFunctionFolderPath, 'dir')==7, 'basisFunctionFolderPath not valid!');

    % Check if basis file name is correctly formatted.
    expectedPath = fullfile(basisFunctionFolderPath, expectedName);
    if ~strcmp(expectedPath, p.basisClassTblPath)
        p.basisClassTblPath = expectedPath;
        disp(['basisClassTblPath has different name for estimated mesh grid spacing (' num2str(gridSpacing) '). ']);
        disp('Now the path is automatically changed to :')
        disp([expectedPath '.']);
    end
        
    % If grid_mat=[], then an optimal hexagonal force mesh is created
    % given the bead locations defined in displField:
    if p.usePaxImg && length(movieData.channels_)>1
        for i=frameSequence
            paxImage=movieData.channels_(2).loadImage(i);
            [pos_f, force, forceMesh, M, ~, ~, ~, ~]=...
                reg_FastBEM_TFM(grid_mat, displField, i, ...
                p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                'thickness',p.thickness/movieData.pixelSize_,'paxImg',paxImage,'pixelSize',movieData.pixelSize_);

            outputFile{3+i,1} = [p.OutputDirectory filesep 'BEMParams ' num2str(i) ' frame.mat'];

            disp(['saving forward map and custom force mesh at ' outputFile{3+i,1} '...'])
            save(outputFile{3+i,1},'forceMesh','M','-v7.3');
            display(['Done: solution for frame: ',num2str(i)]);
            % Fill in the values to be stored:
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            save(outputFile{1},'forceField');
        end
    elseif p.usePaxImg && i>1
        display('Loading BEM parameters... ')
        load(outputFile{3},'forceMesh','M','sol_mats');
        for i=frameSequence
            if p.usePaxImg && length(movieData.channels_)>1
                paxImage=movieData.channels_(2).loadImage(i);
                [pos_f,force,~]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                    displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,p.regParam,[],[], 'paxImg', paxImage, 'useLcurve', p.useLcurve);
            else
                [pos_f,force,~]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                    displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,sol_mats.L,[],[]);
            end
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            % Save each iteration (for recovery of unfinished processes)
            save(outputFile{1},'forceField');
            display(['Done: solution for frame: ',num2str(i)]);
            %Update the waitbar
            if mod(i,5)==1 && ishandle(wtBar)
                ti=toc;
                waitbar(i/nFrames,wtBar,sprintf([logMsg timeMsg(ti*nFrames/i-ti)]));
            end
        end
    else
        if ishandle(wtBar)
            waitbar(0,wtBar,sprintf([logMsg ' for first frame']));
        end
        if p.useLcurve
            if p.divideConquer>1 % divide and conquer
                nOverlap = 10; % the number of grid points to be overlapped
                % sub-divide grid_mat
                nLength = sqrt(p.divideConquer);
                [nRows, nCols,~]=size(grid_mat);
                subGridLimits(nLength,nLength) = struct('rowLim','','colLim','');
                subGrid = cell(nLength);
                subDisplField = cell(nLength);
                nRowBlock = ceil(nRows/nLength);
                nColBlock = ceil(nCols/nLength);
                forceInFullGrid = zeros(size(grid_mat));
                subForceInFullGrid = cell(nLength);% zeros(size(grid_mat));
                fullGrid = zeros(size(grid_mat));
                subFullGrid = cell(nLength); % zeros(size(grid_mat));
                % setting up the limits
                for jj=1:nLength % rows
                    for kk=1:nLength % columns
                        nRowFirst = max(1,1+nRowBlock*(jj-1)-nOverlap);
                        nRowSecond = min(nRows,1+nRowBlock*(jj)+nOverlap);
                        subGridLimits(jj,kk).rowLim=[nRowFirst nRowSecond];
                        nColFirst = max(1,1+nColBlock*(kk-1)-nOverlap);
                        nColSecond = min(nCols,1+nColBlock*(kk)+nOverlap);
                        subGridLimits(jj,kk).colLim=[nColFirst nColSecond];
                    end
                end
                
                nOverlapComb = 1;
                pp=0; %linear index
                % Combining...
                for jj=1:nLength % rows
                    for kk=1:nLength % columns
                        pp=pp+1;
                        curRowRange = subGridLimits(jj,kk).rowLim(1):subGridLimits(jj,kk).rowLim(2);
                        curColRange = subGridLimits(jj,kk).colLim(1):subGridLimits(jj,kk).colLim(2);
                        subGrid{jj,kk} = grid_mat(curRowRange,curColRange,:);
                        % Constructing sub-displacement field
                        subForceMask=zeros(size(firstMask));
                        subDispYFirst = max(1,subGrid{jj,kk}(1,1,2)-2*gridSpacing);
                        subDispYLast = min(size(firstMask,1),subGrid{jj,kk}(end,end,2)+2*gridSpacing);
                        subDispXFirst = max(1,subGrid{jj,kk}(1,1,1)-2*gridSpacing);
                        subDispXLast = min(size(firstMask,2),subGrid{jj,kk}(end,end,1)+2*gridSpacing);
                        subForceMask(subDispYFirst:subDispYLast,subDispXFirst:subDispXLast)=1;
                        subDisplField{jj,kk} = filterDisplacementField(displField,subForceMask);
                        [pos_f_sub{jj,kk}, force_sub{jj,kk}, forceMesh_sub{jj,kk}, M_sub{jj,kk}, ~, ~, ~,  sol_mats{jj,kk}]=...
                            reg_FastBEM_TFM(subGrid{jj,kk}, subDisplField{jj,kk}, i, ...
                            p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                            'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                            'useLcurve',p.useLcurve>0, ...
                            'LcurveFactor',p.LcurveFactor,'thickness',p.thickness/movieData.pixelSize_,...
                            'LcurveDataPath',outputFile{5,1},'LcurveFigPath',outputFile{4,1},'fwdMap',M,...
                            'lcornerOptimal',p.lcornerOptimal);
                        % Now we are adding code to use some shared part
                        % for matching the force coefficient outcome and adjusting regularization parameter
                        % Take the overlapping area
                        clear curForceInGrid
                        clear curGrid
                        curForceInGrid(:,:,1) = reshape(force_sub{jj,kk}(:,1),size(subGrid{jj,kk}(:,:,1),2),size(subGrid{jj,kk}(:,:,1),1)); 
                        curForceInGrid(:,:,2) = reshape(force_sub{jj,kk}(:,2),size(subGrid{jj,kk}(:,:,2),2),size(subGrid{jj,kk}(:,:,2),1));
                        curGrid(:,:,1) = reshape(pos_f_sub{jj,kk}(:,1),size(subGrid{jj,kk}(:,:,1),2),size(subGrid{jj,kk}(:,:,1),1));
                        curGrid(:,:,2) = reshape(pos_f_sub{jj,kk}(:,2),size(subGrid{jj,kk}(:,:,2),2),size(subGrid{jj,kk}(:,:,1),1));
                        % Define the overlapping area in the clock-wise
                        % fashion.
                        nRowFirst = max(1,1+nRowBlock*(jj-1)-nOverlapComb);
                        nRowSecond = min(nRows,1+nRowBlock*(jj)+nOverlapComb);
                        nColFirst = max(1,1+nColBlock*(kk-1)-nOverlapComb);
                        nColSecond = min(nCols,1+nColBlock*(kk)+nOverlapComb);
                        curRowRangeComb = nRowFirst:nRowSecond;
                        curColRangeComb =nColFirst:nColSecond;
                        %Insert in the full grid force
                        tempSubForceInFullGrid = zeros(size(grid_mat));
                        tempSubForceInFullGrid(curColRange,curRowRange,:) = curForceInGrid;
                        subForceInFullGrid{jj,kk}=zeros(size(grid_mat));
                        subForceInFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:)=tempSubForceInFullGrid(curColRangeComb,curRowRangeComb,:);
                        tempSubFullGrid = zeros(size(grid_mat));
                        tempSubFullGrid(curColRange,curRowRange,:) = curGrid;
                        subFullGrid{jj,kk}=zeros(size(grid_mat));
                        subFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:)=tempSubFullGrid(curColRangeComb,curRowRangeComb,:);
                        
                        % Actual combining
                        forceInFullGrid(curColRangeComb,curRowRangeComb,:)=...
                            subForceInFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:);
                        fullGrid(curColRangeComb,curRowRangeComb,:)=...
                            subFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:);
                    end
                end
                % re-format forceInFullGrid into vector format
                forceField(i).pos=[reshape(fullGrid(:,:,1),[],1), reshape(fullGrid(:,:,2),[],1)];
                forceField(i).vec=[reshape(forceInFullGrid(:,:,1),[],1), reshape(forceInFullGrid(:,:,2),[],1)];
                save(outputFile{1},'forceField');
            else
                [pos_f, force, forceMesh, M, pos_u, u, sol_coef,  sol_mats]=...
                    reg_FastBEM_TFM(grid_mat, displField, i, ...
                    p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                    'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                    'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                    'useLcurve',p.useLcurve>0, 'LcurveFactor',p.LcurveFactor,'thickness',p.thickness/movieData.pixelSize_,...
                    'LcurveDataPath',outputFile{5,1},'LcurveFigPath',outputFile{4,1},'fwdMap',M,...
                    'lcornerOptimal',p.lcornerOptimal);
                params = parseProcessParams(forceFieldProc,paramsIn);
                params.regParam = sol_mats.L;
                p.regParam = sol_mats.L;
                forceFieldProc.setPara(params);
                forceField(i).pos=pos_f;
                forceField(i).vec=force;
                save(outputFile{1},'forceField');
            end
        else
            if p.divideConquer>1 % divide and conquer
                nOverlap = 10; % the number of grid points to be overlapped
                % sub-divide grid_mat
                nLength = sqrt(p.divideConquer);
                [nRows, nCols,~]=size(grid_mat);
                subGridLimits(nLength,nLength) = struct('rowLim','','colLim','');
                subGrid = cell(nLength);
                subDisplField = cell(nLength);
                nRowBlock = ceil(nRows/nLength);
                nColBlock = ceil(nCols/nLength);
                forceInFullGrid = zeros(size(grid_mat));
                subForceInFullGrid = cell(nLength);% zeros(size(grid_mat));
                fullGrid = zeros(size(grid_mat));
                subFullGrid = cell(nLength); % zeros(size(grid_mat));
                % setting up the limits
                for jj=1:nLength % rows
                    for kk=1:nLength % columns
                        nRowFirst = max(1,1+nRowBlock*(jj-1)-nOverlap);
                        nRowSecond = min(nRows,1+nRowBlock*(jj)+nOverlap);
                        subGridLimits(jj,kk).rowLim=[nRowFirst nRowSecond];
                        nColFirst = max(1,1+nColBlock*(kk-1)-nOverlap);
                        nColSecond = min(nCols,1+nColBlock*(kk)+nOverlap);
                        subGridLimits(jj,kk).colLim=[nColFirst nColSecond];
                    end
                end
                
                nOverlapComb = 1;
                pp=0; %linear index
                % Combining...
                for jj=1:nLength % rows
                    for kk=1:nLength % columns
                        pp=pp+1;
                        curRowRange = subGridLimits(jj,kk).rowLim(1):subGridLimits(jj,kk).rowLim(2);
                        curColRange = subGridLimits(jj,kk).colLim(1):subGridLimits(jj,kk).colLim(2);
                        subGrid{jj,kk} = grid_mat(curRowRange,curColRange,:);
                        % Constructing sub-displacement field
                        subForceMask=zeros(size(firstMask));
                        subDispYFirst = max(1,subGrid{jj,kk}(1,1,2)-2*gridSpacing);
                        subDispYLast = min(size(firstMask,1),subGrid{jj,kk}(end,end,2)+2*gridSpacing);
                        subDispXFirst = max(1,subGrid{jj,kk}(1,1,1)-2*gridSpacing);
                        subDispXLast = min(size(firstMask,2),subGrid{jj,kk}(end,end,1)+2*gridSpacing);
                        subForceMask(subDispYFirst:subDispYLast,subDispXFirst:subDispXLast)=1;
                        subDisplField{jj,kk} = filterDisplacementField(displField,subForceMask);
                        [pos_f_sub{jj,kk}, force_sub{jj,kk}, forceMesh_sub{jj,kk}, M_sub{jj,kk}, ~, ~, ~,  sol_mats{jj,kk}]=...
                            reg_FastBEM_TFM(subGrid{jj,kk}, subDisplField{jj,kk}, i, ...
                            p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                            'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                            'useLcurve',p.useLcurve>0, 'thickness',p.thickness/movieData.pixelSize_);
                        % Now we are adding code to use some shared part
                        % for matching the force coefficient outcome and adjusting regularization parameter
                        % Take the overlapping area
                        clear curForceInGrid
                        clear curGrid
                        curForceInGrid(:,:,1) = reshape(force_sub{jj,kk}(:,1),size(subGrid{jj,kk}(:,:,1),2),size(subGrid{jj,kk}(:,:,1),1)); 
                        curForceInGrid(:,:,2) = reshape(force_sub{jj,kk}(:,2),size(subGrid{jj,kk}(:,:,2),2),size(subGrid{jj,kk}(:,:,2),1));
                        curGrid(:,:,1) = reshape(pos_f_sub{jj,kk}(:,1),size(subGrid{jj,kk}(:,:,1),2),size(subGrid{jj,kk}(:,:,1),1));
                        curGrid(:,:,2) = reshape(pos_f_sub{jj,kk}(:,2),size(subGrid{jj,kk}(:,:,2),2),size(subGrid{jj,kk}(:,:,1),1));
                        % Define the overlapping area in the clock-wise
                        % fashion.
                        nRowFirst = max(1,1+nRowBlock*(jj-1)-nOverlapComb);
                        nRowSecond = min(nRows,1+nRowBlock*(jj)+nOverlapComb);
                        nColFirst = max(1,1+nColBlock*(kk-1)-nOverlapComb);
                        nColSecond = min(nCols,1+nColBlock*(kk)+nOverlapComb);
                        curRowRangeComb = nRowFirst:nRowSecond;
                        curColRangeComb =nColFirst:nColSecond;
                        %Insert in the full grid force
                        tempSubForceInFullGrid = zeros(size(grid_mat));
                        tempSubForceInFullGrid(curColRange,curRowRange,:) = curForceInGrid;
                        subForceInFullGrid{jj,kk}=zeros(size(grid_mat));
                        subForceInFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:)=tempSubForceInFullGrid(curColRangeComb,curRowRangeComb,:);
                        tempSubFullGrid = zeros(size(grid_mat));
                        tempSubFullGrid(curColRange,curRowRange,:) = curGrid;
                        subFullGrid{jj,kk}=zeros(size(grid_mat));
                        subFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:)=tempSubFullGrid(curColRangeComb,curRowRangeComb,:);
                        
                        % Find the overlapping area
                        if ~(pp==1)
                            % find indices of overlapping grid
                            numOverlap = zeros(pp-1,1);
                            idxOverlap = cell(pp-1,1);
                            for qq=1:(pp-1)
                                jjp = floor((qq-1)/nLength)+1;
                                kkp = qq-(jjp-1)*nLength;
                                idxOverlap{qq}=subFullGrid{jjp,kkp}(:,:,1)== subFullGrid{jj,kk}(:,:,1) & ...
                                    subFullGrid{jjp,kkp}(:,:,2)== subFullGrid{jj,kk}(:,:,2) & ...
                                    subFullGrid{jjp,kkp}(:,:,1) & subFullGrid{jjp,kkp}(:,:,2) & ...
                                    subFullGrid{jj,kk}(:,:,1) & subFullGrid{jj,kk}(:,:,2);
                                numOverlap(qq) = sum(idxOverlap{qq}(:));
                            end
                            % Find the subGrid that has the most overlaps
                            [~,iMaxSubGrid]=max(numOverlap);
                            % Compare force between the current subgrid vs.
                            % iMaxSubGrid
                            jjpSelected = floor((iMaxSubGrid-1)/nLength)+1;
                            kkpSelected = iMaxSubGrid-(jjpSelected-1)*nLength;
                            idxSelectedX = idxOverlap{iMaxSubGrid};
                            idxSelectedY(:,:,2)=idxOverlap{iMaxSubGrid};
                            
                            % force vectors in iMaxSubGrid'th subgrid
                            clear forceInPrevSub gridInPrevSub forceInCurSub gridInCurSub
                            forceInPrevSub(:,1) = subForceInFullGrid{jjpSelected,kkpSelected}(idxSelectedX);
                            forceInPrevSub(:,2) = subForceInFullGrid{jjpSelected,kkpSelected}(idxSelectedY);
                            gridInPrevSub(:,1) = subFullGrid{jjpSelected,kkpSelected}(idxSelectedX);
                            gridInPrevSub(:,2) = subFullGrid{jjpSelected,kkpSelected}(idxSelectedY);
                            
                            % force vectors in the current subgrid
                            forceInCurSub(:,1) = subForceInFullGrid{jj,kk}(idxSelectedX);
                            forceInCurSub(:,2) = subForceInFullGrid{jj,kk}(idxSelectedY);
                            gridInCurSub(:,1) = subFullGrid{jj,kk}(idxSelectedX);
                            gridInCurSub(:,2) = subFullGrid{jj,kk}(idxSelectedY);
                            
                            L=p.regParam;
                            % calculate difference
                            diffNorm=norm(forceInCurSub)-norm(forceInPrevSub);
                            % Inner product between the two vector field
                            innerProdAll = forceInPrevSub*forceInCurSub';
                            innerProd = diag(innerProdAll);
                            innerProdSum = sum(abs(innerProd));
                            prevSelfProdAll = forceInPrevSub*forceInPrevSub';
                            prevSelfProd = diag(prevSelfProdAll);
                            prevSelfProdSum = sum(abs(prevSelfProd));
                            ratioProdSum=innerProdSum/prevSelfProdSum; %<1 if the current subgrid has underestimating force than overlap from previous subgrid.
                            
                            regFactor=10;
%                             oldDiffNorm=0;
                            oldRatio = ratioProdSum;
                            % For debugging
                            sc=0.01;
                            figure,quiver(gridInPrevSub(:,1),gridInPrevSub(:,2),sc*forceInPrevSub(:,1),sc*forceInPrevSub(:,2),0,'k')
                            hold on,quiver(gridInCurSub(:,1),gridInCurSub(:,2),sc*forceInCurSub(:,1),sc*forceInCurSub(:,2),0,'r')
                            disp(['ratioProdSum: ' num2str(ratioProdSum) ])
                            accuFactor=0.1; %deviation by 10 % is allowed
                            while (ratioProdSum>(1+accuFactor) || ratioProdSum<(1-accuFactor)) && regFactor>=1.01 %while the difference in norm is more than 100 Pa,
                                % change regularization parameter appropriately
                                if ratioProdSum<1 % if the current norm is smaller, L should be smaller to yield less-underestimated solution
                                    if oldRatio>1
                                        regFactor=regFactor^0.9;
                                    end
                                    L=L*1/regFactor;
                                elseif ratioProdSum>1 % if the current dot product sum is larger, the larger L should be applied to yield suppressed solution.
                                    if oldRatio<1 
                                        regFactor=regFactor^0.9;
                                    end
                                    L=L*regFactor;
                                end
                                disp(['Testing L= ' num2str(L)])
                                oldRatio=ratioProdSum;
                                [pos_f_sub{jj,kk}, force_sub{jj,kk}, forceMesh_sub{jj,kk}, M_sub{jj,kk}, ~, ~, ~,  sol_mats{jj,kk}]=...
                                    reg_FastBEM_TFM(subGrid{jj,kk}, subDisplField{jj,kk}, i, ...
                                    p.YoungModulus, p.PoissonRatio, L, p.meshPtsFwdSol,p.solMethodBEM,...
                                    'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:}, 'fwdMap', M_sub{jj,kk},...
                                    'useLcurve',false, 'thickness',p.thickness/movieData.pixelSize_);
                                
                                clear curForceInGrid
                                clear curGrid
                                curForceInGrid(:,:,1) = reshape(force_sub{jj,kk}(:,1),size(subGrid{jj,kk}(:,:,1),2),size(subGrid{jj,kk}(:,:,1),1)); 
                                curForceInGrid(:,:,2) = reshape(force_sub{jj,kk}(:,2),size(subGrid{jj,kk}(:,:,2),2),size(subGrid{jj,kk}(:,:,2),1));
                                curGrid(:,:,1) = reshape(pos_f_sub{jj,kk}(:,1),size(subGrid{jj,kk}(:,:,1),2),size(subGrid{jj,kk}(:,:,1),1));
                                curGrid(:,:,2) = reshape(pos_f_sub{jj,kk}(:,2),size(subGrid{jj,kk}(:,:,2),2),size(subGrid{jj,kk}(:,:,1),1));
                                % Define the overlapping area in the clock-wise
                                % fashion.
                                %Insert in the full grid force
                                tempSubForceInFullGrid = zeros(size(grid_mat));
                                tempSubForceInFullGrid(curColRange,curRowRange,:) = curForceInGrid;
                                subForceInFullGrid{jj,kk}=zeros(size(grid_mat));
                                subForceInFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:)=tempSubForceInFullGrid(curColRangeComb,curRowRangeComb,:);
                                tempSubFullGrid = zeros(size(grid_mat));
                                tempSubFullGrid(curColRange,curRowRange,:) = curGrid;
                                subFullGrid{jj,kk}=zeros(size(grid_mat));
                                subFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:)=tempSubFullGrid(curColRangeComb,curRowRangeComb,:);
                                
                                forceInCurSub(:,1) = subForceInFullGrid{jj,kk}(idxSelectedX);
                                forceInCurSub(:,2) = subForceInFullGrid{jj,kk}(idxSelectedY);
                                gridInCurSub(:,1) = subFullGrid{jj,kk}(idxSelectedX);
                                gridInCurSub(:,2) = subFullGrid{jj,kk}(idxSelectedY);
%                                 diffNorm=norm(forceInCurSub)-norm(forceInPrevSub);
                                
                                innerProdAll = forceInPrevSub*forceInCurSub';
                                innerProd = diag(innerProdAll);
                                innerProdSum = sum(abs(innerProd));
                                ratioProdSum=innerProdSum/prevSelfProdSum; %<1 if the current subgrid has underestimating force than overlap from previous subgrid.
                                % For debugging
                                sc=0.01;
                                figure,quiver(gridInPrevSub(:,1),gridInPrevSub(:,2),sc*forceInPrevSub(:,1),sc*forceInPrevSub(:,2),0,'k')
                                hold on,quiver(gridInCurSub(:,1),gridInCurSub(:,2),sc*forceInCurSub(:,1),sc*forceInCurSub(:,2),0,'r')
                                disp(['ratioProdSum: ' num2str(ratioProdSum) ' at L=' num2str(L)])
                            end
                        end
                        
                        % Actual combining
                        forceInFullGrid(curColRangeComb,curRowRangeComb,:)=...
                            subForceInFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:);
                        fullGrid(curColRangeComb,curRowRangeComb,:)=...
                            subFullGrid{jj,kk}(curColRangeComb,curRowRangeComb,:);
                    end
                end
                % re-format forceInFullGrid into vector format
                forceField(i).pos=[reshape(fullGrid(:,:,1),[],1), reshape(fullGrid(:,:,2),[],1)];
                forceField(i).vec=[reshape(forceInFullGrid(:,:,1),[],1), reshape(forceInFullGrid(:,:,2),[],1)];
                save(outputFile{1},'forceField');
            else
                [pos_f, force, forceMesh, M, pos_u, u, sol_coef,  sol_mats]=...
                    reg_FastBEM_TFM(grid_mat, displField, i, ...
                    p.YoungModulus, p.PoissonRatio, p.regParam, p.meshPtsFwdSol,p.solMethodBEM,...
                    'basisClassTblPath',p.basisClassTblPath,wtBarArgs{:},...
                    'imgRows',movieData.imSize_(1),'imgCols',movieData.imSize_(2),...
                    'useLcurve',p.useLcurve>0, 'thickness',p.thickness/movieData.pixelSize_,'fwdMap',M);
                forceField(i).pos=pos_f;
                forceField(i).vec=force;
                save(outputFile{1},'forceField');
            end
        end
        % Error estimation
        % I will use forward matrix to estimate relative uncertainty of
        % calculated displacement field for each force node. - SH
        % 09/08/2015
        % First, get the maxima for each force node from M
        forceNodeMaxima = max(M);
%             [neigh,bounds,bdPtsID]=findNeighAndBds(p,t);
%         forceConfidence.pos = forceMesh.p;
        cellPosition = arrayfun(@(x) x.node, forceMesh.basis,'UniformOutput',false);
        forceConfidence.pos = cell2mat(cellPosition');
        forceConfidence.vec = reshape(forceNodeMaxima,[],2);
        % Make it relative
        maxCfd = max(forceNodeMaxima);
        forceConfidence.vec = forceConfidence.vec/maxCfd;
        u_org = vertcat(displField(i).vec(:,1),displField(i).vec(:,2));
%         if p.divideConquer>1
%             %reconstruct M from M_sub
%             
%             u_predict = M*sol_coef;
%         else
%             u_predict = M*sol_coef;
%         end
%         u_diff = u_org-u_predict;
%         u_diff_vec=reshape(u_diff,[],2);
%         displErrField(i).pos=pos_u;
%         displErrField(i).vec=u_diff_vec;
        % Distance to the closest bead from each force node
        % check if pos_u is already nan-clear in terms of u
        if length(pos_u(:))==length(u)
            beadsWhole = pos_u;
        else
            idNanU = isnan(u_org);
            pos_u = pos_u(~idNanU);
            beadsWhole = pos_u;
        end
        
%         parfor i=1:length(forceField(i).pos(:,1))
%             [~,distToBead(i)] = KDTreeClosestPoint(beadsWhole,forceField(i).pos(:,1));
%         end
        
        %             display('The total time for calculating the FastBEM solution: ')

        % The following values should/could be stored for the BEM-method.
        % In most cases, except the sol_coef this has to be stored only
        % once for all frames!
        if p.saveBEMparams && strcmpi(reuseFwdMap,'No') && p.divideConquer==1
            disp(['saving forward map and force mesh at ' outputFile{3} '...'])
            save(outputFile{3},'forceMesh','M','sol_mats','pos_u','u','-v7.3');
        elseif p.saveBEMparams && strcmpi(reuseFwdMap,'No') && p.divideConquer>1
            disp(['saving forward map and force mesh at ' outputFile{3} '...'])
            save(outputFile{3},'pos_f_sub','force_sub','M_sub','forceMesh_sub','sol_mats','-v7.3');
        end
        for i=frameSequence(2:end)
            % since the displ field has been prepared such
            % that the measurements in different frames are ordered in the
            % same way, we don't need the position information any
            % more. The displ. measurements are enough.
            display('5.) Re-evaluate the solution:... ')
            if p.usePaxImg && length(movieData.channels_)>1
                paxImage=movieData.channels_(2).loadImage(i);
                [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                    displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,p.regParam,[],[], 'paxImg', paxImage, 'useLcurve', p.useLcurve);
            else
                [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,displField(i).pos(:,1),displField(i).pos(:,2),...
                    displField(i).vec(:,1),displField(i).vec(:,2),forceMesh,sol_mats.L,[],[]);
            end
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            % Save each iteration (for recovery of unfinished processes)
            save(outputFile{1},'forceField');
            % Error estimation
            u_org = vertcat(displField(i).vec(:,1),displField(i).vec(:,2));
            u_predict = M*sol_coef;
            u_diff = u_org-u_predict;
            u_diff_vec=reshape(u_diff,[],2);
            displErrField(i).pos=pos_u;
            displErrField(i).vec=u_diff_vec;
            display(['Done: solution for frame: ',num2str(i)]);
            %Update the waitbar
            if mod(i,5)==1 && ishandle(wtBar)
                ti=toc;
                waitbar(i/nFrames,wtBar,sprintf([logMsg timeMsg(ti*nFrames/i-ti)]));
            end
        end
    end
else % FTTC
    reg_corner=p.regParam;
    for i=frameSequence
        [grid_mat,iu_mat, i_max,j_max] = interp_vec2grid(displField(i).pos, displField(i).vec,[],reg_grid);
        if p.useLcurve && i==frameSequence(1)
                [rho,eta,reg_corner,alphas] = calculateLcurveFTTC(grid_mat, iu_mat, p.YoungModulus,...
                    p.PoissonRatio, gridSpacing, i_max, j_max, p.regParam,p.LcurveFactor);
                [reg_corner,ireg_corner,~,hLcurve]=regParamSelecetionLcurve(alphas',eta,alphas,reg_corner,'manualSelection',true);
                save(outputFile{5,1},'rho','eta','reg_corner','ireg_corner');
                saveas(hLcurve,outputFile{4,1});
                close(hLcurve)
        end
        [pos_f,~,force,~,~,~] = reg_fourier_TFM(grid_mat, iu_mat, p.YoungModulus,...
            p.PoissonRatio, movieData.pixelSize_/1000, gridSpacing, i_max, j_max, reg_corner);
        forceField(i).pos=pos_f;
        forceField(i).vec=force;
    end
end
%% For calculation of traction map and prediction error map
% The drift-corrected frames should have independent channel
% ->StageDriftCorrectionProcess
disp('Creating traction map...')
tic
[tMapIn, tmax, tmin, cropInfo,tMapXin,tMapYin] = generateHeatmapShifted(forceField,displField,0);
display(['Estimated traction maximum = ' num2str(tmax) ' Pa.'])
toc
if strcmpi(p.method,'FastBEM')
    [fCfdMapIn, fCmax, fCmin, cropInfoFC] = generateHeatmapShifted(forceConfidence,displField,0);
    fCfdMapIn{1} = fCfdMapIn{1}/max(fCfdMapIn{1}(:));
end
% display(['Displacement error minimum = ' num2str(dEmax) ' pixel.'])

% for ii=frameSequence
%     distBeadField(ii).pos = forceField(ii).pos;
%     distBeadField(ii).vec = [distToBead zeros(size(distToBead))];
% end
% [distBeadMapIn, dBeadmax, dBeadmin] = generateHeatmapShifted(distBeadField,displField,0);

% display(['Distance to closest bead maximum = ' num2str(tmax) ' Pa.'])
%% Insert traction map in forceField.pos 
disp('Writing traction maps ...')
tMap = cell(1,nFrames);
tMapX = cell(1,nFrames);
tMapY = cell(1,nFrames);
fCfdMap = cell(1,1); %force confidence

% Set up the output directories
outputDir = fullfile(p.OutputDirectory,'tractionMaps');
mkClrDir(outputDir);
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);
outFileTMap=@(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];

% distBeadMap = cell(1,nFrames);
for ii=frameSequence
    % starts with original size of beads
    cur_tMap = zeros(size(firstMask));
    cur_tMapX = zeros(size(firstMask));
    cur_tMapY = zeros(size(firstMask));
%     cur_distBeadMap = zeros(size(firstMask));
    cur_tMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapIn{ii};
    cur_tMapX(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapXin{ii};
    cur_tMapY(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapYin{ii};
%     cur_distBeadMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = distBeadMapIn{ii};
    tMap{ii} = cur_tMap;
    tMapX{ii} = cur_tMapX;
    tMapY{ii} = cur_tMapY;
    if ii==1 && strcmpi(p.method,'FastBEM')
        cur_fCfdMap = zeros(size(firstMask));
        cur_fCfdMap(cropInfoFC(2):cropInfoFC(4),cropInfoFC(1):cropInfoFC(3)) = fCfdMapIn{ii};
        fCfdMap = cur_fCfdMap;
    end     
%     distBeadMap{ii} = cur_distBeadMap;
    % Shifted forceField vector field
    curDispVec = displField(ii).vec;
    curDispVec = curDispVec(~isnan(curDispVec(:,1)),:); % This will remove the warning 
    curDispPos = displField(ii).pos;
    curDispPos = curDispPos(~isnan(curDispVec(:,1)),:); % This will remove the warning 
    [grid_mat,iu_mat, ~,~] = interp_vec2grid(curDispPos, curDispVec,[],reg_grid);
%     [grid_mat,iu_mat, ~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
    displ_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 
    
    [forceFieldShiftedpos,forceFieldShiftedvec, ~, ~] = interp_vec2grid(forceField(ii).pos+displ_vec, forceField(ii).vec,[],grid_mat); %1:cluster size
    pos = [reshape(forceFieldShiftedpos(:,:,1),[],1) reshape(forceFieldShiftedpos(:,:,2),[],1)]; %dense
    force_vec = [reshape(forceFieldShiftedvec(:,:,1),[],1) reshape(forceFieldShiftedvec(:,:,2),[],1)]; 

    forceFieldShifted(ii).pos = pos;
    forceFieldShifted(ii).vec = force_vec;

    save(outFileTMap(ii),'cur_tMap','-v7.3');

end
% Fill in the values to be stored:
clear grid_mat;
clear iu;
clear iu_mat;
                
disp('Saving ...')
% save(outputFile{1},'forceField','forceFieldShifted','displErrField');
% save(outputFile{2},'tMap','tMapX','tMapY','dErrMap','distBeadMap'); % need to be updated for faster loading. SH 20141106
save(outputFile{1},'forceField','forceFieldShifted');
if strcmpi(p.method,'FastBEM')
    save(outputFile{2},'tMap','tMapX','tMapY','fCfdMap','-v7.3'); % need to be updated for faster loading. SH 20141106
else
    save(outputFile{2},'tMap','tMapX','tMapY','-v7.3'); % need to be updated for faster loading. SH 20141106
end
forceFieldProc.setTractionMapLimits([tmin tmax])
% forceFieldProc.setDisplErrMapLimits([dEmin dEmax])
% forceFieldProc.setDistBeadMapLimits([dBeadmin dBeadmax])

% Close waitbar
if feature('ShowFigureWindows') || ishandle(wtBar), close(wtBar); end

disp('Finished calculating force field!')
end
