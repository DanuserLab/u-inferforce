function [v,numTrackedPoints,corLength,sigtValues] = trackStackFlowWithHardCandidate(stack,points,minCorL,varargin)
%trackStackFlowWithHardCandidate: Calculate the flow velocity from a stack
%of movie images with hard candidate. It is used for displacement
%correction process in TFM package to fill the untracked bead positions.
%This is a variant from trackStackFlow.
% SYNOPSIS :
%    v = trackStackFlow(stack,points,minCorL)
%    [v,corLen] = trackStackFlow(stack,points,minCorL,maxCorL,varargin)
%    [v,corLen,sigtVal] = trackStackFlow(stack,points,minCorL,maxCorL,varargin)
%
% INPUT :
%    stack : An image stack (i.e. of dimensions n x m x l) to be correlated
% 
%    points : A set of points (size nP x 2) expressed in the xy 
%             coordinate system where the velocity is calculated.
%
%    minCorL : The minimum side length of an image block (or band)
%            that is to be cross correlated over frames to detect flow
%            velocity. Optimal block size will be searched in the range
%            [minCorL maxCorL] for the detection of coherent flow pattern
%            behind noisy data.
%
%    maxCorL : Optional - The maximum side length of an image block (or band)
%            that is to be cross correlated over frames to detect flow
%            velocity. Optimal block size will be searched in the range
%            [minCorL maxCorL] for the detection of coherent flow pattern
%            behind noisy data.
%            If not input, will be set to minCorL
%
%    The following optional parameters can be set as parameter/value pairs:
%
%    'bgAvgImg': A stack of stationary background frames to be substracted 
%                during image correlation. Default is zeros matrix.
%
%    'maxSpd'  : A numerical value that specifies the maximum speed that can
%                be detected (in pixels/frame). The default is 10.
%
%    'bgMask':   A stack of background masks which is used to remove
%                background pixels from being used in correlation.
%
%    'minFeatureSize': The minimum feature size in the image. This is 
%                      measured as the diameter of features.
%                      Default, 11 pixels (typical speckle size).
%
%    'hardCandidates': A set of vectors in cell arrays (nPx1) that is mean
%                      vector of neighboring vectors
%
% OUTPUT :
%    v      : velocity vector of (size nP x2) expressed in the xy
%             coordinate system.
%
%    corLen : The optimal block length in the sense that it is the minimum
%             block length in the range [minCorLen, maxCorLe] that gives a
%             stable coherent flow.
%
%    sigtVal : The 1st, 2nd local maximum and the reference score for
%              significance test can also be output for use in
%              postprocessing.
%
% References:
% J. Li & G. Danuser, J. of microscopy, 220 150-167, 2005.
%
% Copyright (C) 2026, Danuser Lab - UTSouthwestern 
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

% Lin Ji, 2005
% Sebastien Besson, May 2011 (last modified Nov 2011)
% Adapted from imFlowTrack.m
% Sangyoon Han, October 2012 (last modified July 2013)

% Input check
ip= inputParser;
ip.addRequired('stack',@(x) isnumeric(x) && size(x,3)>=2);
ip.addRequired('points',@(x) isnumeric(x) && size(x,2)==2);
ip.addRequired('minCorL',@isscalar);
ip.addOptional('maxCorL',minCorL,@isscalar);
ip.addParamValue('maxSpd',40,@isscalar);
ip.addParamValue('bgMask',true(size(stack)),@(x) isequal(size(x),size(stack)));
ip.addParamValue('bgAvgImg', zeros(size(stack)),@isnumeric);
ip.addParamValue('minFeatureSize',11,@isscalar);
ip.addParamValue('mode','fast',@(x) ismember(x,{'fast','accurate','CCWS','CDWS'})); %This is about interpolation method
ip.addParamValue('scoreCalculation','xcorr',@(x) ismember(x,{'xcorr','difference'}));
ip.addParamValue('hardCandidates',[],@(x) iscell(x) || isempty(x))
% ip.addParamValue('usePIVSuite',false,@islogical);
ip.parse(stack,points,minCorL,varargin{:});
maxCorL=ip.Results.maxCorL;
maxSpd=ip.Results.maxSpd;
minFeatureSize=ip.Results.minFeatureSize;
bgMask=ip.Results.bgMask;
bgAvgImg=ip.Results.bgAvgImg;
mode=ip.Results.mode;
scoreCalculation=ip.Results.scoreCalculation;
closeNeiVecs = ip.Results.hardCandidates;

meanNeiVecs = cellfun(@(x) mean(x,1),closeNeiVecs,'Unif',false);
stdNeiVecs = cellfun(@(x) std(x,1),closeNeiVecs,'Unif',false);
anglesBetweenVecs = cell(numel(closeNeiVecs),1);
disp('Calculating angles between neighboring vectors...')
if feature('ShowFigureWindows'), parfor_progress(numel(closeNeiVecs)); end
parfor k=1:numel(closeNeiVecs) %for angles between vectors
    curNei = closeNeiVecs{k};
    numCurNei = size(curNei,1);
    for p=1:numCurNei-1
        vecP=curNei(p,:);
        for q=p+1:numCurNei
            vecQ = curNei(q,:);
            curAngle = acos(vecP*vecQ'/(norm(vecP)*norm(vecQ)));
            anglesBetweenVecs{k}=[anglesBetweenVecs{k} curAngle];
        end
    end
    if feature('ShowFigureWindows'), parfor_progress; end
end
if feature('ShowFigureWindows'), parfor_progress(0); end
stdAngleAll = cellfun(@std,anglesBetweenVecs);
% usePIVSuite = ip.Results.usePIVSuite;
% contWind = true;

% SH: Poly-fit version

% We automatically update the speed search radius until a high limit is
% reached. If no significant maximum is dclosenessThreshold*maxVNormetected, it means either the image
% quality is bad or the flow velocity is even higher than this limit.
maxSpdLimit = 2*maxSpd;

[imgW,imgL,numFrames] = size(stack);
x=points(:,1);
y=points(:,2);

%Initial maximum speed components in both direction.
initMaxFlowSpd = 20;
initMaxPerpSpd = 20;
closenessThreshold = 0.5; 
closenessThresholdpix = 2.5; %changed from 0.25 to account for not-interpolated maxV from 0.25;

%For isotropic correlation.
maxSpdLimit = max(maxSpdLimit,initMaxFlowSpd);

% Initialize output
nPoints = size(points,1);
v = zeros(nPoints,2);
corLength = minCorL*ones(nPoints,1);
sigtValues = NaN*ones(nPoints,3);

%We calculate a score for each sampling velocity. The score is an average of
% the normalized cross-correlation coefficient of an image block that moves
% over consecutive frames at the sampling velocity. We sample the velocity
% by sampling the components of the velocity in two orthogonal directions and
% in the range [-maxSpd maxSpd]. The two orthogonal directions are parallel to
% the two sides of the square block.
bandDir = [1 0];
perpDir = [-bandDir(2) bandDir(1)];

%We only use odd correlation lengths greater than 3 pixels
minCorL = max(3,minCorL+(1-mod(minCorL,2)));
maxCorL = max(minCorL,maxCorL+(1-mod(maxCorL,2)));

bAreaThreshold = min(0.95*minCorL^2,maxCorL^2*0.5);

% %Options for optimization.
% options = optimset('GradObj','on','Display','off');

% Creates the format string for the numerical indexes
L=length(num2str(nPoints));
strg=sprintf('%%.%dd',L);
backSpc =repmat('\b',1,L);

%Calculate the correlation coefficient for each sampling velocity at
% each point.
startTime = cputime;
fprintf(1,['   Start tracking (total: ' strg ' points): '],nPoints);

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    poolsize = feature('numCores');
else
    poolsize = poolobj.NumWorkers;
end
if isempty(gcp('nocreate'))
    try
        parpool('local',poolsize)
    catch
        disp('Please use matlab version 2015a or higher to use parallel computing') %matlabpool
    end
end % we don't need this any more.

if feature('ShowFigureWindows'), parfor_progress(nPoints); end
% inqryPoint=200;
% for k = inqryPoint
parfor k = 1:nPoints
% for k = 1:nPoints
%     fprintf(1,[strg ' ...'],k);
    
    sigtVal = [NaN NaN NaN];
    
    %Alway start with 'minCorL' for each new point.
    corL = minCorL;
    
    % candidate vector
    curCandVec = meanNeiVecs{k};
    curCandVecStd = stdNeiVecs{k};
    curStdAngle = stdAngleAll(k);
    
    pass = 0;
    while pass == 0 && corL <= maxCorL
        
        %Create kymograph around each point. 'bandDir' is used as the direction
        % of the kymographed line.
        xI = round(x(k));
        yI = round(y(k));
        if xI < 1 || xI > imgL || yI < 1 || yI > imgW
            %The point is outside the image. Mark it untrackable.
            pass = 0;
            corL = 2*maxCorL;
            continue;
        end
        %Always get back the initial max speed for new 'corL'.
        % We devide the max speed by 2 due the use of the while loop below.
        
        %If the quality of the score function is not good enough (pass == 0),
        % we increase the max sampling speed until the limit is reached.
        maxFlowSpd = max(10,round(abs(curCandVec(1)))*4);
        maxPerpSpd = max(10,round(abs(curCandVec(2)))*4);

        %Get sampling speed. Make sure it will not shift the template (block) outside of
        % the image area. We also use bigger stepsize for large speed.
        minSpdF = -min(floor(maxFlowSpd),max(0,xI-(corL+1)/2-1));
        maxSpdF = min(floor(maxFlowSpd),imgL-min(xI+(corL+1)/2+1,imgL));
        minSpdP = -min(floor(maxPerpSpd),max(0,yI-(corL-1)/2)-1);
        maxSpdP = min(floor(maxPerpSpd),imgW-min(yI+(corL-1)/2+1,imgW));
             
        %Debugging
        vF = minSpdF:maxSpdF;
        vP = minSpdP:maxSpdP;
        %End of debugging

        hCLL    = min(xI-1,(corL-1)/2)+max(-vF(1),0);
        hCLR    = min(imgL-xI,(corL-1)/2)+max(vF(end),0);
        hCWL    = min(yI-1,(corL-1)/2)+max(-vP(1),0);
        hCWR    = min(imgW-yI,(corL-1)/2)+max(vP(end),0);
        cropL   = hCLL+hCLR+1;
        cropW   = hCWL+hCWR+1;
        kym     = zeros(cropW,cropL,numFrames);
        for k2 = 1:numFrames
            kym(:,:,k2) = double(stack(yI-hCWL:yI+hCWR,xI-hCLL:xI+hCLR,k2));
        end
        kymMask   = squeeze(bgMask(yI-hCWL:yI+hCWR,xI-hCLL:xI+hCLR,:));
        kymAvgImg = bgAvgImg(yI-hCWL:yI+hCWR,xI-hCLL:xI+hCLR,:);

        %The index of zero velocity.
        zeroI = [find(vP==0) find(vF==0)];

        %The index of the center of image block in the kymographed or cropped
        % image.
        centerI = [hCWL+1,hCLL+1];

        if strcmp(mode,'CDWS')
            [score,blockIsTooSmall] = calScore(kym,centerI,corL,vP,vF,'CDWS',true);
        else
            [score,blockIsTooSmall] = calScore(kym,centerI,corL,vP,vF, ...
                'bAreaThreshold',bAreaThreshold,'kymMask',kymMask,'kymAvgImg',kymAvgImg,'mode',scoreCalculation);
        end

        if blockIsTooSmall
            %Tell the program to increase block size.
            pass = 0;
            maxFlowSpd = Inf;
            maxPerpSpd = Inf;
        else
            %Test the quality of the score function and find the index of the
            % maximum score.
            [~,locMaxI,sigtVal] = findMaxScoreI(score,zeroI,minFeatureSize,0.1);
            if ~isempty(curCandVec)
                % Here I use PIV result to find the best candidate
                % regardless of whether score passed or not from findMaxScoreI
                % Get the candidate vectors
                locMaxV = [vP(locMaxI(:,1)).' vF(locMaxI(:,2)).'];
                % What's the piv result in location closest to [xI, yI]
                options = {1:size(locMaxI,1),'Unif',false};

                curCandVecRot= curCandVec*[perpDir;bandDir]; %Rotated version we should use

                vecDiff = arrayfun(@(x) locMaxV(x,:)-curCandVecRot,options{:});
                orienDiff = arrayfun(@(x) acos(locMaxV(x,:)*curCandVecRot'/(norm(locMaxV(x,:))*norm(curCandVecRot))),options{1});
%                 distToMaxV2 = sqrt(sum((locMaxV-ones(size(locMaxV,1),1)*curCandVec).^2,2));
                magDiff = cellfun(@norm,vecDiff);
                
                [~,indDist]=sort(magDiff.*orienDiff.^2);
                ind = indDist(1);
                minAngle = orienDiff(ind);
                minDistDiff = magDiff(ind);
                if minDistDiff < 2*norm(curCandVecStd) && abs(minAngle) < curStdAngle % || (ind<=3 && sigtVal(ind)>0.8)
                    maxI = locMaxI(ind,:);
                    maxV = maxInterpfromScore(maxI,score,vP,vF,mode);
                    pass = 2;
                else
                    pass = 0;
                end
            else
                error('You have to set hardCandidates!');
            end
        end
        
        if pass == 0
            if corL == maxCorL
                corL = Inf;
            else
                corL = min(maxCorL,odd(corL*3/2));%floor(corL*3/2));
            end
        end
    end
    
    if pass == 0
        maxV = [NaN NaN];
        sigtVal = [NaN NaN NaN];
    elseif pass == 1
        maxV = maxInterpfromScore(maxI,score,vP,vF,mode);
    end
    if pass && strcmp(mode,'accurate')
        % subpixel continuous correlation score. This is more accurate than
        % interpolation from discrete scores - Sangyoon
        %         if norm(maxV)<1 && norm(maxV)>1e-5
        % new vF and vP (around maxV)
        halfCorL=(corL-1)/2;
        refineFactor = 10;%round(10*20/corL); % by this, the pixel value will be magnified.
        refineRange = 1.0; % in pixel
        if pass==1 || pass==2
            incFactor=1;
        end
%         minVP=max(round(maxV(1)*refineFactor) - refineRange*refineFactor,(yI-1-ceil(incFactor*corL)))*refineFactor)
%         maxVP=min(round(maxV(1)*refineFactor) + refineRange*refineFactor,(imgW-(yI+ceil(incFactor*corL)))*refineFactor);
        minVP=round(maxV(1)*refineFactor) - refineRange*refineFactor;
        vpMinTrimmed = yI*refineFactor+minVP-incFactor*halfCorL*refineFactor-refineFactor;
        if vpMinTrimmed<0
            minVP=incFactor*halfCorL*refineFactor+refineFactor-yI*refineFactor;
        end
        maxVP=round(maxV(1)*refineFactor) + refineRange*refineFactor;
        minVF=round(maxV(2)*refineFactor) - refineRange*refineFactor;
        vfMinTrimmed = xI*refineFactor+minVF-incFactor*halfCorL*refineFactor-refineFactor;
        if vfMinTrimmed<0
            minVF=incFactor*halfCorL*refineFactor+refineFactor-xI*refineFactor;
        end
        maxVF=round(maxV(2)*refineFactor) + refineRange*refineFactor;
        if minVP>=maxVP
            minVF=round(maxV(2)*refineFactor) - refineRange*refineFactor;
        end
        newvP = minVP:maxVP;
        newvF = minVF:maxVF;

        curXL = max(1,floor(xI-incFactor*halfCorL)); trimmedXL=max(incFactor*halfCorL-xI+1,0);
        curXR = min(imgL,ceil(xI+incFactor*halfCorL)); trimmedXR=max(xI+incFactor*halfCorL-imgL,0);
        curYL = max(1,floor(yI-incFactor*halfCorL)); trimmedYL=max(incFactor*halfCorL-yI+1,0);
        curYR = min(imgW,ceil(yI+incFactor*halfCorL)); trimmedYR=max(yI+incFactor*halfCorL-imgW,0);

        curXLfine = max(1,(xI-incFactor*halfCorL));
        curXRfine = min(imgL,(xI+incFactor*halfCorL));
        curYLfine = max(1,(yI-incFactor*halfCorL));
        curYRfine = min(imgW,(yI+incFactor*halfCorL));

        [curXI,curYI] = meshgrid(curXL:curXR,curYL:curYR);
        [fineXI,fineYI] = meshgrid(curXLfine:1/refineFactor:curXRfine,curYLfine:1/refineFactor:curYRfine);
        
        curXL2 = floor(xI+newvF(1)/refineFactor-incFactor*halfCorL+trimmedXL);
        curXR2 = ceil(xI+newvF(end)/refineFactor+incFactor*halfCorL-trimmedXR);
        curYL2 = floor(yI+newvP(1)/refineFactor-incFactor*halfCorL+trimmedYL);
        curYR2 = ceil(yI+newvP(end)/refineFactor+incFactor*halfCorL-trimmedYR);

        curXLfine2 = (xI+newvF(1)/refineFactor-incFactor*halfCorL+trimmedXL);
        curXRfine2 = (xI+newvF(end)/refineFactor+incFactor*halfCorL-trimmedXR);
        curYLfine2 = (yI+newvP(1)/refineFactor-incFactor*halfCorL+trimmedYL);
        curYRfine2 = (yI+newvP(end)/refineFactor+incFactor*halfCorL-trimmedYR);
                
        [curXI2,curYI2] = meshgrid(curXL2:curXR2,curYL2:curYR2);
        [fineXI2,fineYI2] = meshgrid(curXLfine2:1/refineFactor:curXRfine2,curYLfine2:1/refineFactor:curYRfine2);
        fineKym1 = interp2(curXI,curYI,stack(curYL:curYR,curXL:curXR,1),fineXI,fineYI);
        fineKym2 = interp2(curXI2,curYI2,stack(curYL2:curYR2,curXL2:curXR2,2),fineXI2,fineYI2);
        if strcmp(scoreCalculation,'xcorr')
            score_nxc2 = normxcorr2(fineKym1,fineKym2);
            K1=size(fineKym1,1); N1=size(fineKym2,1);
            K2=size(fineKym1,2); N2=size(fineKym2,2);
            score3 = score_nxc2(K1:N1,K2:N2); % normalized
        else
            newhCLL    = min(xI-1,halfCorL*refineFactor)-newvF(1);
            newhCWL    = min(yI-1,halfCorL*refineFactor)-newvP(1);
            %The index of the center of image block in the kymographed or cropped
            % image.
            newcenterI = [newhCWL+1,newhCLL+1];
            bI1 = round(newcenterI(1)-(size(fineKym1,1)-1)/2:newcenterI(1)+(size(fineKym1,1)-1)/2);
            bI2 = round(newcenterI(2)-(size(fineKym1,2)-1)/2:newcenterI(2)+(size(fineKym1,2)-1)/2);
            
            score3 = zeros(length(newvP),length(newvF));
            for j1 = 1:length(newvP)
                v1 = newvP(j1);
                for j2 = 1:length(newvF)
                    v2 = newvF(j2);
                    corrM = abs(fineKym1- fineKym2(bI1+v1,bI2+v2));

                    %Normalize the correlation coefficients.
                    score3(j1,j2) = sum(corrM(:));
                end
            end
            score3=ones(size(score3))-score3;
        end
        [~,maxI3ind] = max(score3(:));
        [ind3x,ind3y] = ind2sub(size(score3),maxI3ind);
        maxI3 = [ind3x,ind3y];
        maxVmagnified = maxInterpfromScore(maxI3,score3,newvP,newvF);
        
        maxV = maxVmagnified/refineFactor;
    end
    if pass && strcmp(mode,'CCWS')
        refineRange = 1; % in pixel
        maxIterCCWS = 4;
        oldmaxV = [1e6 1e6];
        k3 = 0;
        refineFactor = 10;% by this, the pixel value will be magnified.
        while norm(maxV-oldmaxV)>1e-6 && k3<maxIterCCWS
            k3 = k3+1;
            oldmaxV = maxV;
            maxV = contWindShift(maxV,kym,centerI,corL,refineFactor,refineRange);
            
            prev_refineFactor=refineFactor;
            refineFactor = prev_refineFactor *10;
            refineRange = refineRange/prev_refineFactor; % in pixel
        end
    end
    if ~isnan(maxV(1)) && ~isnan(maxV(2))
        rotv= maxV*[perpDir;bandDir];
        v(k,:) = [rotv(1) rotv(2)];
    else
        v(k,:) = [NaN NaN];
    end
    corLength(k) = corL;
    sigtValues(k,:) = sigtVal;
    
%     fprintf(1,[backSpc '\b\b\b\b']);
    if feature('ShowFigureWindows'), parfor_progress; end
end
if feature('ShowFigureWindows'), parfor_progress(0); end

nanInd = find(isnan(v(:,1)));
endTime = cputime;
fprintf(1,[strg '.\n'],nPoints);
fprintf(1,'   Tracking is done in %f sec (%f sec per point).\n', ...
    endTime-startTime,(endTime-startTime)/nPoints);
numTrackedPoints = nPoints-length(nanInd);
fprintf(1,'   Total tracked points: %d (out of %d).\n', ...
    numTrackedPoints,nPoints);

function [maxV] = contWindShift(maxV,kym,centerI,corL,refineFactor,refineRange)
newvP = round(maxV(1)*refineFactor)/refineFactor - refineRange:1/refineFactor:round(maxV(1)*refineFactor)/refineFactor + refineRange;
newvF = round(maxV(2)*refineFactor)/refineFactor - refineRange:1/refineFactor:round(maxV(2)*refineFactor)/refineFactor + refineRange;
[score4] = calScore(kym,centerI,corL,newvP,newvF,'Continuous',true);
[~,maxI4ind] = max(score4(:));
[ind4x,ind4y] = ind2sub(size(score4),maxI4ind);
maxV = [newvP(ind4x),newvF(ind4y)];

function [pass,locMaxI,sigtVal] = findMaxScoreI(score,zeroI,minFeatureSize,sigThreshold)
%
% INPUT:
%    score : The cross-correlation score function.
%    zeroI : The index of 'score' that corresponds to zero velocity.
% OUTPUT:
%    pass  : If an unambiguous global maximum is found, pass = 1 is returned.
%            Otherwise, pass = 0 is returned indicating that the quality of the
%            score function is not good.
%    locMaxI : Index of local maximum whose scores pass the significant test.
%    sigtVal : A 1x3 vector that contains the scores of the 1st, 2nd local
%              maximum and the reference score.
if nargin < 4
    sigThreshold = 0.5; %0.5;
end

[m,n] = size(score);
numSamples = m*n;

% avg = sum(score(:))/numSamples;
% dev = sum(abs(score(:)-avg))/numSamples;

%Some threshold used for quality test.
minFeatureRadius   = max(1,ceil((minFeatureSize-1)/2)); %Unit: pixel.
closenessThreshold = 0.2;
maxNbDist          = 0.5;

%We divide the score into 5x5 pixels blocks centered at 'zeroI' which is
% the index for zero velocity and find the local maximum in each block
% and then do a significant test on them in the sense that they have to be
% significantly bigger than the average. If there are more than one
% significant local maximums, return 0. To return 1, the unique significant
% local maximum also needs to be close the center of the sample region.

ind1 = 1:5:m;
if ind1(end) < m
    if ind1(end) == m-1
        ind1(end) = m;
    else
        if length(ind1)>1
            ind1(end) = floor((ind1(end-1)+m)/2);
            ind1(end+1) = m;
        else % length(ind1)==1
            ind1(end+1) = m;
        end
    end
end
ind2 = 1:5:n;
if ind2(end) < n
    if ind2(end) == n-1
        ind2(end) = n;
    else
        if length(ind2)>1
            ind2(end) = floor((ind2(end-1)+n)/2);
            ind2(end+1) = n;
        else
            ind2(end+1) = n;
        end
    end
end

locMaxI    = [];
locMaxS    = [];
locAvgMinS = [];
locCount   = 0;
for k = 1:length(ind1)-1
    for j = 1:length(ind2)-1
        locScore = score(ind1(k):ind1(k+1),ind2(j):ind2(j+1));
        [tmp,index] = max(locScore,[],1);
        [maxS,i2]   = max(tmp);
        i1 = index(i2);
        %maxI = [ind1(k,i1) ind2(j,i2)];
        maxI = [ind1(k)+i1-1 ind2(j)+i2-1];
        
        %Check if it is true local maximum in the sense that it is bigger than
        % its surrounding pixels or it is a boundary maximum.
        if maxI(1) == 1
            yOffset = maxI(1):maxI(1)+2;
        elseif maxI(1) == m
            yOffset = maxI(1)-2:maxI(1);
        else
            yOffset = maxI(1)-1:maxI(1)+1;
        end
        if maxI(2) == 1
            xOffset = maxI(2):maxI(2)+2;
        elseif maxI(2) == n
            xOffset = maxI(2)-2:maxI(2);
        else
            xOffset = maxI(2)-1:maxI(2)+1;
        end
        if maxS >= max(max(score(yOffset,xOffset)))
            %maxS = sum(sum(score(yOffset,xOffset)))/9; %This caused flow
            %underestimation. We should use a single maximum value at the
            %maximum velocity position rather than averaging with neiboring
            %points. This can prevent a value at the border of the score
            %from not being captured as a miximum, which will lead to
            %expansion of correlation length. BTW, what was the reason of
            %averaging maximum score with neighboring scores? To prevent
            %very narrow peak from being true maximum? I don't think
            %that'll happen. - Sangyoon 3/2/2013
            
            %The following 'avgMinS' is used in the calibration of the
            % 'baseS' below. It is the averge of scores around the local
            % maximum (excluding the maximum).
            avgMinS = (sum(sum(score(yOffset,xOffset)))-score(maxI(1),maxI(2)))/8;
            
            %Further check if it is close to any previouse
            % local maximum.
            %Distance to previouse local maximum.
            dist = zeros(locCount,1);
            for jj = 1:locCount
                dist(jj) = norm(maxI-locMaxI(jj,:));
            end
            [minD,minDI] = min(dist);
            if locCount == 0 || minD >= ...
                    max(2,norm(locMaxI(minDI,:)-zeroI)*closenessThreshold)
                locCount = locCount+1;
                
                locMaxI(locCount,:)  = maxI;
                locMaxS(locCount)    = maxS;
                locAvgMinS(locCount) = avgMinS;
            elseif maxS > locMaxS(minDI)
                locMaxI(minDI,:)  = maxI;
                locMaxS(minDI)    = maxS;
                locAvgMinS(minDI) = avgMinS;
            end
        end
    end
end

[maxS,maxInd] = max(locMaxS);
maxI = locMaxI(maxInd,:);
[locMaxS,desI] = sort(locMaxS,'descend');
locMaxI    = locMaxI(desI,:);
maxI       = locMaxI(1,:);
locAvgMinS = locAvgMinS(desI);
avgMinS    = locAvgMinS(1); %Note: this is not global 'minS'. It is the
% minimum score among the 9 elements around
% 'maxI'. It is used in the calibartion of
% 'baseS' below.
maxINorm = norm(maxI-zeroI);

if size(locMaxI,1) == 1
    sigtVal = [maxS 0 maxS];
    if maxI(1) < m/4 || maxI(1) > 3*m/4 || ...
            maxI(2) < n/4 || maxI(2) > 3*n/4
        pass = 0;
    else
        pass = 1;
    end
    return;
end

%Calculate the distance from all local maximum to the global maximum. It
%will be used to determine the neiborhood for calculating local average
%around the global maximum point.
dist = sqrt((locMaxI(2:end,1)-maxI(1)).^2+(locMaxI(2:end,2)-maxI(2)).^2);
minD = min(dist);

%Calculate the local averge around the maximun point. This local average
% will be used to calculate 'baseS', the reference score.
offset = max(minFeatureRadius,ceil(min(minD,maxINorm)*maxNbDist));
if maxI(1) > offset
    y0 = maxI(1)-offset;
else
    y0 = 1;
end
yNeighbor = y0:min(m,y0+2*offset);

if maxI(2) > offset
    x0 = maxI(2)-offset;
else
    x0 = 1;
end
xNeighbor = x0:min(n,x0+2*offset);
locAvg = sum(sum(score(yNeighbor,xNeighbor)))/ ...
    length(yNeighbor)/length(xNeighbor);

baseS  = locMaxS;
for k = 2:length(locMaxS)
    %Anisotropically adapted reference score:
    % Calculate local minimum along the line between the two local maximums
    % for the testing of local maximum significance.
    %The distance between the two local maximum.
    xD = locMaxI(k,2)-maxI(2);
    yD = locMaxI(k,1)-maxI(1);
    if abs(xD) > abs(yD)
        xShift = 0:sign(xD):xD;
        yShift = floor(xShift*yD/xD+0.5);
    else
        yShift = 0:sign(yD):yD;
        xShift = floor(yShift*xD/yD+0.5);
    end
    lineMinS = min(score(m*(maxI(2)+xShift-1)+maxI(1)+yShift));
    baseS(k)   = max(lineMinS,locAvg);
    
    %In very rare cases, 'baseS' can be even bigger than 'maxS' since our
    % 'maxS' is the average of the 9 elements around 'maxI'. So, we need
    % the following control to make sure 'baseS' is always less than
    % 'maxS'.
    if baseS(k) > maxS
        baseS(k) = min(avgMinS,baseS(k));
    end
end
if length(baseS) == 1
    sigtVal = [maxS 0 maxS];
elseif baseS(2) == maxS
    sigtVal = [maxS locMaxS(2) maxS/2];
else
    sigtVal = [maxS locMaxS(2) baseS(2)];
end
inSigI = find(maxS-baseS > 0 & locMaxS-baseS < (maxS-baseS)*sigThreshold);
locMaxI(inSigI,:) = [];
locMaxS(inSigI)   = [];

if isempty(locMaxS)
    pass = 0;
    return;
elseif length(locMaxS) > 1
    pass = 0;
    return;
end

if maxI(1) < min(m/20,2*minFeatureRadius) || maxI(1) > max(19*m/20,m-2*minFeatureRadius) || ...
        maxI(2) < min(n/20,2*minFeatureRadius) || maxI(2) > max(19*n/20, n-2*minFeatureRadius)
    pass = 0;
    return;
end

pass = 1;


function minV = minInterpfromScore(maxI2,score,vP,vF)
% Sangyoon: I made a change for this refining process to
% use parabola approximation. Once parabola fit is
% too much apart from integer maxV (maxVorg), I
% started to use the fmincon again for more correct refining process.
% parabola approximation
% input:    maxI2       :index for maxV in score
%           score       :score
%           vP,vF       :velocity range
% output:   maxV2       :refined velocity

subv = 1; % radius of subgroup for subscore
maxVorg  = [vP(maxI2(1)) vF(maxI2(2))];

if (maxI2(1)-subv)>=1 && (maxI2(1)+subv)<=size(score,1)...
   && (maxI2(2)-subv)>=1 && (maxI2(2)+subv)<=size(score,2)
    subv = 1; % radius of subgroup for subscore
    sub_score = score(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv),...
                        max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));
    % my field of interest
    subvP = vP(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv));
    subvF = vF(max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));

    [subvFG,subvPG]=meshgrid(subvF,subvP);
    subvF1D = reshape(subvFG,[],1);
    subvP1D = reshape(subvPG,[],1);
    sub_score1D = reshape(sub_score,[],1);

    % starting point estimation SH based on discretized maxV (-b/2a =
    % maxVorg(2)) in quadratical expression to avoid the random starting point warning SH
    asp = 4.9; %decided empirically
    bsp = -2*asp*maxVorg(2);
    csp = asp;
    dsp = -2*csp*maxVorg(1);
    esp = -6.7; %arbitrary number
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint', [asp,bsp,csp,dsp,esp]); 
    f = fittype('a*x^2+b*x+c*y^2+d*y+e','independent', {'x', 'y'}, 'dependent', 'z','option',s);
    sf = fit( [subvF1D, subvP1D], sub_score1D, f);

    px = [sf.a sf.b sf.e]; py = [sf.c sf.d sf.e];
    minV = [roots(polyder(py)) roots(polyder(px)) ];
end

