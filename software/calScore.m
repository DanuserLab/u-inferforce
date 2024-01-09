function [score,blockIsTooSmall,vP,vF] = calScore(kym,centerI,corL,vP,vF,varargin)
% centerI : The coordinate index of the image block center in the 'kym' image.
% SH: this local function was updated for bead tracking instead of speckle
% tracking. The first slice of kym contains the reference frame and the
% second one has the current bead image.
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

%For boundary points, the cutoff of background can make the effective image area too small for
% stable correlation. We check this and report it back as output.
blockIsTooSmall = 0;

if mod(corL,2) == 0, corL = corL+1; end

numFrames = size(kym,3);
kymLen    = size(kym,2);
kymWidth  = size(kym,1);

%Check additional parameters
ip =inputParser;
ip.addParamValue('bAreaThreshold',0.5*corL^2,@isscalar);
ip.addParamValue('kymMask',[],@islogical)
ip.addParamValue('Continuous',false,@islogical)
ip.addParamValue('CDWS',false,@islogical)
ip.addParamValue('kymAvgImg',zeros(size(kym)),@isnumeric)
ip.addParamValue('mode','xcorr',@(x) ismember(x,{'xcorr','difference'}));
ip.parse(varargin{:});
bAreaThreshold=ip.Results.bAreaThreshold;
kymMask=ip.Results.kymMask;
kymAvgImg=ip.Results.kymAvgImg;
bCont = ip.Results.Continuous;
bCDWS = ip.Results.CDWS;
mode = ip.Results.mode;

% score = zeros(length(vP),length(vF));
bI1 = round(centerI(1)-(corL-1)/2:centerI(1)+(corL-1)/2);
bI2 = round(centerI(2)-(corL-1)/2:centerI(2)+(corL-1)/2);

% elseif isempty(kymMask) && (numFrames==2)
%     %The index of the correlating image block in the big cropped image.
%     bI1e = centerI(1)-(corL-1)/2+vP(1):centerI(1)+(corL-1)/2+vP(end);
%     bI2e = centerI(2)-(corL-1)/2+vF(1):centerI(2)+(corL-1)/2+vF(end);
%     % fft-based cross-correlation. This can reduce computation time 
%     % especially for large velocity range - Sangyoon
% =======
bI1e = round(centerI(1)-(corL-1)/2+vP(1):centerI(1)+(corL-1)/2+vP(end));
bI2e = round(centerI(2)-(corL-1)/2+vF(1):centerI(2)+(corL-1)/2+vF(end));
% m = length(bI1e);%length(vP);
% n = length(bI2e);%length(vF);
% mi = length(vP);
% ni = length(vF);

%Find the part of the image block that is outside the cropped image and cut it off from the template.
% vP([find(bI1e<1) find(bI1e>kymWidth)-corL+1])=[];
% vF([find(bI2e<1) find(bI2e>kymLen)-corL+1])=[];

% this part needs to be a bit smarter. If bI1 or bI2 are cut, it should not
% affect vP or vF because only template size was changed. If bI1e or bI2e are cut, it should affect vP or vF.
% bI1(bI1<1 | bI1>kymWidth) = [];
% bI2(bI2<1 | bI2>kymLen) = [];
% bI1e(bI1e<1 | bI1e>kymWidth) = [];
% bI2e(bI2e<1 | bI2e>kymLen) = [];
idxBI1cutFirst = find(bI1<1);
idxBI1cutLast = find(bI1>kymWidth);
idxBI2cutFirst = find(bI2<1);
idxBI2cutLast = find(bI2>kymLen);
if ~isempty(idxBI1cutFirst)
    bI1(idxBI1cutFirst)=[];
    bI1e(idxBI1cutFirst)=[];
end
if ~isempty(idxBI2cutFirst)
    bI2(idxBI2cutFirst)=[];
    bI2e(idxBI2cutFirst)=[];
end
if ~isempty(idxBI1cutLast)
    bI1(idxBI1cutLast)=[];
    bI1e(end-length(idxBI1cutLast)+1:end)=[];
end
if ~isempty(idxBI2cutLast)
    bI2(idxBI2cutLast)=[];
    bI2e(end-length(idxBI2cutLast)+1:end)=[];
end

% For bI1e or bI2e, vP and vF are affected
idxBI1EcutFirst = find(bI1e<1);
idxBI1EcutLast = find(bI1e>kymWidth);
idxBI2EcutFirst = find(bI2e<1);
idxBI2EcutLast = find(bI2e>kymLen);
if ~isempty(idxBI1EcutLast)
    bI1e(idxBI1EcutLast)=[];
    vP(1:1+length(idxBI1EcutLast)-1)=[];
end
if ~isempty(idxBI2EcutLast)
    bI2e(idxBI2EcutLast)=[];
    vP(1:1+length(idxBI2EcutLast)-1)=[];
end
if ~isempty(idxBI1EcutFirst)
    bI1e(idxBI1EcutFirst)=[];
    vP(idxBI1EcutFirst)=[];
end
if ~isempty(idxBI2EcutFirst)
    bI2e(idxBI2EcutFirst)=[];
    vF(idxBI2EcutFirst)=[];
end

score = zeros(length(vP),length(vF));
if strcmp(mode,'difference')
    for j1 = 1:length(vP)
        v1 = vP(j1);
        for j2 = 1:length(vF)
            v2 = vF(j2);
            corrM = abs(kym(bI1,bI2,1:numFrames-1)- ...
                kym(bI1+v1,bI2+v2,2:numFrames));

            %Normalize the correlation coefficients.
            score(j1,j2) = sum(corrM(:));
        end
    end
    % Invert
    score=ones(size(score))-score;
elseif bCont && (numFrames==2)
    % normalized cross-correlation with continuous window shift (by
    % interpolation)
    % Here vP and vF are decimal numbers. 
    g1 = kym(bI1,bI2,1:numFrames-1);
    g1m = mean(g1(:)); % mean g1
    g1n = g1-g1m; % normalized g1
    for j1 = 1:length(vP)
        v1 = vP(j1);
        for j2 = 1:length(vF)
            v2 = vF(j2);
            %integer part of v1 and v2
            v1i = floor(v1);
            v2i = floor(v2);
            % decimal parts of v1 and v2
            x = v1-v1i;
            y = v2-v2i;
            
            g2 = (1-x)*(1-y)*kym(bI1+v1i,bI2+v2i,numFrames)...
                + x*(1-y)*kym(bI1+v1i+1,bI2+v2i,numFrames)...
                + y*(1-x)*kym(bI1+v1i,bI2+v2i+1,numFrames)...
                + x*y*kym(bI1+v1i+1,bI2+v2i+1,numFrames);
            g2m = mean(g2(:));
            g2n = g2 - g2m;
            corrM = g1n.* g2n;
            %Normalize the correlation coefficients.
            score(j1,j2) = sum(corrM(:));
        end
    end
    score = score/max(score(:));
elseif bCDWS && (numFrames==2)
    for j1 = 1:length(vP)
        v1 = vP(j1);
        for j2 = 1:length(vF)
            v2 = vF(j2);
            kym1 = kym(bI1,bI2,1);
            kym2 = kym(bI1+v1,bI2+v2,2);
            corrM = (kym1-mean(kym1(:))).*(kym2-mean(kym2(:)));

            %Normalize the correlation coefficients.
            bnormNorm1 = sqrt(sum(sum((kym1-mean(kym1(:))).^2)));
            bnormNorm2 = sqrt(sum(sum((kym2-mean(kym2(:))).^2)));
            score(j1,j2) = sum(corrM(:))/(bnormNorm1*bnormNorm2);
        end
    end
elseif isempty(kymMask)
    % normalized cross-correlation. This can make the correlation less sensitive to bright region in template window - Sangyoon
    K1 = length(bI1); K2 = length(bI2);
    N1 = length(bI1e); N2 = length(bI2e);
    % the correlation less sensitive to bright region in template window - Sangyoon
    score_nxc2 = normxcorr2(kym(bI1,bI2,1:numFrames-1),kym(bI1e,bI2e,2:numFrames));
    score = score_nxc2(K1:N1,K2:N2); % normalized
    
elseif min(min(kymMask(:,:,1))) == 1 || min(min(kymMask(:,:,1))) == 0
    %Background intensities are set to be zero. So, if there is no zero intensities, there is no
    % background.
    kymP2 = kym.*kym;
    
    %The norm of the kymographed image band at each frame.
    bNorm1 = sqrt(sum(sum(sum(kymP2(bI1,bI2,1:numFrames-1)))));
    
%     % fft-based cross-correlation. This can reduce computation time 
%     % especially for large velocity range - Sangyoon
    K1 = length(bI1); K2 = length(bI2);
    N1 = length(bI1e); N2 = length(bI2e);
%     LEN1 = 2^nextpow2(K1+N1-1);
%     LEN2 = 2^nextpow2(K2+N2-1);
%     ref_stack = zeros(K1+(N1-1),K2+(N2-1),numFrames-1); %zero padding is needed
%     ref_stack(N1:N1+K1-1,N2:N2+K2-1,:) = kym(bI1,bI2,1:numFrames-1);
%     cur_stack = zeros(K1+(N1-1),K2+(N2-1),numFrames-1); %zero padding is needed
%     cur_stack(1:N1,1:N2) = kym(bI1e,bI2e,2:numFrames);
%     
%     tic;
%     score_ffte = real(ifft2( fft2(ref_stack,LEN1,LEN2).*...
%         conj(fft2(cur_stack,LEN1,LEN2)) ));
%     % normalize
%     bNorm2 = ifft2((fft2(cur_stack)).*conj(fft2(cur_stack)));
%     
%     score_fft = score_ffte(1:K1+N1-1,1:K2+N2-1);
%     score_fft = score_fft/bNorm1./sqrt(bNorm2);
%     score_un = score_fft(K1:N1,K2:N2); % un-normalized
%     score_n = score_un(length(vP):-1:1,length(vF):-1:1); %counterindexing
%     toc;

    if numFrames == 2 % For image stack of 2, fft-based crosscorrelation is faster and make 
        % the correlation less sensitive to bright region in template window - Sangyoon
        score_nxc2 = normxcorr2(kym(bI1,bI2,1:numFrames-1),kym(bI1e,bI2e,2:numFrames));
        score = score_nxc2(K1:N1,K2:N2); % normalized
    else
        for j1 = 1:length(vP)
            v1 = vP(j1);
            for j2 = 1:length(vF)
                v2 = vF(j2);
                corrM = kym(bI1,bI2,1:numFrames-1).* ...
                    kym(bI1+v1,bI2+v2,2:numFrames);

                %The norm of the shifted image band at each frame.
                bNorm2 = sqrt(sum(sum(sum(kymP2(bI1+v1,bI2+v2,2:numFrames)))));

                %Normalize the correlation coefficients.
                score(j1,j2) = sum(corrM(:))/bNorm1/bNorm2;
            end
        end
    end
else
    kym       = reshape(kym,kymLen*kymWidth,numFrames);
    kymMask   = reshape(kymMask,kymLen*kymWidth,numFrames);
    
    numAvgImgs    = size(kymAvgImg,3);
    kymAvgImg     = reshape(kymAvgImg,kymLen*kymWidth,numAvgImgs);
    lastKymAvgImg = kymAvgImg(:,numAvgImgs);
    
    %Extend 'kymAvgImg' to have 'numFrames' columns by repeating 'lastKymAvgImg'.
    kymAvgImg = [kymAvgImg lastKymAvgImg*ones(1,numFrames-numAvgImgs)];
    %kymAvgImg = kymAvgImg(:)*ones(1,numFrames);
    
    %%%%% Debugging %%%%%%%%%%%%%%%
    %kymAvgImg(:) = 0;
    %%%%% Debugging %%%%%%%%%%%%%%%
    
    kym0   = (kym-kymAvgImg).*kymMask;
    %kym0P2 = kym0.*kym0;
    
    %We only consider frames whose texture area (after cutting off background
    % area) is bigger than 'bAreaThreshold'. Note: Background pixel values are
    % zero.
    %First, Get the linear index of the image block in the big cropped image.
    [BI1,BI2] = ndgrid(bI1,bI2);
    bI = (BI2(:)-1)*kymWidth+BI1(:);  
    allFrames = 1:numFrames-1;
    nzInd = arrayfun(@(x)sum(kymMask(bI,x)),allFrames);    
    validFrames = allFrames(nzInd >= bAreaThreshold);
    
    %If the number of valid frames is less than half of the number of
    % correlating frames, we reject the tracking for this point with the default
    % zero score.
    if length(validFrames) < numFrames/2;
        blockIsTooSmall = 1;
        return;
    end
    
    %We consider template shift in both the positive (to next frame) flow direction and
    % negative (from previous frame) flow direction.
    kymValid   = kym0(bI,validFrames);
    kymValidP2 = kymValid.*kymValid;
    bNorm      = sqrt(sum(kymValidP2,1));
    
    % SB:old score calculation function from imFlowTrack
    %     for j1 = 1:length(vP)
    %       v1 = vP(j1);
    %       for j2 = 1:length(vF)
    %          v2 = vF(j2);
    %          v = v2*kymWidth+v1;
    %
    %          kymShift = kym(bI+v,validFrames+1)-kymAvgImg(bI,validFrames);
    %          bNormS   = sqrt(sum(kymShift.^2.*kymMask(bI,validFrames),1));
    %
    %          corrM   = -ones(1,length(validFrames));
    %          nzInd   = find(bNorm~=0);
    %          nzInd   = nzInd(bNormS(nzInd)~=0);
    %          zeroInd = find(bNorm==0);
    %          zeroInd = zeroInd(bNormS(zeroInd)==0);
    %
    %          corrM(zeroInd) = 1;
    %          corrM(nzInd)   = sum(kymValid(:,nzInd).*kymShift(:,nzInd),1);
    %
    %          zeroInd = find(bNorm==0);
    %          if ~isempty(zeroInd)
    %             bNorm(zeroInd) = 1;
    %          end
    %          zeroInd = find(bNormS==0);
    %          if ~isempty(zeroInd)
    %             bNormS(zeroInd) = 1;
    %          end
    %
    %          score(j1,j2) = mean(corrM./bNorm./bNormS);
    %       end
    %     end
    
    % SB: beginning of vectorized score calculation
    [v1,v2]=ndgrid(vP,vF);
    v = v2*kymWidth+v1;    
    [bI2,v2]=ndgrid(bI,v(:));
    allbI=bI2+v2;
    
    % Create a matrix of size (size(bI,1)xsize(validFrames)xsize(v))
    kymShiftMatrix= reshape(kym(allbI(:),validFrames+1),...
        length(bI),numel(v),length(validFrames))-...
        permute(repmat(kymAvgImg(bI,validFrames),[1 1 numel(v)]),[1 3 2]);
    kymShiftMatrix= permute(kymShiftMatrix,[2 3 1]);
    bNormS=squeeze(sqrt(sum(kymShiftMatrix.^2.*...
        permute(repmat(kymMask(bI,validFrames),[1 1 numel(v)]),[3 2 1]),3)));
    validCorrM = squeeze(sum(kymShiftMatrix.* ...
        permute(repmat(kym0(bI,validFrames),[1 1 numel(v)]),[3 2 1]),3));
    
    % Initialize the correlation matrix
    corrM = -ones(numel(v),length(validFrames));
    % Set the correlation value of zero correlation elements to 1
    corrM(repmat(bNorm==0,numel(v),1) & bNormS==0)=1;
    nzInd = repmat(bNorm~=0,numel(v),1) & bNormS~=0;
    corrM(nzInd) = validCorrM(nzInd);
    
    % Set bNorm and bNormS null components to 1 (to avoid division by zero)
    bNorm(bNorm==0)=1;
    bNorm=repmat(bNorm,numel(v),1);
    bNormS(bNormS==0)=1;
    score = mean(corrM./bNorm./bNormS,2);
    score = reshape(score,size(v));
    
    minusOnesI = find(score(:)==-1);
    nMOnesI    = (score(:)~=-1);
    [~,minScoreI] = min(score(nMOnesI));
    ind = 1:length(score(:));
    ind([minusOnesI minScoreI]) = [];
    score([minusOnesI minScoreI]) = min(score(ind));
end
