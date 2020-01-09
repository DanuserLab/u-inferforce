function efficientSubPixelRegistration(movieData, varargin)
%EFFICIENTSUBPIXELREGISTRATION Method for registering MovieData frames to correct for stage drift  
%
% efficientSubPixelRegistration 
%
% SYNOPSIS Process Wrapper to execute registering of MovieData frames to correct for stage drift
%          based. Core algorithm/code by :% [1] Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
%          "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
%          156-158 (2008). (s)
%
% INPUT   
%   MovieData - A MovieData object describing the movie to be processed
%                    
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   
%  
%   
% See also: dftregistration.m 
%
% Andrew R. Jamieson Feb. 2017
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


% ----------- Input ----------- %%
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData, varargin{:});
paramsIn=ip.Results.paramsIn;


%Get the indices of any previous stage drift processes  
[movieData1, thisProcess, iProc] = getOwnerAndProcess(movieData,'EfficientSubpixelRegistrationProcess',true);

assert(movieData==movieData1);

%Parse input, store in parameter structure
p = parseProcessParams(thisProcess, paramsIn);

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
    wtBar = waitbar(0,'Initializing...','Name', thisProcess.getName());
else
    wtBar = -1;
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
thisProcess.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(3,numel(movieData.channels_));
mkClrDir(p.OutputDirectory);
for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i)];
    mkClrDir(outFilePaths{1,i});
end


% Loading reference channel images and bead image stack
if ~isempty(p.referenceFramePath) 
  [~,refName,refExt] = fileparts(thisProcess.funParams_.referenceFramePath);
else
  refName = ['refernceFrame' num2str(p.referenceFrameNum)];
  refExt  ='.tiff';
  % [~,refName,refExt] = fileparts(thisProcess.funParams_.referenceFramePath);
  % refFrame = double(imread(p.referenceFramePath));
end

outFilePaths{2,p.ChannelIndex(1)} = [p.OutputDirectory filesep refName refExt];
outFilePaths{3,p.ChannelIndex(1)} = [p.OutputDirectory filesep 'transformationParameters.mat'];

thisProcess.setOutFilePaths(outFilePaths);


%% --------------- Stage drift correction ---------------%%% 

disp('Starting correcting stage drift [EfficientSubpixelRegistrationbyCrossCorrelation]...')

% Anonymous functions for reading input/output
outFile=@(chan,frame) [outFilePaths{1,chan} filesep imageFileNames{chan}{frame}];


ImStack = zeros([movieData.imSize_ nFrames]);
beadsChannel = movieData.channels_(p.ChannelIndex(1));
for j = 1:nFrames, ImStack(:,:,j) = double(beadsChannel.loadImage(j)); end


% Loading reference channel images and bead image stack
if ~isempty(p.referenceFramePath) 
  refFrame = double(imread(p.referenceFramePath));
else
  disp(['Using frame ' num2str(p.referenceFrameNum) ' as reference']);
  refFrame = ImStack(:,:,p.referenceFrameNum);
  % refFrame = double(imread(p.referenceFramePath));
end


disp('Calculating subpixel-wise registration...')
logMsg = @(t) ['Performing sub-pixel registration on frame: ' num2str(t)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;
 % Perform sub-pixel registration
if ishandle(wtBar), waitbar(0,wtBar,sprintf(logMsg(1))); end


nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;
Tout = zeros(nFrames, 4);

for frame_num = 1:nFrames

   moving = ImStack(:,:,frame_num);

   % Core Algorithm
   [DFT_output, Greg_beads] = dftregistration(fft2(refFrame), fft2(moving), p.usfac);
   DFTout(frame_num) = DFT_output;
   disp(['-Frame: ' num2str(frame_num)]);
   disp(['x shift: ' num2str(DFTout(frame_num).row_shift)]);
   disp(['y shift: ' num2str(DFTout(frame_num).col_shift)]);
   I_beads = abs(ifft2(Greg_beads));
   imwrite(uint16(I_beads), outFile(1, frame_num));

   % Apply transform to each selected channel
   for i = 2:numel(p.ChannelIndex)

      % Load frame from appropriate channel
      iChan = p.ChannelIndex(i);
      imChan = double(movieData.channels_(iChan).loadImage(frame_num));
      % disp('Results saved under:');
      % disp(outFilePaths{1, iChan});

      % Apply DFT-based transformation
      fft_imChan = DFT_apply(fft2(imChan), DFT_output);
      I2 = abs(ifft2(fft_imChan));
      
      imwrite(uint16(I2), outFile(iChan, frame_num));

   end
% Update the waitbar
%      if mod(frame_num, 5)==1 && ishandle(wtBar)
     if ishandle(wtBar)
         tj = toc;
    %      nj = (i-1)*nFrames+ frame_num;
         nj = frame_num;
         waitbar(nj/nFrames, wtBar, sprintf([logMsg(frame_num) timeMsg(tj*nTot/nj-tj)]));
     end
end


% Loading reference channel images and bead image stack
if ~isempty(p.referenceFramePath) 
  refFrame = double(imread(p.referenceFramePath));
else
  disp(['Using frame ' num2str(p.referenceFrameNum) ' as reference']);
  refFrame = ImStack(:,:,p.referenceFrameNum);
  % refFrame = double(imread(p.referenceFramePath));
end

imwrite(uint16(refFrame), outFilePaths{2, p.ChannelIndex(1)});

T = [DFTout.row_shift; DFTout.col_shift]';
save(outFilePaths{3, p.ChannelIndex(1)},'DFTout', 'T');
% Close waitbar
if ishandle(wtBar), close(wtBar); end
disp('Finished correcting stage drift!')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output, Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007
%
% Rewrote all code not authored by either Manuel Guizar or Jim Fienup
% Manuel Guizar - May 13, 2016
%
% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).
%
% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)
%
% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.
%
%
% Copyright (c) 2016, Manuel Guizar Sicairos, James R. Fienup, University of Rochester
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Rochester nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if ~exist('usfac','var')
    usfac = 1;
end

[nr,nc]=size(buf2ft);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);

if usfac == 0
    % Simple computation of error and phase difference without registration
    CCmax = sum(buf1ft(:).*conj(buf2ft(:)));
    row_shift = 0;
    col_shift = 0;
elseif usfac == 1
    % Single pixel registration
    CC = ifft2(buf1ft.*conj(buf2ft));
    CCabs = abs(CC);
    [row_shift, col_shift] = find(CCabs == max(CCabs(:)));
    CCmax = CC(row_shift,col_shift)*nr*nc;
    % Now change shifts so that they represent relative shifts and not indices
    row_shift = Nr(row_shift);
    col_shift = Nc(col_shift);
elseif usfac > 1
    % Start with usfac == 2
    CC = ifft2(FTpad(buf1ft.*conj(buf2ft),[2*nr,2*nc]));
    CCabs = abs(CC);
    [row_shift, col_shift] = find(CCabs == max(CCabs(:)),1,'first');
    CCmax = CC(row_shift,col_shift)*nr*nc;
    % Now change shifts so that they represent relative shifts and not indices
    Nr2 = ifftshift(-fix(nr):ceil(nr)-1);
    Nc2 = ifftshift(-fix(nc):ceil(nc)-1);
    row_shift = Nr2(row_shift)/2;
    col_shift = Nc2(col_shift)/2;
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;     
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac));
        % Locate maximum and map back to original pixel grid 
        CCabs = abs(CC);
        [rloc, cloc] = find(CCabs == max(CCabs(:)),1,'first');
        CCmax = CC(rloc,cloc);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;    
    end

    % If its only one row or column the shift along that dimension has no
    % effect. Set to zero.
    if nr == 1,
        row_shift = 0;
    end
    if nc == 1,
        col_shift = 0;
    end
    
end  

rg00 = sum(abs(buf1ft(:)).^2);
rf00 = sum(abs(buf2ft(:)).^2);
error = 1.0 - abs(CCmax).^2/(rg00*rf00);
error = sqrt(abs(error));
diffphase = angle(CCmax);

% output=[error,diffphase,row_shift,col_shift, Nr, Nc, nr, nc];

output.error = error;
output.diffphase = diffphase;
output.row_shift = row_shift;
output.col_shift = col_shift;
output.Nr = Nr;
output.Nc = Nc;
output.nr = nr;
output.nc = nc;
output.usfac = usfac;

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(1i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(1i*diffphase);
end
return


function Greg = DFT_apply(buf2ft, p)
     
   % Compute registered version of buf2ft
   if p.usfac > 0
       [p.Nc, p.Nr] = meshgrid(p.Nc, p.Nr);
       Greg = buf2ft.*exp(1i*2*pi*(-p.row_shift*p.Nr/p.nr-p.col_shift*p.Nc/p.nc));
       Greg = Greg*exp(1i*p.diffphase);
   elseif p.usfac == 0
       Greg = buf2ft*exp(1i*p.diffphase);
   end


function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff', 'var')~=1, roff=0;  end
if exist('coff', 'var')~=1, coff=0;  end
if exist('usfac','var')~=1, usfac=1; end
if exist('noc',  'var')~=1, noc=nc;  end
if exist('nor',  'var')~=1, nor=nr;  end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-1i*2*pi/(nc*usfac))*( ifftshift(0:nc-1).' - floor(nc/2) )*( (0:noc-1) - coff ));
kernr=exp((-1i*2*pi/(nr*usfac))*( (0:nor-1).' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return


function [ imFTout ] = FTpad(imFT,outsize)
% imFTout = FTpad(imFT,outsize)
% Pads or crops the Fourier transform to the desired ouput size. Taking 
% care that the zero frequency is put in the correct place for the output
% for subsequent FT or IFT. Can be used for Fourier transform based
% interpolation, i.e. dirichlet kernel interpolation. 
%
%   Inputs
% imFT      - Input complex array with DC in [1,1]
% outsize   - Output size of array [ny nx] 
%
%   Outputs
% imout   - Output complex image with DC in [1,1]
% Manuel Guizar - 2014.06.02

if ~ismatrix(imFT)
    error('Maximum number of array dimensions is 2')
end
Nout = outsize;
Nin = size(imFT);
imFT = fftshift(imFT);
center = floor(size(imFT)/2)+1;

imFTout = zeros(outsize);
centerout = floor(size(imFTout)/2)+1;

% imout(centerout(1)+[1:Nin(1)]-center(1),centerout(2)+[1:Nin(2)]-center(2)) ...
%     = imFT;
cenout_cen = centerout - center;
imFTout(max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)),max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2))) ...
    = imFT(max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1),Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)));

imFTout = ifftshift(imFTout)*Nout(1)*Nout(2)/(Nin(1)*Nin(2));
return


