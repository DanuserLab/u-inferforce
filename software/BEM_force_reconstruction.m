function [fx, fy, x_out, y_out, M, pos_u, u, sol_coef, sol_mats] = ...
    BEM_force_reconstruction(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,varargin)
% Synopsis  [fx fy x_out y_out M pos_u u sol_coef sol_mats] = BEM_force_reconstruction(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,meshPtsFwdSol,solMethodBEM)
%
% Input :  x,y,ux,uy have to be in the same units, namely
%         pixels. 
% Output: The output fx,fy is actually a surface stress with the same units
%         as the input E! In particular, the unit of the output force is
%         independent of the units of the input x,y,ux,uy.
%         The reason for this is essentially that the elastic stress is
%         only dependent on the non-dimensional strain which is given by
%         spatial derivatives of the displacements, that is du/dx. If u and
%         dx (essentially cluster_size) are in the same units, then the
%         resulting force has the same dimension as the input E.
%         u: is the measured displacement! (not the model u!)
% Achim Besser 2011
% Sangyoon Han 2013
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
ip.addRequired('x',@isnumeric);
ip.addRequired('y',@isnumeric);
ip.addRequired('ux',@isnumeric);
ip.addRequired('uy',@isnumeric);
ip.addRequired('forceMesh',@isstruct);
ip.addRequired('E',@isscalar);
ip.addRequired('L',@isscalar);
ip.addRequired('x_out',@isnumeric);%@(x)isscalar(x)||isempty(x));
ip.addRequired('y_out',@isnumeric);%@(x)isscalar(x)||isempty(x));
ip.addRequired('method',@(x)ischar(x)||isempty(x)); % updated in case for BEM
ip.addOptional('meshPtsFwdSol',[],@(x)isscalar(x) ||isempty(x));
ip.addOptional('solMethodBEM','QR',@ischar);
ip.addParamValue('basisClassTblPath','',@ischar);
ip.addParamValue('LcurveDataPath','',@ischar);
ip.addParameter('LcurveFigPath','',@ischar);
ip.addParamValue('LcurveFactor','',@isscalar);
ip.addParamValue('wtBar',-1,@isscalar);
ip.addParamValue('imgRows',[],@isscalar);
ip.addParamValue('imgCols',[],@isscalar);
ip.addParamValue('fwdMap',[],@isnumeric);
ip.addParamValue('thickness',472,@isscalar); % default assuming 34 um with 72 nm/pix resolution
ip.addParamValue('PoissonRatio',0.5,@isscalar); % default assuming 34 um with 72 nm/pix resolution
ip.addParamValue('useLcurve',false,@islogical); % default assuming 34 um with 72 nm/pix resolution
ip.addParamValue('paxImg',[],@ismatrix);
ip.addParamValue('strictBEM',false,@islogical); 
ip.addParamValue('lcornerOptimal','optimal',@ischar);
ip.parse(x,y,ux,uy,forceMesh,E,L,x_out,y_out,method,varargin{:});
meshPtsFwdSol=ip.Results.meshPtsFwdSol;
solMethodBEM=ip.Results.solMethodBEM;
basisClassTblPath=ip.Results.basisClassTblPath;
LcurveDataPath=ip.Results.LcurveDataPath;
LcurveFigPath=ip.Results.LcurveFigPath;
LcurveFactor=ip.Results.LcurveFactor;
wtBar=ip.Results.wtBar;
imgRows = ip.Results.imgRows;
imgCols = ip.Results.imgCols;
M = ip.Results.fwdMap;
lcornerOptimal = ip.Results.lcornerOptimal;
thickness = ip.Results.thickness;    
paxImage = ip.Results.paxImg;
useLcurve = ip.Results.useLcurve;    
strictBEM = ip.Results.strictBEM;    
v = ip.Results.PoissonRatio;

if nargin < 12 || isempty(solMethodBEM)
    solMethodBEM='QR';
end

[~, cols]=size(x);

if cols>1
    x_vec=reshape(x,[],1);
    y_vec=reshape(y,[],1);
    ux_vec=reshape(ux,[],1);
    uy_vec=reshape(uy,[],1);
    u=vertcat(ux_vec,uy_vec);
else
    x_vec=x;
    y_vec=y;
    ux_vec=ux;
    uy_vec=uy;
    u=vertcat(ux,uy);
end
pos_u=horzcat(x_vec,y_vec);

%construction of forward map, this takes a long time!

display('2.) Building up forward map:...');
inputFwdMap = false;
tic;
if nargin >= 10 && strcmp(method,'fast') && isempty(M) && ~strictBEM
    if isempty(imgRows) || isempty(imgCols)
        M=calcFwdMapFastBEM(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
            'basisClassTblPath',basisClassTblPath,'wtBar',wtBar,'PoissonRatio',v);
    else
%         M=calcFwdMapFastBEM(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
%             'basisClassTblPath',basisClassTblPath,'wtBar',wtBar,'imgRows',imgRows,'imgCols',imgCols,'thickness',thickness,'PoissonRatio',v);    
            % for old calcFwdMapFastBEM
            M=calcFwdMapFastBEM(x_vec, y_vec, forceMesh, E, meshPtsFwdSol,...
            'basisClassTblPath',basisClassTblPath,'wtBar',wtBar,'imgRows',imgRows,'imgCols',imgCols,'thickness',thickness);    
    end
elseif isempty(M)
    span = 1:length(forceMesh.bounds);
    M=calcFwdMap(x_vec, y_vec, forceMesh, E, span, meshPtsFwdSol/2^3,'fast');
else
    display('Using input Forward Map');
    inputFwdMap = true;
end
toc;
display('Done: forward map!');



% x = A\B is the solution to the equation Ax = B. The equation we have to
% solve is:
% (L*eyeWeights+ MpM)*sol_coef=(M'*u);

tic;
display('3.) Solve for coefficients, this is memory intensive [~5min]:... ')
if nargin >= 10 && strcmp(method,'fast')
    % If there are more than one basis function class, then correct
    % weighting of basis function according to their volume is of course
    % important when taking the norm!!! If there is only one basis function
    % class, then the weights are all one (e.g. square lattice where
    % boundary nodes are skipped) 
    % See refine_BEM_force_reconstruction for a nice explanation of the next
    % steps: <-Achim's comment
    
    % It is more correct to consider the weight to be Cholesky factor of /int hi(x) hj(x) dx
    
%     [normWeights,listNormWeights]=getNormWeights(forceMesh);
%     eyeWeights =diag(normWeights);
    
%     if length(listNormWeights)==1
        needGSVD =0;
%     else
%         needGSVD =1;
%     end
    % Checked (g)SVD against Matlab inversion! Found no advantage of
    % (g)SVD over Matlab inversion but (g)SVD is much slower. Differences
    % seem to arise only for irregular meshes or at mesh boundaries.
    % For size(M'*M)=7688*7688 I find:
    % M'*M\ =  19sec
    % csvd  = 340sec
    % cgsvd = 561sec
    % Therefore we force the Matlab back slash operator:
%     forceQR=1;
%     forceBackSlash=0;
    
    if ~needGSVD && strcmpi(solMethodBEM,'svd')
        [U,s,V] = csvd(M);
        [sol_coef,~,~] = tikhonov(U,s,V,u,sqrt(L));
        % store these matrices for next frames:
        sol_mats.U=U;
        sol_mats.s=s;
        sol_mats.V=V;
        sol_mats.tool='svd';
    elseif strcmpi(solMethodBEM,'svd') || strcmpi(solMethodBEM,'gsvd')
        % gSVD takes about twice as long as SVD
        [eyeWeights,~] =getGramMatrix(forceMesh);
        [U,sm,X,~] = cgsvd(M,eyeWeights);
        [sol_coef,~,~] = tikhonov(U,sm,X,u,sqrt(L));
        % store these matrices for next frames:
        sol_mats.U =U;
        sol_mats.sm=sm;
        sol_mats.X =X;
        sol_mats.tool='gsvd';
    elseif strcmpi(solMethodBEM,'QR')
        % for a force field with 2*6400 basis function, the residual
        % between the QR-sol and the sol obtained from the backslash
        % operator was: 2.0057e-06 for a mean force magnitude of
        % 85.7. Thus they seem to be numerical identical!
        % accounting for badly scaled linear system - Sangyoon 02/20/13
        % sol_coef = (M'*M+L*D^2)\(M'*u_sol); where D = scaling matrix
        % reference: p11 in Neumaier, Solving ill-conditioned and singular
        % linear systems: a tutorial on regularization
        [eyeWeights,~] =getGramMatrix(forceMesh);
        
        %Check for nan in u
        idxNonan = ~isnan(u);
        if any(~idxNonan)
            u = u(idxNonan);
            M = M(idxNonan,:);
            MpM = M'*M;
        else
            MpM=M'*M;
        end
        
        % make matrix with paxImage at the basis nodes
        if ~isempty(paxImage)
            paxWeights = getPaxWeights(forceMesh,paxImage,x_vec,y_vec,ux_vec,uy_vec);
            [Q,R] = qr((MpM+L*eyeWeights.*paxWeights));
            sol_coef=R\(Q'*(M'*u));
            sol_mats.eyeWeights=eyeWeights;
            sol_mats.L=L;
        else
            if useLcurve
                [~,L] = calculateLcurve(L,M,MpM,u,eyeWeights,LcurveDataPath,LcurveFigPath,LcurveFactor);
            end
            [Q,R] = qr((MpM+L*eyeWeights));
            sol_coef=R\(Q'*(M'*u));
            sol_mats.Q=Q;
            sol_mats.R=R;
            sol_mats.L=L;
            sol_mats.eyeWeights=eyeWeights;
        end            
        sol_mats.tool='QR';
    elseif strcmpi(solMethodBEM,'LaplacianReg')
        % second order tikhonov regularization (including diagonal)
        % make Lap matrix
        Lap = buildLaplacian(forceMesh);
        MpM=M'*M;
        % For L-curve
%         [sol_coef,L] = calculateLfromLcurve(M,MpM,u,Lap,solMethodBEM);
        if useLcurve
            [~,L] = calculateLcurve(L,M,MpM,u,-Lap,LcurveDataPath,LcurveFigPath,LcurveFactor);
        end
        sol_coef=(-L*Lap+ MpM)\(M'*u);
        % store these matrices for next frames:
        sol_mats.M=M;
        sol_mats.MpM=MpM;
        sol_mats.Lap=Lap;
        sol_mats.L=L;
        sol_mats.tool='LaplacianReg';
    elseif strcmpi(solMethodBEM,'1NormReg')
        % Now, perform the sparse deconvolution.
        disp('Performing sparse deconvolution; adoped from Aster et. al.')

        if strictBEM
            eyeWeights = eye(numel(forceMesh.base));
            tolx =  numel(forceMesh.base)*2.5e-6; % This will make tolx sensitive to overall number of nodes. (rationale: the more nodes are, 
        else
            [eyeWeights,~] =getGramMatrix(forceMesh); % possibly this is a
%         culprit that makes diagonalized force map. Now switching to
%         Achim's approach.
%         [normWeights]=getNormWeights(forceMesh);
%         eyeWeights =diag(normWeights);    
%             tolx =  forceMesh.numBasis*5e-6; % This will make tolx sensitive to overall number of nodes. (rationale: the more nodes are, 
%             % the larger tolerance should be, because misfit norm can be larger out of more nodes).
%             tolx =  sqrt(forceMesh.numBasis)*1e-3; % This will make tolx sensitive to overall number of nodes. (rationale: the more nodes are,the larger tolerance should be, because misfit norm can be larger out of more nodes). 
            % Filter out u by forceMesh boundary
%             bwstackImg=zeros(ceil(max(y)),ceil(max(x)));
%             bwstackImg(min(forceMesh.p(:,2)):max(forceMesh.p(:,2)),min(forceMesh.p(:,1)):max(forceMesh.p(:,1)))=1;
%             [insideIdx] = maskVectors(x,y,bwstackImg);
%             cur_ux=ux(insideIdx);
%             cur_uy=uy(insideIdx);
%             cur_u = [cur_ux; cur_uy];
%             cur_umag = (cur_u(:,1).^2+cur_u(:,1).^2).^0.5;
%             tolxScale=1;
%             tolxEstimate = tolxScale*quantile((cur_umag),0.95)*quantile((cur_umag),0.2)*(forceMesh.numBasis)^(49/60)/E; %empirical based on cur_u
%             if tolxEstimate<1e-5 % this is when u is almost none. We have to be generous
%                 tolx = 0.05;
%             elseif tolxEstimate<1e-3
%                 tolx = 0.03;
%             else
%                 tolx =  max(0.02,tolxEstimate); % based on u to estimate
%             end
            tolx = 0.2; % based on u to estimate
        end
        % plot the solution for the corner
        disp(['tolerance value: ' num2str(tolx)])
        %Check for nan in u
        idxNonan = ~isnan(u);
        if any(~idxNonan) && ~inputFwdMap
            u = u(idxNonan);
            Mreal = M(idxNonan,:);
            MpM = Mreal'*Mreal;
        elseif inputFwdMap
            u = u(idxNonan);
            Mreal = M(idxNonan,:);
            MpM = Mreal'*Mreal;
%             MpM=M'*M;
%             Mreal = M;
        else
            MpM=M'*M;
            Mreal = M;
       end
%         if any(~idxNonan)
%             u = u(idxNonan);
%             Mreal = M(idxNonan,:);
%             MpM = Mreal'*Mreal;
%         else
%             MpM=M'*M;
%             Mreal = M;
%         end
        maxIter = 50;
        tolr = 10;
        if useLcurve
            disp('L-curve ...')
            [sol_coef,L] = calculateLcurveSparse(L,Mreal,MpM,u,eyeWeights,maxIter,tolx,tolr,LcurveDataPath,LcurveFigPath,LcurveFactor,'lcornerOptimal',lcornerOptimal);
        else
            sol_coef = iterativeL1Regularization(Mreal,MpM,u,eyeWeights,L,maxIter,tolx,tolr); 
%             sol_coef = iterativeL1Regularization(Mreal,MpM,u,eyeWeights,L,maxIter,tolx,tolr,1,forceMesh); 
%             sol_coef = l1_ls(M,u,L,tolx); 
        end
%         sol_mats.nW=normWeights;
        sol_mats.eyeWeights=eyeWeights;
        sol_mats.L=L;
        sol_mats.M = M;
        sol_mats.MpM = MpM;
        sol_mats.maxIter = maxIter;
        sol_mats.tolx = tolx;
        sol_mats.tolr = tolr;
        sol_mats.tool='1NormReg';
    elseif strcmpi(solMethodBEM,'1NormRegLaplacian')
        % Now, perform the sparse deconvolution.
        disp('Performing sparse deconvolution with laplacian')

        disp('Building laplacian operator')
        Lap = buildLaplacian(forceMesh);
        % plot the solution for the corner
        %Check for nan in u
        idxNonan = ~isnan(u);
        if any(~idxNonan)
            u = u(idxNonan);
            M = M(idxNonan,:);
            MpM = M'*M;
        else
            MpM=M'*M;
        end

        maxIter = 10;
        tolx =  forceMesh.numBasis*3e-6; % This will make tolx sensitive to overall number of nodes. (rationale: the more nodes are, 
        % the larger tolerance should be, because misfit norm can be larger out of more nodes).
        disp(['tolerance value: ' num2str(tolx)])
        tolr = 1e-7;
%         [sol_coef,L] = calculateLfromLcurveSparse(M,MpM,u,Lap,maxIter,tolx,tolr,solMethodBEM);
        if useLcurve
            disp('L-curve ...')
            [sol_coef,L] = calculateLcurveSparse(L,M,MpM,u,-Lap,maxIter,tolx,tolr,LcurveDataPath,LcurveFigPath,LcurveFactor,'lcornerOptimal',lcornerOptimal);
        else
            sol_coef = iterativeL1Regularization(M,MpM,u,L,-Lap,maxIter,tolx,tolr); %400=maximum iteration number
        end
        sol_mats.L=L;
        sol_mats.Lap = Lap;
        sol_mats.M = M;
        sol_mats.MpM = MpM;
        sol_mats.maxIter = maxIter;
        sol_mats.tolx = tolx;
        sol_mats.tolr = tolr;
        sol_mats.tool='1NormRegLaplacian';
    elseif strcmpi(solMethodBEM,'backslash') || strcmpi(solMethodBEM,'\')
        % This matrix multiplication takes most of the time. Therefore we
        % store it for later use:
        [eyeWeights,~] =getGramMatrix(forceMesh);
        %Check for nan in u
        idxNonan = ~isnan(u);
        if any(~idxNonan) && ~inputFwdMap
            u = u(idxNonan);
            Mreal = M(idxNonan,:);
            MpM = Mreal'*Mreal;
            % for pos_u
            pos_u=reshape(pos_u,[],1);
            pos_u = pos_u(idxNonan);
            pos_u=reshape(pos_u,[],2);
        elseif inputFwdMap
            u = u(idxNonan);
            Mreal = M(idxNonan,:);
            MpM = Mreal'*Mreal;
            % for pos_u
            pos_u=reshape(pos_u,[],1);
            pos_u = pos_u(idxNonan);
            pos_u=reshape(pos_u,[],2);
        else
            MpM=M'*M;
            Mreal = M;
        end
        if useLcurve
            [sol_coef,L] = calculateLcurve(L,Mreal,MpM,u,eyeWeights,LcurveDataPath,LcurveFigPath,LcurveFactor,'lcornerOptimal',lcornerOptimal);
        else
            sol_coef=(L*eyeWeights+ MpM)\(Mreal'*u);
        end
        sol_mats.eyeWeights=eyeWeights;
        sol_mats.MpM=MpM;
        sol_mats.L=L;
        sol_mats.tool='backslash';

%         [fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,sol_coef,[],[],'new');
%         generateHeatmapFromGridData(x_out,y_out,fx,fy,'./backslash',0,1600)
%         sol_coef=(L*eyeWeights+ MpM)\(M'*u);
        % store these matrices for next frames:
    elseif strcmpi(solMethodBEM,'fourier')
        % Test Fourier-based solution
        % pad M to 2*(N-1,O-1) points to avoid wrap-around effects
        % beta_alpha = (M'*M+L*K^2)^(-1)*M'*D
        y_out = forceMesh.p(:,2);
        x_out = forceMesh.p(:,1);
        x_vec = pos_u(:,1);
        y_vec = pos_u(:,2);
        ux_vec = u(1:length(x_vec),1);
        uy_vec = u(length(x_vec)+1:end,1);
        colsOut = sum(y_out==y_out(1));
        reg_grid(:,:,1) = reshape(x_out,[],colsOut);
        reg_grid(:,:,2) = reshape(y_out,[],colsOut);
        [grid_mat,iu_mat, i_max, j_max] = interp_vec2grid([x_vec, y_vec], [ux_vec, uy_vec],[],reg_grid);
        iuvec = iu_mat(:); % make it 1D

        % I am trying to invert M in frequency domain.
        % generate spectra
        rowsOut = size(grid_mat,1);
        
        Mgrid = calcFwdMapFastBEM(grid_mat(:,:,1),grid_mat(:,:,2), forceMesh, E, meshPtsFwdSol,'basisClassTblPath',basisClassTblPath,'PoissonRatio',v);
        Mgridfirst = Mgrid(:,1); % this might need to be shifted back some how.
        % scaling (by grid spacing) and shifting (grid spacing) should be fixed
%         MgridfirstMax = max(abs(Mgridfirst));
        gridSpacing = grid_mat(1,2,1) - grid_mat(1,1,1);
        MgridfirstNor = Mgridfirst*gridSpacing;%/MgridfirstMax;
%         MgridMid = Mgrid(:,colsOut/2*rowsOut+rowsOut/2);
        % zero padding to avoid wrap-around effect of circular convolution
        Nu = length(iuvec);
        N_G=2*Nu;
        % These might not be a power of 2, make sure that they are:
        N_pad=pow2(nextpow2(N_G));
        
        iuvecPadded=padarray(iuvec,N_pad-Nu,0,'post');
        MgridfirstPadded  = padarray(MgridfirstNor,N_pad-Nu,0,'post');
        
        uspec = fft(iuvecPadded);
        Mspec = fft(MgridfirstPadded);
        % regularization
        Fspec = conj(Mspec)./(conj(Mspec).*Mspec+L*ones(size(Mspec))).*uspec;
        sol_coef_fftPadded = ifft(Fspec,'symmetric');
        sol_coef_fft = sol_coef_fftPadded(1:Nu);
        [fxFT,fyFT,x_outFT,y_outFT]=calcForcesFromCoef(forceMesh,sol_coef_fft,x_out,y_out,'new');
        generateHeatmapFromGridData(x_outFT,y_outFT,fxFT,fyFT,'./fft',0,1600);
        
        sol_mats.eyeWeights=eyeWeights;
        sol_mats.MpM=MpM;
        sol_mats.L=L;
        sol_mats.tool='fourier';
    else
        error(['I don''t understand the input for the solution method: ',solMethodBEM])
    end
    % Here we use the identity matrix (all basis classes have equal weight):
    % sol_coef=(L*eye(2*forceMesh.numBasis)+M'*M)\(M'*u);
else
    % normalization of basis function will be important when taking the norm!!!
    % This has not been considered yet!
%     [normWeights,~]=getNormWeights(forceMesh);
%     eyeWeights =diag(normWeights);

    MpM=M'*M;
    sol_coef=(L*eye(2*forceMesh.numNodes)+MpM)\(M'*u);
%     sol_coef=(L*eyeWeights+MpM)\(M'*u);
    % store these matrices for next frames:
    sol_mats.MpM=MpM;
    sol_mats.tool='backslash';
end
toc;
display('Done: coefficients!');


%if no points of interests are specified, i.e x_out, y_out, then forces are 
%calculated on the nodes of the force mesh:                                                      
if nargin<9 || isempty(x_out) && ~strictBEM
    x_out=forceMesh.p(:,1);
    y_out=forceMesh.p(:,2);
else
    x_min=min(forceMesh.p(:,1));
    x_max=max(forceMesh.p(:,1));
    y_min=min(forceMesh.p(:,2));
    y_max=max(forceMesh.p(:,2));
    x_out_vec = x_min:x_max;
    y_out_vec = y_min:y_max;
    [x_out,y_out] = meshgrid(x_out_vec,y_out_vec);
    x_out = x_out(:);
    y_out = y_out(:);
end

%Evaluation of the solution:
display('4.) Evaluate solution:... ')
tic;
if nargin >= 10 && strcmp(method,'fast') && ~strictBEM
    [fx,fy,x_out,y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,'new');
else
    fx = zeros(size(x_out));
    fy = zeros(size(y_out));
    t=cputime;
    for j=1:2*forceMesh.numNodes
        fx = fx+sol_coef(j)*forceMesh.base(j).f_intp_x(x_out,y_out);
        fy = fy+sol_coef(j)*forceMesh.base(j).f_intp_y(x_out,y_out);
        if mod(j,20)==0
            display([num2str(j) 'th basis function was processed so far ... time passed: ' num2str(cputime-t)])
            t=cputime;
        end
    end
end
toc;
display('Done: solution!')

function Lap = buildLaplacian(forceMesh)
nBasis = forceMesh.numBasis;
basisx=zeros(nBasis,1);
basisy=zeros(nBasis,1);
for k=1:nBasis
    basisx(k) = forceMesh.basis(k).node(1);
    basisy(k) = forceMesh.basis(k).node(2);
end
nBasisx = size(unique(basisx),1);
nBasisy = size(unique(basisy),1);

% this is for assuring the boundary to be considered for being
% penalized for regularization
display('Building Laplacian Map...')
tic;
Lap = zeros(2*nBasis,2*nBasis);
k=1;
for ii=1:nBasisx
%     tempLapx = zeros(nBasisy,nBasis); %for parfor constraints
%     tempLapy = zeros(nBasisy,nBasis); %for parfor constraints
    for jj=1:nBasisy
        Lap2D = zeros(nBasisy,nBasisx);
        if ii==1 || jj==1 || ii==nBasisx || jj==nBasisy
            % this can be improved with basis function convolution
            if ii==1 
                if jj==1
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj+1,ii)=1;
                    Lap2D(jj,ii+1)=1;
                elseif jj==nBasisy
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj-1,ii)=1;
                    Lap2D(jj,ii+1)=1;
                else
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj,ii+1)=1;
                    Lap2D(jj+1,ii)=1;
                end
            elseif ii==nBasisx
                if jj==1
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj+1,ii)=1;
                    Lap2D(jj,ii-1)=1;
                elseif jj==nBasisy
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj,ii-1)=1;
                    Lap2D(jj-1,ii)=1;
                else
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj,ii-1)=1;
                    Lap2D(jj-1,ii)=1;
                end
            elseif jj==1
                if ii~=nBasisx && ii~=1
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj+1,ii)=1;
                    Lap2D(jj,ii+1)=1;
                end
            elseif jj==nBasisy                        
                if ii~=nBasisx && ii~=1
                    Lap2D(jj,ii)=-2;
                    Lap2D(jj-1,ii)=1;
                    Lap2D(jj,ii+1)=1;
                end
            end                        
        else
            % diagonal laplacian
            Lap2D(jj,ii) = -6;
            Lap2D(jj,ii+1) = 1;
            Lap2D(jj,ii-1) = 1;
            Lap2D(jj+1,ii) = 1;
            Lap2D(jj-1,ii) = 1;
            Lap2D(jj+1,ii+1) = 0.5;
            Lap2D(jj-1,ii+1) = 0.5;
            Lap2D(jj+1,ii-1) = 0.5;
            Lap2D(jj-1,ii-1) = 0.5;
        end
        Lap(k,1:nBasis) = reshape(Lap2D,nBasis,1)';
        Lap(k+nBasis,nBasis+1:2*nBasis) = reshape(Lap2D,nBasis,1)';
%         Lap(k+nBasis,1:nBasis) = reshape(Lap2D,nBasis,1)';
%         Lap(k,nBasis+1:2*nBasis) = reshape(Lap2D,nBasis,1)';
        k=k+1;
%         tempLapx(jj,:) = reshape(Lap2D,nBasis,1)';
%         tempLapy(jj,:) = reshape(Lap2D,nBasis,1)';
    end
%     Lap((ii-1)*nBasisy+1:(ii-1)*nBasisy+nBasisy,1:nBasis) = tempLapx;
%     Lap((ii-1)*nBasisy+nBasis+1:(ii-1)*nBasisy+nBasis+nBasisy,nBasis+1:2*nBasis) = tempLapy;
end
toc

function [sol_coef,reg_corner] = calculateLcurve(L,M,MpM,u,eyeWeights,LcurveDataPath,LcurveFigPath,LcurveFactor,varargin)
% Input check
ip =inputParser;
ip.addParamValue('lcornerOptimal','optimal',@ischar);
ip.parse(varargin{:});
lcornerOptimal = ip.Results.lcornerOptimal;

%examine a logarithmically spaced range of regularization parameters
alphas=10.^(log10(L)-4.5:1.25/LcurveFactor:log10(L)+2.5);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
mtik=zeros(size(M,2),length(alphas));

for i=1:length(alphas);
    mtik(:,i)=(MpM+alphas(i)*eyeWeights)\(M'*u);
    rho(i)=norm(M*mtik(:,i)-u);
    eta(i)=norm(mtik(:,i));
end
% Find the corner of the Tikhonov L-curve
try
    if strcmp(lcornerOptimal,'optimal')
        disp('Inflection point smaller than L-corner will be chosen')
        [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L,'inflection',2);
        if isempty(reg_corner)
            disp('Inflection point larger than L-corner was not identified. L-corner will be chosen')
            [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L,'inflection',0); %L-corner
        end
    elseif strcmp(lcornerOptimal,'lcorner')
        disp('L-corner will be chosen')
        [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L,'inflection',0); %L-corner
    end
    %     ireg_corner=[];
    %     [reg_corner,rhoC,etaC]=l_corner(rho,eta,alphas);
catch
    ireg_corner=[];
    [reg_corner,rhoC,etaC]=l_corner(rho,eta,alphas);
end
    
% Plot the sparse deconvolution L-curve.
hLcurve = figure;
set(hLcurve, 'Position', [50 300 200 200])

loglog(rho,eta,'k-');
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Norm ||Lm||_{2}');
hold on
% mark and label the corner
if ~isempty(ireg_corner) && mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    rho_corner = rho(floor(ireg_corner))+mod(ireg_corner,1)*(rho(floor(ireg_corner)+1)-rho(floor(ireg_corner)));
    eta_corner = eta(floor(ireg_corner))+mod(ireg_corner,1)*(eta(floor(ireg_corner)+1)-eta(floor(ireg_corner)));
elseif ~isempty(ireg_corner) 
    rho_corner = rho(ireg_corner);
    eta_corner = eta(ireg_corner);
else
    rho_corner = rhoC;
    eta_corner = etaC;
end    
H=loglog(rho_corner,eta_corner,'ro');
set(H,'markersize',6)
H=text(rho_corner,1.1*eta_corner,...
    ['    ',num2str(reg_corner,'%5.1e')]);
set(H,'Fontsize',7);
% axis([1e-2 100 0.001 1e8])
disp('Printing L-curve...')
% print -deps2 nameSave
% print(hLcurve,'Lcurve.eps','-depsc')
saveas(hLcurve,LcurveFigPath);
close(hLcurve)

save(LcurveDataPath,'rho','eta','reg_corner','ireg_corner','alphas','mtik','eyeWeights','-v7.3');

if isempty(ireg_corner) || mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    disp(['L-corner regularization parmater L = ' num2str(reg_corner) '... final solution calculation ...'])
    sol_coef=(MpM+reg_corner*eyeWeights)\(M'*u);
else
    sol_coef = mtik(:,ireg_corner);
end

%old code
% display('Calculating L-curve ...')
% b  = M'*u;
% bfSigmaRange2 = L*[1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5];
% bfSigmaRange1 = bfSigmaRange2/5;
% bfSigmaRange3 = bfSigmaRange2*5;
% bfSigmaCandidate = [bfSigmaRange1 bfSigmaRange2 bfSigmaRange3];
% bfSigmaCandidate = sort(bfSigmaCandidate);
% coefNorm    = zeros(1,length(bfSigmaCandidate));
% residueNorm = zeros(1,length(bfSigmaCandidate));
% Lap2 = Lap'*Lap;
% for kk = 1:length(bfSigmaCandidate)
%  B_i = MpM + bfSigmaCandidate(kk)*Lap2;
%  coef_i = B_i\b;
%  coefNorm(kk)    = sqrt(sum((Lap*coef_i).^2)); % Solution norm ||Lm||
%  residueNorm(kk) = sqrt(sum((M*coef_i-u).^2)); % Residual norm ||Gm-d||
% end
% 
% %Plot the L-curve.
% figure; hold on;
% plot(residueNorm,coefNorm,'.'); xlabel('Residual norm ||M*coef-u||'); 
% ylabel('Solution seminorm ||L*coef||');
% for kk = 1:length(bfSigmaCandidate)
%  text(residueNorm(kk),coefNorm(kk),num2str(bfSigmaCandidate(kk)));
% end
% 
% options.WindowStyle='normal';
% answer = inputdlg('Please identify the corner:','Input for corner',1,{num2str(L)},options);
% Lout = str2double(answer{1});

function [sol_coef,reg_corner] = calculateLcurveSparse(L,M,MpM,u,eyeWeights,maxIter,tolx,tolr,LcurveDataPath,LcurveFigPath,LcurveFactor,varargin)
% Input check
ip =inputParser;
ip.addParamValue('lcornerOptimal','optimal',@ischar);
ip.parse(varargin{:});
lcornerOptimal = ip.Results.lcornerOptimal;
%examine a logarithmically spaced range of regularization parameters
alphas=10.^(log10(L)-2.5:1.25/LcurveFactor:log10(L)+2.5);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
% eta0=zeros(length(alphas),1);
msparse=zeros(size(M,2),length(alphas));
% if matlabpool('size')==0
%     matlabpool open
% end
tolFactor = 1; % make the L-curve calculation faster with generous tolx with this factor
for i=1:length(alphas);
    disp(['testing L = ' num2str(alphas(i)) '... '])
    msparse(:,i)=iterativeL1Regularization(M,MpM,u,eyeWeights,alphas(i),maxIter,tolx*tolFactor,tolr);
%     msparse(:,i) = l1_ls(M,u,alphas(i),tolx);
    rho(i)=norm(M*msparse(:,i)-u);
    eta(i)=norm(msparse(:,i),1);
%     eta0(i)=sum(abs(msparse(:,i))>1);
end

% Find the L-corner
% [reg_corner,ireg_corner,~]=l_curve_corner(rho,eta,alphas);
save(LcurveDataPath,'rho','eta','alphas','L','msparse','-v7.3'); % saving before selection.
if strcmp(lcornerOptimal,'optimal')
    disp('Inflection point larger than L-corner will be chosen')
    [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L,'inflection',1); %inflection point larger than l-corner
    if isempty(reg_corner)
        disp('Inflection point larger than L-corner was not identified. L-corner will be chosen')
        [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L,'inflection',0); %L-corner
    %     ireg_corner=[];
    %     [reg_corner,rhoC,etaC]=l_corner(rho,eta,alphas);
    end    
elseif strcmp(lcornerOptimal,'lcorner')
    disp('Inflection point larger than L-corner was not identified. L-corner will be chosen')
    [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L,'inflection',0); %L-corner    
end
% Also, I can use L0 norm information to choose regularization parameter

% Plot the sparse deconvolution L-curve.
hLcurve = figure;
set(hLcurve, 'Position', [100 100 500 500])

loglog(rho,eta,'k-');
ylim([min(eta) max(eta)])
xlabel('Residual Norm ||Gm-d||_{2}');
ylabel('Solution Norm ||m||_{1}');
hold on
% mark and label the corner
if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    rho_corner = rho(floor(ireg_corner))+mod(ireg_corner,1)*(rho(floor(ireg_corner)+1)-rho(floor(ireg_corner)));
    eta_corner = eta(floor(ireg_corner))+mod(ireg_corner,1)*(eta(floor(ireg_corner)+1)-eta(floor(ireg_corner)));
else
    rho_corner = rho(ireg_corner);
    eta_corner = eta(ireg_corner);
end    
H=loglog(rho_corner,eta_corner,'ro');
set(H,'markersize',6)
H=text(rho_corner,1.1*eta_corner,...
    ['    ',num2str(reg_corner,'%5.1e')]);
set(H,'Fontsize',7);
% axis([1e-2 100 0.001 1e8])
disp('Displaying the 1-norm L-curve')
% print -deps2 nameSave
% print(hLcurve,strcat(nameSave,'.eps'),'-depsc')
saveas(hLcurve,LcurveFigPath);
save(LcurveDataPath,'rho','eta','reg_corner','ireg_corner','alphas','rho_corner','eta_corner','msparse','-v7.3');
close(hLcurve)

if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    disp(['L-corner regularization parmater L = ' num2str(reg_corner) '... final solution calculation ...'])
    sol_coef=iterativeL1Regularization(M,MpM,u,eyeWeights,reg_corner,maxIter,tolx,tolr);
else
    sol_coef = msparse(:,ireg_corner);
end


