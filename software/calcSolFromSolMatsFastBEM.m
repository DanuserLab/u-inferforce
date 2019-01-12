function [pos_f,force,sol_coef]=calcSolFromSolMatsFastBEM(M,sol_mats,x,y,ux,uy,forceMesh,L,x_out,y_out,varargin)
% Input check
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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
ip =inputParser;
ip.addRequired('M',@isnumeric);
ip.addRequired('sol_mats',@isstruct);
ip.addRequired('x',@isnumeric);
ip.addRequired('y',@isnumeric);
ip.addRequired('ux',@isnumeric);
ip.addRequired('uy',@isnumeric);
ip.addRequired('forceMesh',@isstruct);
ip.addRequired('L',@isscalar);
ip.addRequired('x_out',@isnumeric);%@(x)isscalar(x)||isempty(x));
ip.addRequired('y_out',@isnumeric);%@(x)isscalar(x)||isempty(x));
ip.addParamValue('LcurveDataPath','',@ischar);
ip.addParamValue('LcurveFigPath','',@ischar);
ip.addParamValue('LcurveFactor','',@isscalar);
ip.addParamValue('useLcurve',false,@islogical); % default assuming 34 um with 72 nm/pix resolution
ip.addParamValue('paxImg',[],@ismatrix);
ip.parse(M,sol_mats,x,y,ux,uy,forceMesh,L,x_out,y_out,varargin{:});
LcurveDataPath=ip.Results.LcurveDataPath;
LcurveFigPath=ip.Results.LcurveFigPath;
LcurveFactor=ip.Results.LcurveFactor;
paxImage = ip.Results.paxImg;
useLcurve = ip.Results.useLcurve;    

% ip =inputParser;
% ip.addParamValue('paxImg',[],@ismatrix);
% paxImage = ip.Results.paxImg;
if nargin <11
    paxImage =[];
    useLcurve = false;
elseif nargin <12
    useLcurve = false;
end
% If M has size mxn then u is a column vector of length m. u is the 
% displacement data on the mesh, that has been used to determine the
% forward map M.

%if no points of interests are specified, i.e x_out, y_out, then forces are 
%calculated on the nodes of the force mesh:                                                      
if nargin<6 || isempty(x_out)
    x_out=forceMesh.p(:,1);
    y_out=forceMesh.p(:,2);
end

[~, cols]=size(x_out);

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
    u=vertcat(ux_vec,uy_vec);
end
    
% See BEM_force_reconstruction for a nice explanation of the next
% steps:
if strcmp(sol_mats.tool,'svd')
    U=sol_mats.U;
    s=sol_mats.s;
    V=sol_mats.V;
    [sol_coef,~,~] = tikhonov(U,s,V,u,sqrt(L));
elseif strcmp(sol_mats.tool,'gsvd')
    % gSVD takes about twice as long as SVD
    U =sol_mats.U;
    sm=sol_mats.sm;
    X =sol_mats.X;
    [sol_coef,~,~] = tikhonov(U,sm,X,u,sqrt(L));
elseif strcmp(sol_mats.tool,'QR')
%     [normWeights]=getNormWeights(forceMesh);
%     eyeWeights = sol_mats.eyeWeights;
    sol_L =sol_mats.L;
    % check that regularization parameter and weights have not changed
    % (since Q,R have been calculated for a certain set of reg. par. and
    % weights!). But this should always be the case:
    if sol_L==L %&& sum(sol_nW~=normWeights)==0
        eyeWeights = sol_mats.eyeWeights;
        L = sol_mats.L;
        MpM = M'*M;
        if ~isempty(paxImage)
            paxWeights = getPaxWeights(forceMesh,paxImage,x_vec,y_vec,ux_vec,uy_vec);
            [Q,R] = qr((MpM+L*eyeWeights.*paxWeights));
            sol_coef=R\(Q'*(M'*u));
        else
            Q=sol_mats.Q;
            R=sol_mats.R;
            sol_coef=R\(Q'*(M'*u));
        end            
    else
        error('Weights or regularization parameter have been changed. QR cannot be reused!')
    end
elseif strcmp(sol_mats.tool,'1NormReg')
    eyeWeights =sol_mats.eyeWeights;
%     MpM=sol_mats.MpM;
    M=sol_mats.M;
    L=sol_mats.L;
    maxIter = sol_mats.maxIter;
    tolx = sol_mats.tolx;
    tolr = sol_mats.tolr;
    idxNonan = ~isnan(u);
    if any(~idxNonan)
        u = u(idxNonan);
        Mreal = M(idxNonan,:);
        MpM = Mreal'*Mreal;
    else
        MpM=M'*M;
        Mreal = M;
        clear M
    end

    sol_coef = iterativeL1Regularization(Mreal,MpM,u,eyeWeights,L,maxIter,tolx,tolr); 
elseif strcmpi(sol_mats.tool,'1NormRegLaplacian')
    % Now, perform the sparse deconvolution.
    Lap = sol_mats.Lap;
    % plot the solution for the corner
    MpM=sol_mats.MpM;
    maxIter = sol_mats.maxIter ;
    tolx = sol_mats.tolx;
    tolr = sol_mats.tolr;
    L=sol_mats.L;
%         [sol_coef,L] = calculateLfromLcurveSparse(M,MpM,u,Lap,maxIter,tolx,tolr,solMethodBEM);
    sol_coef = iterativeL1Regularization(M,MpM,u,L,-Lap,maxIter,tolx,tolr); %400=maximum iteration number
elseif strcmpi(sol_mats.tool,'LaplacianReg')
    % second order tikhonov regularization (including diagonal)
    % make Lap matrix
%         nBeads = round(size(M,1)/2);
    Lap = sol_mats.Lap;
    MpM=sol_mats.MpM;
    M=sol_mats.M;
    L=sol_mats.L;
    % For L-curve
%         [sol_coef,L] = calculateLfromLcurve(M,MpM,u,Lap,solMethodBEM);
    sol_coef=(-L*Lap+ MpM)\(M'*u);
elseif strcmp(sol_mats.tool,'backslash')
    % This matrix multiplication takes most of the time. Therefore we
    % store it for later use:
%     [normWeights]=getNormWeights(forceMesh);
%     eyeWeights =diag(normWeights);
    eyeWeights = sol_mats.eyeWeights;
    idxNonan = ~isnan(u);
    if any(~idxNonan) 
        u = u(idxNonan);
        Mreal = M(idxNonan,:);
        MpM = Mreal'*Mreal;
    else
        MpM=sol_mats.MpM;
        Mreal = M;
    end
    sol_coef=(L*eyeWeights+ MpM)\(Mreal'*u);
else
    error('No solution matrices have been found')
end

%Evaluation of the solution:
[fx, fy, x_out, y_out]=calcForcesFromCoef(forceMesh,sol_coef,x_out,y_out,'new');
pos_f=horzcat(x_out,y_out);
force=horzcat(   fx,   fy);

function [sol_coef,reg_corner] = calculateLfromLcurveSparse(L,M,MpM,u,eyeWeights,maxIter,tolx,tolr,LcurveDataPath,LcurveFigPath,LcurveFactor)
%examine a logarithmically spaced range of regularization parameters
alphas=10.^(log10(L)-0.5:1.25/LcurveFactor:log10(L)+0.5);
rho=zeros(length(alphas),1);
eta=zeros(length(alphas),1);
% eta0=zeros(length(alphas),1);
msparse=zeros(size(M,2),length(alphas));
% if matlabpool('size')==0
%     matlabpool open
% end
tolFactor = 20; % make the L-curve calculation faster with generous tolx with this factor
for i=1:length(alphas);
    disp(['testing L = ' num2str(alphas(i)) '... '])
    msparse(:,i)=iterativeL1Regularization(M,MpM,u,eyeWeights,alphas(i),maxIter,tolx*tolFactor,tolr);
    rho(i)=norm(M*msparse(:,i)-u);
    eta(i)=norm(msparse(:,i),1);
%     eta0(i)=sum(abs(msparse(:,i))>1);
end

% Find the L-corner
% [reg_corner,ireg_corner,~]=l_curve_corner(rho,eta,alphas);
save(LcurveDataPath,'rho','eta','alphas','L','msparse','-v7.3'); % saving before selection.
[reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L);

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

if mod(ireg_corner,1)>0 % if ireg_corner is interpolated
    disp(['L-corner regularization parmater L = ' num2str(reg_corner) '... final solution calculation ...'])
    sol_coef=iterativeL1Regularization(M,MpM,u,eyeWeights,reg_corner,maxIter,tolx,tolr);
else
    sol_coef = msparse(:,ireg_corner);
end
