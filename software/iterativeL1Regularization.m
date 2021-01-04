function mreg=iterativeL1Regularization(G,GTG,d,L,alpha,maxiter,tolx,tolr,doPlot,forceMesh)
% mreg=iterativeL1Regularization(G,GTG,d,L,alpha,maxiter,tolx,tolr,m_diff)
% solves L1 regularization problem with forward matrix G, displacement
% vector d, regularization parameter alpha, semi-norm matrix L.
% tolx and m_diff are used for convergence criteria
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

% Default for tolr=1.0e-6
if (nargin < 9)
  doPlot=false; % plot all iterations
end

if (nargin < 8)
  tolr=1.0e-6;
end

% Default for tolx=1.0e-4;
if (nargin < 7)
  tolx=1.0e-3;
end

% Default for maxiter=100
if (nargin < 6)
  maxiter=100;
end

%% Use gpuArray if possible to accelerate inversion process (assisted by Yi Du 20150629)
nGPU = gpuDeviceCount;
usedGPU = false;
if nGPU>0
    gD = gpuDevice();
    memGTG = whos('GTG');
    memL = whos('L');
    if gD.AvailableMemory>2*(memGTG.bytes + memL.bytes)
%     Ggpu=gpuArray(G);
%     dgpu=gpuArray(d);
%     Lgpu = gpuArray(L);
    % Start with an initial unweighted solution.
        iter=1;
        display(['L1 Norm regularization (iteration:' num2str(iter) ')....'])
        tic
        % process this step with CPU owing to memory limits
        Agpu = gpuArray(2*GTG+alpha*(L'*L));
        Bgpu = gpuArray(2*G'*d);
        m=Agpu\Bgpu;
        toc
        while (iter < maxiter)
            mold=gather(m);

            iter=iter+1;

            % get get the magnitude of Lm, but don't let any element be less than tolr
            absLm=abs(L*mold);
            absLm(absLm<tolr)=tolr;

            % build the diagonal weighting matrix for this iteration
            Rgpu=diag(1./absLm); % R is also a gpuArray
            display(['(iteration:' num2str(iter) ')'])

            tic
            R = gather(Rgpu);
            Agpu = gpuArray(2*GTG+alpha*(L'*R*L));
            m=Agpu\Bgpu;
            m=gather(m);
            clear Agup;
            toc

            if (norm(m-mold)/(1+norm(mold)) < tolx) %|| norm(m-mold)<m_diff
                mreg=m;
                display(['norm(m-mold)=' num2str(norm(m-mold)) ', 1+norm(mold)=' num2str(1+norm(mold)) ', norm(m-mold)/(1+norm(mold))=' ...
                  num2str(norm(m-mold)/(1+norm(mold)))])
                return
            else
                display(['norm(m-mold)=' num2str(norm(m-mold)) ', 1+norm(mold)=' num2str(1+norm(mold)) ', norm(m-mold)/(1+norm(mold))=' ...
                num2str(norm(m-mold)/(1+norm(mold)))])
            end
        end   
        usedGPU=true;
    end
end
if ~usedGPU
    % unchanging constants in the system that is solved repeatedly
    % GTG=G'*G;
    GTd=G'*d;
    % Start with an initial unweighted solution.
    iter=1;
    display(['L1 Norm regularization (iteration:' num2str(iter) ')....'])
    tic
    m=(2*GTG+alpha*(L'*L))\(2*G'*d);
    toc
    % iterate until maxiter or we converge and return
    sparFac = 1;
    % accFac = 0.9;
    while (iter < maxiter)
      iter=iter+1;

      % get get the magnitude of Lm, but don't let any element be less than tolr
      absLm=abs(L*m);
      if doPlot
          height = sum(forceMesh.p(:,1)==forceMesh.p(1,1));
          figure, imshow(reshape(m,height,[]),[])
          colorbar
          colormap jet
          title(['iteration number=' num2str(iter-1)])
    %       sparFac = sparFac*accFac;
      end
      absLm(absLm<tolr)=tolr*sparFac;


      % build the diagonal weighting matrix for this iteration
      R=diag(1./absLm);

      mold=m;

    %   % use LSQR to get m that minimizes argmin||[            M              ]f - [d]||2
    %   %                                           sqrt(alpha/2)*sqrt(absLm)*L      0
    %   display('LSQR...')
    %   tic
    %   A_bottom = sqrt(alpha/2)*sqrt(R)*L;
    %   A = [G;A_bottom];
    %   b = [d;zeros(size(A_bottom,1),1)];
    %   m = lsqr(A,b,tolx,maxiter);
    %   toc
    %   % get the new iterate and check for convergance
    %   display('QR...')
    %   tic
    %   [Q,Rr] = qr(2*GTG+alpha*L'*R*L);
    %   m=Rr\(Q'*(2*GTd));
    %   toc
    %   display(norm(m))
      display(['(iteration:' num2str(iter) ')'])
      tic
      m=(2*GTG+alpha*L'*R*L)\(2*GTd);
      toc

      if (norm(m-mold)/(1+norm(mold)) < tolx) %|| norm(m-mold)<m_diff
        mreg=m;
        display(['norm(m-mold)=' num2str(norm(m-mold)) ', 1+norm(mold)=' num2str(1+norm(mold)) ', norm(m-mold)/(1+norm(mold))=' ...
          num2str(norm(m-mold)/(1+norm(mold)))])
        return
      else
          display(['norm(m-mold)=' num2str(norm(m-mold)) ', 1+norm(mold)=' num2str(1+norm(mold)) ', norm(m-mold)/(1+norm(mold))=' ...
          num2str(norm(m-mold)/(1+norm(mold)))])
      end
    end

    % Give a warning, if desired, but return best solution.
    display('L1 norm regularization: maximum iterations exceeded.')
    display(['Final ratio= norm(m-mold)/(1+norm(mold))='  num2str(norm(m-mold)/(1+norm(mold)))]);
    mreg=m;
end