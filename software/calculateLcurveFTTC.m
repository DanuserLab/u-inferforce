function  [rho,eta,reg_corner,alphas] = calculateLcurveFTTC(grid_mat,u,E,s, cluster_size, i_max, j_max, L, LcurveFactor)
% this function calculateLcurveFTTC calculates L-curve by obtaining residual norm and self-norm
% Input : grid_mat, u, cluster_size have to be in the same units, namely
%         pixels. If they are given in different units e.g. meters, then a
%         different scaling factor for the elastic energy has to be
%         applied! The force value remains, however, unaffected, see below.
% Output: The output force is actually a surface stress with the same units
%         as the input E! In particular, the unit of the output force is
%         independent of the units of the input grid_mat,u and cluster_size
%         The reason for this is essentially that the elastic stress is
%         only dependent on the non-dimensional strain which is given by
%         spatial derivatives of the displacements, that is du/dx. If u and
%         dx (essentially cluster_size) are in the same units, then the
%         resulting force has the same dimension as the input E.
%     LcurveFactor = 10;
% 
% Instead of relying on L-curve, we'll use L-curve between alphas vs eta, which is more
% appropriate for FTTC. Stricker's method doesn't makes sense either.
% Sangyoon Han, Oct 2017
    %% set up
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
    alphas=10.^(log10(L)-5:1.25/LcurveFactor:log10(L)+2.5);
    rho=zeros(length(alphas),1);
    eta=zeros(length(alphas),1);
%     etaMax=zeros(length(alphas),1);
    grid_mat1=grid_mat(:,:,1);
    grid_mat2=grid_mat(:,:,2);
    u1=u(:,:,1);
    u2=u(:,:,2);
    firstGridPointX=grid_mat(1,1,1);
    endGridPointX=grid_mat(end,end,1);
    firstGridPointY=grid_mat(1,1,2);
    endGridPointY=grid_mat(end,end,2);
    %% loop
    parfor i=1:length(alphas)
        disp(['testing L = ' num2str(alphas(i)) '... '])
        curL=alphas(i);
        V = 2*(1+s)/E;

        kx_vec = 2*pi/i_max/cluster_size.*[0:(i_max/2-1) (-i_max/2:-1)];
        ky_vec = 2*pi/j_max/cluster_size.*[0:(j_max/2-1) (-j_max/2:-1)];
        kx = repmat(kx_vec',1,j_max);
        ky = repmat(ky_vec,i_max,1);
        Rx=ones(size(kx));
        Ry=ones(size(ky));

        kx(1,1) = 1;
        ky(1,1) = 1;

        X = i_max*cluster_size/2;
        Y = j_max*cluster_size/2; 

        g0x = pi.^(-1).*V.*((-1).*Y.*log((-1).*X+sqrt(X.^2+Y.^2))+Y.*log( ...
          X+sqrt(X.^2+Y.^2))+((-1)+s).*X.*(log((-1).*Y+sqrt(X.^2+Y.^2) ...
          )+(-1).*log(Y+sqrt(X.^2+Y.^2))));
        g0y = pi.^(-1).*V.*(((-1)+s).*Y.*(log((-1).*X+sqrt(X.^2+Y.^2))+( ...
          -1).*log(X+sqrt(X.^2+Y.^2)))+X.*((-1).*log((-1).*Y+sqrt( ...
          X.^2+Y.^2))+log(Y+sqrt(X.^2+Y.^2))));
        % kx, ky: wave vector
        k = (kx.^2+ky.^2).^(-1/2);
        Ginv_xx =k.*V.*(kx.^2.*curL.*Rx+ky.^2.*curL.*Ry+V.^2).^(-1).*(kx.^2.* ...
                  curL.*Rx+ky.^2.*curL.*Ry+((-1)+s).^2.*V.^2).^(-1).*(kx.^4.*(curL.*Rx+(-1).*curL.*s.*Rx)+ ...
                  kx.^2.*((-1).*ky.^2.*curL.*Ry.*((-2)+s)+(-1).*((-1)+s).*V.^2)+ky.^2.*( ...
                  ky.^2.*curL.*Ry+((-1)+s).^2.*V.^2));
        Ginv_yy = k.*V.*(kx.^2.*curL+ky.^2.*curL.*Ry+V.^2).^(-1).*(kx.^2.* ...
                  curL.*Rx+ky.^2.*curL.*Ry+((-1)+s).^2.*V.^2).^(-1).*(kx.^4.*curL+(-1).*ky.^2.*((-1)+ ...
                  s).*(ky.^2.*curL.*Rx+V.^2)+kx.^2.*((-1).*ky.^2.*curL.*Ry.*((-2)+s)+((-1)+s).^2.* ...
                  V.^2));
        Ginv_xy = (-1).*kx.*ky.*k.*s.*V.*(kx.^2.*curL.*Rx+ky.^2.*curL.*Ry+ ...
                  V.^2).^(-1).*(kx.^2.*curL.*Rx+ky.^2.*curL.*Ry+((-1)+s).*V.^2).*(kx.^2.*curL.*Rx+ky.^2.* ...
                  curL.*Ry+((-1)+s).^2.*V.^2).^(-1);


        Ginv_xx(1,1) = 1/g0x;
        Ginv_yy(1,1) = 1/g0y;
        Ginv_xy(1,1) = 0;

        Ginv_xy(i_max/2+1,:) = 0;
        Ginv_xy(:,j_max/2+1) = 0;

        Ftu1 = fft2(u1);
        Ftu2 = fft2(u2);

        Ftf1 = Ginv_xx.*Ftu1 + Ginv_xy.*Ftu2;
        Ftf2 = Ginv_xy.*Ftu1 + Ginv_yy.*Ftu2;

        f1 = ifft2(Ftf1,'symmetric');
        f2 = ifft2(Ftf2,'symmetric');

        force1 = reshape(f1,i_max*j_max,1);
        force2 = reshape(f2,i_max*j_max,1);     

        fnorm = (force1.^2 + force2.^2).^0.5;
        eta(i) = sum(fnorm(:));
%         etaMax(i) = max(fnorm(:));
        % residual norm calculation
%         G_xx=V*((1-s).*k.^2+s*ky.^2)./(k.^3);
%         G_xy=V*s*kx.*ky./(k.^3);
%         G_yy=V*((1-s).*k.^2+s*kx.^2)./(k.^3);
% 
% %         Ginv_xx =k*(-kx.^2+ky.^2.*(-1+s))/(V*(-1+s));
% % 
% %         Ginv_xx =k.*(kx.^2+(1-s)*ky.^2)./(V*(1-s));
% %         G_xxfor =V.*(1-s)./(k.*(k.^2-ky.^2.*s));
%         
% %         G_xx(1,1) = g0x;
% %         G_yy(1,1) = g0y;
%         G_xx(1,1) = 0;
%         G_yy(1,1) = 0;
%         G_xy(1,1) = 0;
%         G_xy(i_max/2+1,:) = 0;
%         G_xy(:,j_max/2+1) = 0;
% 
%         Ftu_estm(:,:,1) = G_xx.*Ftf(:,:,1)+G_xy.*Ftf(:,:,2);
%         Ftu_estm(:,:,2) = G_xy.*Ftf(:,:,1)+G_yy.*Ftf(:,:,2);
% %         ue(:,:,1) = ifft2(Ftu_estm(:,:,1),'symmetric');
% %         ue(:,:,2) = ifft2(Ftu_estm(:,:,2),'symmetric');
%         ue(:,:,1) = real(ifft2(Ftu_estm(:,:,1),'symmetric'));
%         ue(:,:,2) = real(ifft2(Ftu_estm(:,:,2),'symmetric'));
%         
%         figure, quiver(grid_mat(:,:,1),grid_mat(:,:,2),u(:,:,1),u(:,:,2),0)     
%         hold on
%         quiver(grid_mat(:,:,1),grid_mat(:,:,2),ue(:,:,1),ue(:,:,2),0,'r')
        
        [ue1,ue2]=fwdSolution(grid_mat1,grid_mat2,E,...
            firstGridPointX,endGridPointX,firstGridPointY,endGridPointY,...
            f1,f2,'fft','Intp',2^10,[],0.5,false,false);
%         [ue(:,:,1),ue(:,:,2)]=fwdSolution(grid_mat(:,:,1),grid_mat(:,:,2),E,...
%             grid_mat(1,1,1),grid_mat(end,end,1),grid_mat(1,1,2),grid_mat(end,end,2),...
%             f(:,:,1),f(:,:,2),'fft','noIntp',[],[],0.5,false,true);
        residueU = ((ue1-u1).^2+(ue2-u2).^2).^0.5;
        rho(i) = sum(residueU(:));
    end
    %% Find the L-corner
%     save(LcurveDataPath,'rho','eta','alphas','L','msparse','-v7.3'); % saving before selection.
%     [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(rho,eta,alphas,L);
%     [reg_corner,~,~]=regParamSelecetionLcurve(rho,eta,alphas,L);
    [reg_corner,ireg_corner,~]=regParamSelecetionLcurve(alphas',eta,alphas,L);
%     figure, plot(alphas,etaMax,'-'); hold on, plot(alphas(ireg_corner),etaMax(ireg_corner),'ro')
end