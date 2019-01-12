function [ux,uy,x_grid,y_grid,meshPtsFwdSol]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,method,opt,meshPtsFwdSol,h,v,refine,useSameSampling)
% This forward solution is only valid for a Poisson's ratio v=0.5 if not
% specified.
% Input: No matter what the dimension of x0 and y0 is (pix, or um), the
%        dimension of the surface stresses (force_x, force_y) must have the
%        same dimension as the elastic modulus E, usually Pa.
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

% Output: The calculated ux and uy have the same dimension as the input
%         x0, y0.
% Achim Besser 2012
% Updated with allowing Poisson's ratio other than 0.5.
% Sangyoon Han, August 2014
if nargin <14
    v=0.5;
    refine = true;
    useSameSampling = false;
elseif nargin <15
    refine = true;
    useSameSampling = false;
elseif nargin <16
    useSameSampling = false;
end
if strcmpi(method,'conv_free')
    tic;
    display('Calulate the convolution explicitely in free triangulated mesh')
    [nRow,~]=size(x0);

    ux = zeros(nRow,1);
    uy = zeros(nRow,1);
    for i=1:nRow
        integrandx = @(x,y) boussinesqGreens(1,1,x0(i)-x,y0(i)-y,E,v).*force_x(x,y) + boussinesqGreens(1,2,x0(i)-x,y0(i)-y,E,v).*force_y(x,y);
        integrandy = @(x,y) boussinesqGreens(2,1,x0(i)-x,y0(i)-y,E,v).*force_x(x,y) + boussinesqGreens(2,2,x0(i)-x,y0(i)-y,E,v).*force_y(x,y);

        ux(i) = quad2d(integrandx,xmin,xmax,ymin,ymax,'MaxFunEvals',10^5,'AbsTol',5e-6);% RelTol sucks! 'RelTol',5e-13);
        uy(i) = quad2d(integrandy,xmin,xmax,ymin,ymax,'MaxFunEvals',10^5,'AbsTol',5e-6);% RelTol sucks! 'RelTol',5e-13);
    end
    toc;
elseif nargin<10 || strcmpi(method,'conv')
    tic;
    display('Calulate the convolution explicitely')
    [nRow,nCol]=size(x0);

    for i=1:nRow
        for j=1:nCol  
            integrandx = @(x,y) boussinesqGreens(1,1,x0(i,j)-x,y0(i,j)-y,E,v).*force_x(x,y) + boussinesqGreens(1,2,x0(i,j)-x,y0(i,j)-y,E,v).*force_y(x,y);
            integrandy = @(x,y) boussinesqGreens(2,1,x0(i,j)-x,y0(i,j)-y,E,v).*force_x(x,y) + boussinesqGreens(2,2,x0(i,j)-x,y0(i,j)-y,E,v).*force_y(x,y);

            ux(i,j) = quad2d(integrandx,xmin,xmax,ymin,ymax,'MaxFunEvals',10^10,'AbsTol',5e-10);% RelTol sucks! 'RelTol',5e-13);
            uy(i,j) = quad2d(integrandy,xmin,xmax,ymin,ymax,'MaxFunEvals',10^10,'AbsTol',5e-10);% RelTol sucks! 'RelTol',5e-13);
        end
    end
    toc;
elseif strcmpi(method,'fft')
    %display('Use fast convolution')    

    %***************************************************************
    % Here starts the calculation using the fast fourier transform *
    %***************************************************************
    
    % Number of points to calculate the force field and the Greensfunction.
    % Since below Nx_G and Ny_G are odd, the Greensfunctions will be
    % evaluated at zero. Since the Greensfunctions at x=0 diverges, this
    % will in general cause a problem. For this I have found a work around.
    % The Greensfunction will be evaluated as is and the divergent value
    % at x=0 will be set to zero, this part of the support will be
    % integrated seperately. In order to do this, I assume that for dense
    % sampling, the force field doesn't very strongly over one gridsize.
    % Assuming it to be constant around x=0, allows to integrate the
    % Greensfunction around a domain x=+-r and y=+-r. This yields a
    % correction term which is of particular importance for sparse
    % sampling, meaning that Nx_F is small. This alogorithm performs very
    % well and has been cross-validated with the results obtained using the
    % 'conv' option. If you want to repeat the test, use the few lines of
    % code at the very end of this function.
    
    %tic;
    
    % This determines the sampling of the force field:
    if (nargin < 12 || isempty(meshPtsFwdSol)) && ~useSameSampling
        display('Use meshPtsFwdSol=2^10. This value should be given with the function call!!!');
        meshPtsFwdSol=2^10;
    end
    
    if useSameSampling
        Nx_F=size(x0,2); % 2^10 is the densest sampling possible.
        Ny_F=size(y0,1);
    else
        Nx_F=meshPtsFwdSol; % 2^10 is the densest sampling possible.
        Ny_F=Nx_F;
    end
    
    % To account for dx*dy in the convolution integral one has to finally
    % rescale the result by the following scaling factor:
    xRange=(max(max(x0))-min(min(x0)));
    yRange=(max(max(y0))-min(min(y0)));
    scalingFactor=(xRange*yRange)/(Nx_F*Ny_F);
%     scalingFactor=1;
    
    % To cover the whole support of the force field, the domain over which
    % the Greensfunctions have to be calculated need to be at least of size:
    % (2Nx-1)x(2Ny-1).

    Nx_G=2*Nx_F-1;
    Ny_G=2*Ny_F-1;    
    
    % Subsequently, these have to be padded with zeros according to:
    Nx_pad=Nx_F+Nx_G-1;
    Ny_pad=Ny_F+Ny_G-1;
    
    % These might not be a power of 2, make sure that they are:
    Nx_pad=pow2(nextpow2(Nx_pad));
    Ny_pad=pow2(nextpow2(Ny_pad));

    % First determine the boundaries of the mesh:
    leftUpperCorner =[min(min(x0)) min(min(y0))];
    rightLowerCorner=[max(max(x0)) max(max(y0))];

    % create a regular mesh with Nx*Ny meshpoints where the force field is
    % calculated. This need not to be a power of 2 yet:
    xvec_F=linspace(leftUpperCorner(1),rightLowerCorner(1),Nx_F);
    yvec_F=linspace(leftUpperCorner(2),rightLowerCorner(2),Ny_F);
    [xgrid_F,ygrid_F]=meshgrid(xvec_F,yvec_F);
    
    % create a mesh centered at zero with Nx_G*Ny_G meshpoints, where the
    % Greensfunctions are calculated.
    xvec_G=linspace(-xRange,xRange,Nx_G);
    yvec_G=linspace(-yRange,yRange,Ny_G);
    
%     [xgrid_G,ygrid_G]=meshgrid(yvec_G,xvec_G);
    [xgrid_G,ygrid_G]=meshgrid(xvec_G,yvec_G);
      
    %calculate the force values at the grid_F positions:
    if useSameSampling
        discrete_Force_x_unPadded=force_x; %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
        discrete_Force_y_unPadded=force_y;
    else
    % Make force_x and force_y a function handle if it is a matrix
        if ismatrix(force_x) && ismatrix(force_y) && ~isa(force_x,'TriScatteredInterp') && ~isa(force_x,'function_handle')
        %     [xmat,ymat]=meshgrid(xmin:xmax,ymin:ymax);
        %     xvec=xmat(:);
        %     yvec=ymat(:);
            xvec=x0(:);
            yvec=y0(:);
            force_x_vec=force_x(:);
            force_y_vec=force_y(:);
            force_x = scatteredInterpolant(xvec,yvec,force_x_vec);
            force_y = scatteredInterpolant(xvec,yvec,force_y_vec);
        end
        discrete_Force_x_unPadded=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
        discrete_Force_y_unPadded=force_y(xgrid_F,ygrid_F);
        % taking care of nans
        checkVec=isnan(discrete_Force_x_unPadded);
        discrete_Force_x_unPadded(checkVec)=0;
        checkVec=isnan(discrete_Force_y_unPadded);
        discrete_Force_y_unPadded(checkVec)=0;
    end
    % Calculate the Greens-function values at the grid_G positions. This can
    % be improved since the Greensfunction never change for a given grid
    % size. When the Basis functions are calculated this has to be done
    % only once (as well as the FFT for these fields!!!):
    discrete_boussinesqGreens11=boussinesqGreens(1,1,xgrid_G,ygrid_G,E,v);
    discrete_boussinesqGreens12=boussinesqGreens(1,2,xgrid_G,ygrid_G,E,v);
   %discrete_boussinesqGreens21=discrete_boussinesqGreens12;
    discrete_boussinesqGreens22=boussinesqGreens(2,2,xgrid_G,ygrid_G,E,v);
    
    % Pad the calculated fields with zero to the next power larger than 
    % (2*N-1), see above. For this setup, the FFT is fastest.
%     discrete_Force_x=padarray(discrete_Force_x_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');%'symmetric','post');
%     discrete_Force_y=padarray(discrete_Force_y_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
    discrete_Force_x=padarray(discrete_Force_x_unPadded,[Ny_pad-Ny_F Nx_pad-Nx_F],0,'post');%'symmetric','post');
    discrete_Force_y=padarray(discrete_Force_y_unPadded,[Ny_pad-Ny_F Nx_pad-Nx_F],0,'post');
    
%     discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%     discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
%    %discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post'); 
%     discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    
    discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Ny_pad-Ny_G Nx_pad-Nx_G],0,'post');
    discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Ny_pad-Ny_G Nx_pad-Nx_G],0,'post');
    discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Ny_pad-Ny_G Nx_pad-Nx_G],0,'post');
    
    % Now calculate the fourier transforms:
    dFT_Force_x=fft2(discrete_Force_x);
    clear discrete_Force_x;
    dFT_Force_y=fft2(discrete_Force_y);
    clear discrete_Force_y;
    
    % This has to be calculated only once for all basis functions!
    dFT_boussinesqGreens11=fft2(discrete_boussinesqGreens11);
    clear discrete_boussinesqGreens11;
    dFT_boussinesqGreens12=fft2(discrete_boussinesqGreens12);
    clear discrete_boussinesqGreens12;
    dFT_boussinesqGreens21=dFT_boussinesqGreens12;
    % nothing to clear here!
    dFT_boussinesqGreens22=fft2(discrete_boussinesqGreens22);
    clear discrete_boussinesqGreens22;
    
    % Now calculate the solution:                
    ux_grid=ifft2(dFT_boussinesqGreens11.*dFT_Force_x+dFT_boussinesqGreens12.*dFT_Force_y);
    clear dFT_boussinesqGreens11 dFT_boussinesqGreens12;
    uy_grid=ifft2(dFT_boussinesqGreens21.*dFT_Force_x+dFT_boussinesqGreens22.*dFT_Force_y);
    clear dFT_boussinesqGreens21 dFT_Force_x dFT_boussinesqGreens22 dFT_Force_y;
    
    % Plot the solution:
%     figure(10)
%     imshow(ux_grid,[])
%     
%     figure(11)
%     surf(uy_grid)
    
    % Now extract the essential part from the solution. It is located in
    % the center of the padded field.    
    % I really don't understand why to cut it out like this, but it works!
    startIndex_x=abs(Nx_G-Nx_F)+1; % Or is it just: startIndex_x=Nx_F
    startIndex_y=abs(Ny_G-Ny_F)+1;
    
    endIndex_x=startIndex_x+Nx_F-1;
    endIndex_y=startIndex_y+Ny_F-1;
    

    % Remove imaginary part caused by round off errors:
%     ux_grid=real(ux_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
%     uy_grid=real(uy_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
    ux_grid=real(ux_grid(startIndex_y:endIndex_y,startIndex_x:endIndex_x));
    uy_grid=real(uy_grid(startIndex_y:endIndex_y,startIndex_x:endIndex_x));

%!!! This could be improved by using the analytical solution for the Fourie
%!!! Transform of the Greensfunction!
    % Add the solution for G(0,0). This is a correction term which becomes
    % irrelevant for very dense sampling. But for small Nx_F it is REALLY
    % essential!
    % Set the Poisson's ratio to 0.5:
    v=0.5;
    dx=abs(xvec_G(2)-xvec_G(1));
    dy=abs(yvec_G(2)-yvec_G(1));
    
    int_x2_over_r3=2*dy*log((dy^2+2*dx*(dx+sqrt(dx^2+dy^2)))/(dy^2));    
    int_y2_over_r3=2*dx*log((dx^2+2*dy*(dy+sqrt(dx^2+dy^2)))/(dx^2));    
    int_1_over_r  =int_x2_over_r3 + int_y2_over_r3;
        
    corrTerm_11=(1+v)/(pi*E)*((1-v)*int_1_over_r+v*int_x2_over_r3);
    corrTerm_22=(1+v)/(pi*E)*((1-v)*int_1_over_r+v*int_y2_over_r3);
    
    ux_grid=ux_grid+discrete_Force_x_unPadded*corrTerm_11;
    clear discrete_Force_x_unPadded;
    uy_grid=uy_grid+discrete_Force_y_unPadded*corrTerm_22;
    clear discrete_Force_y_unPadded;
    
    % scale the solution appropriately!
    ux_grid=scalingFactor*ux_grid;
    uy_grid=scalingFactor*uy_grid;
    
    
    
    % Recursive call to fwdSolution:   
%     refine =true;
    if ~isempty(xmin) && ~isempty(xmax) && ~isempty(ymin) && ~isempty(ymax) && refine
        % and the support of the force is much small than ROI for the
        % displacement, this check should be included otherwise one calculate
        % the same thing twice
        % extend the size of the region a little bit, here by a factor of 1.2,
        % but this is arbitrary.
        % The range of the support is:
        xsupp=xmax-xmin;
        ysupp=ymax-ymin;
        
        % Now expand the x- and y-range:
        xminExp=xmin-xsupp/10;
        xmaxExp=xmax+xsupp/10;
        yminExp=ymin-ysupp/10;
        ymaxExp=ymax+ysupp/10;
        
        [ux_fine uy_fine x_grid_fine y_grid_fine]=fwdSolution([xminExp xmaxExp],[yminExp ymaxExp],E,[],[],[],[],force_x,force_y,method,'noIntp',meshPtsFwdSol);
        
        % Later on we want to have ux and uy defined on a regular
        % grid. for this reason we now interpolate the fine solution onto
        % the regular xgrid_F and ygrid_F. This grid is usually so fine
        % that it already corresponds to subpixel sampling: (e.g. 2^10=1024
        % positions along each image dimension)
        
        % find the positions that are within
        
        iux_fine = griddata(x_grid_fine,y_grid_fine,ux_fine,xgrid_F,ygrid_F,'linear');%'*cubic'
        iuy_fine = griddata(x_grid_fine,y_grid_fine,uy_fine,xgrid_F,ygrid_F,'linear');%'*linear'
        
        % those contain a lot of NaNs since most of the points are outside
        % of [xminExp xmaxExp] and [yminExp ymaxExp]. This will give us two
        % masks that we can now use to update ux_grid and uy_grid:        
        valMat=~isnan(iux_fine);
        
        diff=2*abs(uy_grid(valMat)-iuy_fine(valMat))./abs(uy_grid(valMat)+iuy_fine(valMat));
        display(['max. correction: ',num2str(max(diff(:)))]);
        
        % Update the values in the coarse solution:
        ux_grid(valMat)=iux_fine(valMat);
        uy_grid(valMat)=iuy_fine(valMat);
    end
    
        
    % Now interpolate the displacement field from the regular grid to the irregular
    % measured grid. Since we have the force field defined on a regular grid
    % we can use the fast *option. 'linear' is about a factor of two faster than
    % 'cubic'. Hard to tell if cubic performs better than linear.    
    if nargin>10 && strcmp(opt,'noIntp')
        ux = ux_grid;%'*cubic'
        uy = uy_grid;%'*linear'
        x_grid=xgrid_F;
        y_grid=ygrid_F;
    else
        ux = interp2(xgrid_F,ygrid_F,ux_grid,x0,y0,'*cubic');%'*cubic'
        uy = interp2(xgrid_F,ygrid_F,uy_grid,x0,y0,'*cubic');%'*linear'
        x_grid=x0;
        y_grid=y0;
    end
    
    %toc;
    
    %x0_vec=reshape(x0,[],1);
    %y0_vec=reshape(y0,[],1);

    %x_vec=reshape(xgrid_F,[],1);
    %y_vec=reshape(ygrid_F,[],1);
    %ux_vec=reshape(ux_grid,[],1);
    %uy_vec=reshape(uy_grid,[],1);


    %tic;
    %[~, ~, ux] = griddata(x_vec,y_vec,ux_vec,x0,y0,'cubic');
    %[~, ~, uy] = griddata(x_vec,y_vec,uy_vec,x0,y0,'cubic');
    %toc;
elseif strcmpi(method,'fft_finite')
    %display('Use fast convolution')    

    % This determines the sampling of the force field:
    if nargin < 12 || isempty(meshPtsFwdSol)
        display('Use meshPtsFwdSol=2^10. This value should be given with the function call!!!');
        meshPtsFwdSol=2^10;
    end
        
    Nx_F=meshPtsFwdSol; % 2^10 is the densest sampling possible.
    Ny_F=Nx_F;
    
    % To account for dx*dy in the convolution integral one has to finally
    % rescale the result by the following scaling factor:
    xRange=(max(max(x0))-min(min(x0)));
    yRange=(max(max(y0))-min(min(y0)));
    scalingFactor=(xRange*yRange)/(Nx_F*Ny_F);
    
    % To cover the whole support of the force field, the domain over which
    % the Greensfunctions have to be calculated need to be at least of size:
    % (2Nx-1)x(2Ny-1).

    Nx_G=2*Nx_F-1;
    Ny_G=2*Ny_F-1;    
    
    % Subsequently, these have to be padded with zeros according to:
    Nx_pad=Nx_F+Nx_G-1;
    Ny_pad=Ny_F+Ny_G-1;
    
    % These might not be a power of 2, make sure that they are:
    Nx_pad=pow2(nextpow2(Nx_pad));
    Ny_pad=pow2(nextpow2(Ny_pad));

    % First determine the boundaries of the mesh:
    leftUpperCorner =[min(min(x0)) min(min(y0))];
    rightLowerCorner=[max(max(x0)) max(max(y0))];

    % create a regular mesh with Nx*Ny meshpoints where the force field is
    % calculated. This need not to be a power of 2 yet:
    xvec_F=linspace(leftUpperCorner(1),rightLowerCorner(1),Nx_F);
    yvec_F=linspace(leftUpperCorner(2),rightLowerCorner(2),Ny_F);
    [xgrid_F,ygrid_F]=meshgrid(xvec_F,yvec_F);
    
    % create a mesh centered at zero with Nx_G*Ny_G meshpoints, where the
    % Greensfunctions are calculated.
    xvec_G=linspace(-xRange,xRange,Nx_G);
    yvec_G=linspace(-yRange,yRange,Ny_G);
    
    [xgrid_G,ygrid_G]=meshgrid(xvec_G,yvec_G);
      
    %calculate the force values at the grid_F positions:
    discrete_Force_x_unPadded=force_x(xgrid_F,ygrid_F); %this has only to be calculated over the support xmin,xmax,ymin,ymax rest is zero
    discrete_Force_y_unPadded=force_y(xgrid_F,ygrid_F);

    % This part if what's different from 'fft' method that uses
    % boussinesque greens function
    discrete_boussinesqGreens11=finiteThicknessGreens(1,1,xgrid_G,ygrid_G,E,h);
    discrete_boussinesqGreens12=finiteThicknessGreens(1,2,xgrid_G,ygrid_G,E,h);
   %discrete_boussinesqGreens21=discrete_boussinesqGreens12;
    discrete_boussinesqGreens22=finiteThicknessGreens(2,2,xgrid_G,ygrid_G,E,h);
    
    % Pad the calculated fields with zero to the next power larger than 
    % (2*N-1), see above. For this setup, the FFT is fastest.
    discrete_Force_x=padarray(discrete_Force_x_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');%'symmetric','post');
    discrete_Force_y=padarray(discrete_Force_y_unPadded,[Nx_pad-Nx_F Ny_pad-Ny_F],0,'post');
    
    discrete_boussinesqGreens11=padarray(discrete_boussinesqGreens11,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    discrete_boussinesqGreens12=padarray(discrete_boussinesqGreens12,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
   %discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post'); 
    discrete_boussinesqGreens22=padarray(discrete_boussinesqGreens22,[Nx_pad-Nx_G Ny_pad-Ny_G],0,'post');
    
    % Now calculate the fourier transforms:
    dFT_Force_x=fft2(discrete_Force_x);
    clear discrete_Force_x;
    dFT_Force_y=fft2(discrete_Force_y);
    clear discrete_Force_y;
    
    % This has to be calculated only once for all basis functions!
    dFT_boussinesqGreens11=fft2(discrete_boussinesqGreens11);
    clear discrete_boussinesqGreens11;
    dFT_boussinesqGreens12=fft2(discrete_boussinesqGreens12);
    clear discrete_boussinesqGreens12;
    dFT_boussinesqGreens21=dFT_boussinesqGreens12;
    % nothing to clear here!
    dFT_boussinesqGreens22=fft2(discrete_boussinesqGreens22);
    clear discrete_boussinesqGreens22;
    
    % Now calculate the solution:                
    ux_grid=ifft2(dFT_boussinesqGreens11.*dFT_Force_x+dFT_boussinesqGreens12.*dFT_Force_y);
    clear dFT_boussinesqGreens11 dFT_boussinesqGreens12;
    uy_grid=ifft2(dFT_boussinesqGreens21.*dFT_Force_x+dFT_boussinesqGreens22.*dFT_Force_y);
    clear dFT_boussinesqGreens21 dFT_Force_x dFT_boussinesqGreens22 dFT_Force_y;
    
    % Now extract the essential part from the solution. It is located in
    % the center of the padded field.    
    % I really don't understand why to cut it out like this, but it works!
    startIndex_x=abs(Nx_G-Nx_F)+1; % Or is it just: startIndex_x=Nx_F
    startIndex_y=abs(Ny_G-Ny_F)+1;
    
    endIndex_x=startIndex_x+Nx_F-1;
    endIndex_y=startIndex_y+Ny_F-1;
    

    % Remove imaginary part caused by round off errors:
    ux_grid=real(ux_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));
    uy_grid=real(uy_grid(startIndex_x:endIndex_x,startIndex_y:endIndex_y));

%     figure
%     imshow(ux_grid,[])
    % scale the solution appropriately!
    ux_grid=scalingFactor*ux_grid;
    uy_grid=scalingFactor*uy_grid;
    
    % Now interpolate the displacement field from the regular grid to the irregular
    % measured grid. Since we have the force field defined on a regular grid
    % we can use the fast *option. 'linear' is about a factor of two faster than
    % 'cubic'. Hard to tell if cubic performs better than linear.    
    if nargin>10 && strcmp(opt,'noIntp')
        ux = ux_grid;%'*cubic'
        uy = uy_grid;%'*linear'
        x_grid=xgrid_F;
        y_grid=ygrid_F;
    else
        ux = interp2(xgrid_F,ygrid_F,ux_grid,x0,y0,'*cubic');%'*cubic'
        uy = interp2(xgrid_F,ygrid_F,uy_grid,x0,y0,'*cubic');%'*linear'
        x_grid=x0;
        y_grid=y0;
    end
    
    %toc;
    
    %x0_vec=reshape(x0,[],1);
    %y0_vec=reshape(y0,[],1);

    %x_vec=reshape(xgrid_F,[],1);
    %y_vec=reshape(ygrid_F,[],1);
    %ux_vec=reshape(ux_grid,[],1);
    %uy_vec=reshape(uy_grid,[],1);


    %tic;
    %[~, ~, ux] = griddata(x_vec,y_vec,ux_vec,x0,y0,'cubic');
    %[~, ~, uy] = griddata(x_vec,y_vec,uy_vec,x0,y0,'cubic');
    %toc;
end

return;

% to test the example:

% strange!!!:
%x0_vec=linspace(-5,20,30);
%y0_vec=linspace(1,6,15);

x0_vec=linspace(-50,200,51);
y0_vec=linspace(0,600,61);

[x0 y0]=meshgrid(x0_vec,y0_vec);

E=10000;
meshPtsFwdSol=2^10;

xmin=min(x0_vec);
xmax=max(x0_vec);
ymin=min(y0_vec);
ymax=max(y0_vec);

xmin=0;
xmax=2;
ymin=1;
ymax=3;

force_x=@(x,y) assumedForce(1,x,y);
force_y=@(x,y) assumedForce(2,x,y);

[ux_conv uy_conv]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'conv');

figure(1)
quiver(x0,y0,ux_conv,uy_conv);


[ux_fft uy_fft]=fwdSolution(x0,y0,E,xmin,xmax,ymin,ymax,force_x,force_y,'fft',[],meshPtsFwdSol);

% figure(11)
% quiver(x0,y0,ux_fft,uy_fft,'r');

%compare the two results:
scalePlot=0.3*sqrt(max(max(ux_conv.^2+uy_conv.^2)));
figure(40)
quiver(x0,y0,ux_fft/scalePlot,uy_fft/scalePlot,0,'r');
hold on
quiver(x0,y0,ux_conv/scalePlot,uy_conv/scalePlot,0,'g');
hold off

%figure(5)
%quiver(x0,y0,assumedForce(1,x0,y0),assumedForce(2,x0,y0),0,'g');

corr_x=2*abs(ux_fft-ux_conv)./abs(ux_fft+ux_conv);
corr_y=2*abs(uy_fft-uy_conv)./abs(uy_fft+uy_conv);

figure(3)
surf(x0,y0,corr_x)

figure(4)
surf(x0,y0,corr_y)

display(['mean rel. deviation in %: ',num2str(mean(corr_y(:)))]);
display([' max rel. deviation in %: ',num2str(max(corr_y(:)))]);
%uncorr_x=2*abs(ux_fft-ux_conv)./abs(ux_fft+ux_conv)
%uncorr_y=2*abs(uy_fft-uy_conv)./abs(uy_fft+uy_conv)

% uncorr_y-corr_y

display('This should be 0')



% in case of:
% x0_vec=linspace(-10,10,25);
% y0_vec=linspace(-2,2,15);
% Important note: for a sampling of 2^8 in the fast solution,
% mean(corr_y(:)) decreases for increasing numerical precision in the
% direct integration of the fwd solution (up to 'AbsTol' 10^-(8)). Thus, it seems 
% like as if 2^8 sampling in the Fourier solution is almost as precise as what
% can be achieved by direct numerical intergeation (up to 'AbsTol' 10^-(10))... which is amazing!!!
% Fourier sampling: 2^(8)
% AbsTol:       mean(corr_y(:))
% 10^(-10)      0.0067
% 10^(- 9)      0.0070
% 10^(- 8)      0.0074  (here precision: fast sol ~= direct conv)
% 10^(- 7)      0.0210
% 10^(- 6)      0.0442
% 10^(- 5)      0.1079
% 10^(- 4)      0.1179

% Fourier sampling: 2^(10)
% AbsTol:       mean(corr_y(:))
% 10^(-10)      0.0011
% 10^(- 9)      0.0013  (here precision: fast sol ~= direct conv)
% 10^(- 8)      0.0023
% 10^(- 7)      0.0162
% 10^(- 6)      0.0391
% 10^(- 5)      0.1024
% 10^(- 4)      0.1123

% Note also that the best precision is where the displacement is predicted
% to be highest, that is at the force center! The mean rel. difference
% between sampling 2^(8) and 2^(10) is 0.0056 which roughly indicates the
% accuracy of the 2^(8) sampling, which is inline with the results above.

















