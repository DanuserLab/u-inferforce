function [h2,uMap]=generateHeatmapFromGridData(x_mat_u,y_mat_u,ux,uy,dataPath,band,umin,umax,quiverTrue,w,h)
imSizeX = x_mat_u(end,end)-x_mat_u(1,1);
imSizeY = y_mat_u(end,end)-y_mat_u(1,1);
if nargin<5
    dataPath=[];
    band=0;
    quiverTrue=true;
    w = imSizeX;
    h = imSizeY;
elseif nargin<9
    quiverTrue=true;
    w = imSizeX;
    h = imSizeY;
elseif nargin<10
    w = imSizeX;
    h = imSizeY;
end
centerX = ((x_mat_u(end,end)+x_mat_u(1,1))/2);
centerY = ((y_mat_u(end,end)+y_mat_u(1,1))/2);
xmin = centerX-w/2+band;
xmax = centerX+w/2-band;
ymin = centerY-h/2+band;
ymax = centerY+h/2-band;
[XI,YI] = meshgrid(xmin:xmax,ymin:ymax);

% [XI,YI]=meshgrid(x_mat_u(1,1):x_mat_u(1,1,1)+imSizeX,y_mat_u(1,1):y_mat_u(1,1)+imSizeY);
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
unorm = (ux.^2 + uy.^2).^0.5;
uMap = griddata(x_mat_u,y_mat_u,unorm,XI,YI,'linear');
if nargin<7 || isempty(umin)
    umin = min(uMap(:));
end
if nargin<8 || isempty(umax)
    umax = max(uMap(:));
end

h2=figure('color','w');
set(h2, 'Position', [100 100 w*1.25 h])
subplot('Position',[0 0 0.8 1])
imshow(uMap,[umin umax]), colormap jet;
%quiver plot
%quiver
% unit vector plot
hold on
cfactor = 2;
grid_mat(:,:,2) = y_mat_u;
grid_mat(:,:,1) = x_mat_u;
xmax = size(y_mat_u,1);%y_mat_u(end,end);
ymax =size(x_mat_u,2);% x_mat_u(end,end);
grid_mat_coarse = grid_mat(cfactor-round(cfactor/2):cfactor:xmax,cfactor-round(cfactor/2):cfactor:ymax,:);
umat_vecx = reshape(ux(cfactor-round(cfactor/2):cfactor:xmax,cfactor-round(cfactor/2):cfactor:ymax),[],1);
umat_vecy = reshape(uy(cfactor-round(cfactor/2):cfactor:xmax,cfactor-round(cfactor/2):cfactor:ymax),[],1);
pos_vecx = reshape(grid_mat_coarse(:,:,1),[],1);
pos_vecy = reshape(grid_mat_coarse(:,:,2),[],1);
dispScale=0.04*umax;%max(sqrt(umat_vecx.^2+umat_vecy.^2));

xmin = centerX-w/2;
xmax = centerX+w/2;
ymin = centerY-h/2;
ymax = centerY+h/2;
Npoints = length(umat_vecx);
inIdx = false(Npoints,1);

for ii=1:Npoints
    if pos_vecx(ii)>xmin+1 && pos_vecx(ii)<xmax-1 ...
            && pos_vecy(ii)>ymin+1 && pos_vecy(ii)<ymax-1
        inIdx(ii) = true;
    end
end

if quiverTrue
    quiver(pos_vecx(inIdx)-xmin+1,pos_vecy(inIdx)-ymin+1,umat_vecx(inIdx)./dispScale,umat_vecy(inIdx)./dispScale,0,'Color',[75/255 0/255 130/255]);
end

% quiver(pos_vecx-grid_mat(1,1,1),pos_vecy-grid_mat(1,1,2), umat_vecx./dispScale,umat_vecy./dispScale,0,'Color',[75/255 0/255 130/255]);

subplot('Position',[0.8 0.1 0.1 0.8])
axis tight
caxis([umin umax]), axis off
hc = colorbar('West');
set(hc,'Fontsize',12)

% saving
% Set up the output file path
if ~isempty(dataPath)
    outputFilePath = [dataPath filesep 'heatMap'];
    tifPath = [outputFilePath filesep 'tifs'];
    figPath = [outputFilePath filesep 'figs'];
    epsPath = [outputFilePath filesep 'eps'];
    if ~exist(tifPath,'dir') || ~exist(epsPath,'dir')
        mkdir(tifPath);
        mkdir(figPath);
        mkdir(epsPath);
    end

    I = getframe(h2);
    imwrite(I.cdata, strcat(tifPath,'/hMapTif','.tif'));
    hgsave(h2,strcat(figPath,'/hMapFig.fig'),'-v7.3')
    print(h2,strcat(epsPath,'/hMapEps.eps'),'-depsc2')
end
end
    