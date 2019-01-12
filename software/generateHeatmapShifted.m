function [tMap, tmax, tmin, cropInfo,tMapX,tMapY,reg_grid1] = generateHeatmapShifted(forceField,displField,band)
%[tMap, tmax, tmin, cropInfo] = generateHeatmapShifted(forceField,displField,band)
% generates an image of traction in the place of deformed position defined
% by displField. 
% input: 
%           forceField: traction field with pos and vec
%           displField: displacement field with pos and vec
%           band: pixel band that you want to exclude from the map from the
%           edge
% output:
%           tMap: image of traction magnitude contained in cell array
%           tmax: max value of traction magnitude
%           tmin: min value of traction magnitude
%           cropInfo: pos min and max that is used in creating tMap [xmin,ymin,xmax,ymax]
% Sangyoon Han, Nov, 2014
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

%% tmax and tmin determination
tmax = -1;
tmin = 1e10;
reg_grid1=createRegGridFromDisplField(forceField,1,0); %2=2 times fine interpolation
for ii=1:numel(forceField)
    %Load the saved body force map.
    [~,fmat, ~, ~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid1); %1:cluster size
    fnorm = (fmat(:,:,1).^2 + fmat(:,:,2).^2).^0.5;
    % Boundary cutting - I'll take care of this boundary effect later
    fnorm(end-round(band/2):end,:)=[];
    fnorm(:,end-round(band/2):end)=[];
    fnorm(1:1+round(band/2),:)=[];
    fnorm(:,1:1+round(band/2))=[];
    fnorm_vec = reshape(fnorm,[],1); 

    tmax = max(tmax,max(fnorm_vec));
    tmin = min(tmin,min(fnorm_vec));
end
tmax = 0.8*tmax;
% display(['Estimated force maximum = ' num2str(tmax) ' Pa.'])
%% tMap creation    
% account for if displField contains more than one frame
imSizeX = reg_grid1(end,end,1)-reg_grid1(1,1,1)+1;
imSizeY = reg_grid1(end,end,2)-reg_grid1(1,1,2)+1;
w = imSizeX;
h = imSizeY;
centerX = ((reg_grid1(end,end,1)+reg_grid1(1,1,1))/2);
centerY = ((reg_grid1(end,end,2)+reg_grid1(1,1,2))/2);
xmin = ceil(centerX-w/2+band);
xmax = floor(centerX+w/2-band);
ymin = ceil(centerY-h/2+band);
ymax = floor(centerY+h/2-band);
cropInfo = [xmin,ymin,xmax,ymax];
[reg_grid(:,:,1),reg_grid(:,:,2)] = meshgrid(xmin:xmax,ymin:ymax);
tMap = cell(1,numel(forceField));
tMapX = cell(1,numel(forceField));
tMapY = cell(1,numel(forceField));
progressText(0,'Traction map creation:') % Create text
for ii=1:numel(forceField)
    curDispVec = displField(ii).vec;
    curDispPos = displField(ii).pos;
    curDispVec = curDispVec(~isnan(curDispVec(:,1)),:); % This will remove the warning 
    curDispPos = curDispPos(~isnan(curDispVec(:,1)),:); % This will remove the warning 
    try
        [~,iu_mat,~,~] = interp_vec2grid(curDispPos, curDispVec,[],reg_grid1);
    catch
        if ii>1
            disp(['Too many NaNs in the field in ' num2str(ii) 'th frame. Assigining values in the ' num2str(ii-1) 'th frame ...'])
            tMap{ii} = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;
            tMapX{ii} = tmat(:,:,1);
            tMapY{ii} = tmat(:,:,2);
            continue
        else
            disp(['Too many NaNs in the field in ' num2str(ii) 'th frame. Assigining zeros ...'])
            tMap{ii} = zeros(size(reg_grid(:,:,1)));
            tMapX{ii} = zeros(size(reg_grid(:,:,1)));
            tMapY{ii} = zeros(size(reg_grid(:,:,1)));
            continue
        end
    end
%     [~,iu_mat,~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid1);
    [~,if_mat,~,~] = interp_vec2grid(forceField(ii).pos, forceField(ii).vec,[],reg_grid1);
    pos = [reshape(reg_grid1(:,:,1),[],1) reshape(reg_grid1(:,:,2),[],1)]; %dense
    disp_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 
    force_vec = [reshape(if_mat(:,:,1),[],1) reshape(if_mat(:,:,2),[],1)]; 
    [~,tmat, ~, ~] = interp_vec2grid(pos+disp_vec, force_vec,[],reg_grid); %1:cluster size
    tMap{ii} = (tmat(:,:,1).^2 + tmat(:,:,2).^2).^0.5;
    tMapX{ii} = tmat(:,:,1);
    tMapY{ii} = tmat(:,:,2);
    progressText(ii/numel(forceField))
end

