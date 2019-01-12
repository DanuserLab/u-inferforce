function h= quiverColormap(x,y,u,v,varargin)
% QUIVERCOLORMAP display velocity vectors as color-coded arrows
%
% quiverColormap uses the basic infrastructure of Matlab built-in quiver
% plot. In top of that, vectors are indexed by magnitude and displayed as
% color-coded according to a colormap. Arrowheads are rescaled to be
% consistent with a simple quiver.
% Note:
% 1- the autoscale function has to be disabled by this function.
% 2- varargin must be input as param/value paris. To include the linespec,
% the ideal case would be to have access to the quiverparseargs private
% function.
%
% Input:
%
%    x,y,u,v: see quiver documentation
%
%    Additional arguments: see quiver documentation
%    Two additional parameter/values can be supplied
%
%    'Colormap' (parameter/value): the colormap to use to draw the vectors.
%    Default is the jet colormap.
%
%    'CLim' : the magnitude limits to classify the vectors by magnitude. If
%    not input, the minimum and maximum of the vectors amplitude will be
%    used.
%     
% Output:
%      
%    h : array of handles for the quivergroup objects
%
% Inspired by quiverc.m, adjust_quiver_arrowhead_size.m, quiverclr.m
% Sebastien Besson, March 2012
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

% Parse input parameters
ip=inputParser;
ip.addParamValue('Colormap',jet,@(x) isnumeric(x) && size(x,2)==3);
ip.addParamValue('CLim',[],@(x) isvector(x) || isempty(x));
ip.KeepUnmatched=true;
ip.parse(varargin{:});

% Create non-color coded quiverplot to save the good arrow heads
quiverArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
h=quiver(x,y,u,v,0,quiverArgs{:});

% Retrieve the x and y data for the arrow head and delete the quivergroup
l=get(h,'Children');
if isempty(l)
    headData=get(h,{'XData','YData'});
else
    headData=get(l(2),{'XData','YData'});
end
delete(h);
hold on;

% Removing NaNs
nanIndex=isnan(u);
x(nanIndex)=[]; y(nanIndex)=[]; u(nanIndex)=[]; v(nanIndex)=[]; 

% Index vectors per magnitude
nColors = size(ip.Results.Colormap,1);
intensity= (u.^2+v.^2).^(1/2);
vColor=floor(scaleContrast(intensity,ip.Results.CLim,[1 nColors]));
vColor(vColor<1)=1;
vColor(vColor>nColors)=nColors;
vIndex= unique(vColor);

% Create array of quiverplots
for i=1:numel(vIndex)
    idx = find(vColor==vIndex(i));
    h(i) = quiver(x(idx),y(idx),u(idx),v(idx),0,quiverArgs{:},...
        'Color',ip.Results.Colormap(vIndex(i),:));
    
    % Set arrowhead x and y-data
    l=get(h(i),'Children');
    if isempty(l)
        headIdx = arrayfun(@(x) 4*(x-1)+1:4*x,idx,'Unif',false);
        headIdx=horzcat(headIdx{:});
%         set(h(i),'XData',headData{1}(headIdx),'YData',headData{2}(headIdx));
    else
        headIdx = arrayfun(@(x) 4*(x-1)+1:4*x,idx,'Unif',false);
        headIdx=horzcat(headIdx{:});
        set(l(2),'XData',headData{1}(headIdx),'YData',headData{2}(headIdx));
    end
end

% SB: the following piece of code is not satisfying for two reasons
% 1- any change of the figure colormap is going to affect the colorbar. 
% Need to implement something similar to freezeColors there.
% 2 - If vectors have been scaled for display, CLim needs to be scaled as 
% well but not the colormap. I think the colorbar should be implemented 
% upstream of this function.

% ip.addParamValue('Colorbar','on',@(x) ismember(x,{'on','off'}));
% if ~strcmp(ip.Results.Colorbar,'on'), hCbar=[]; return; end
% % Re-format the colorbar
% hCbar=colorbar;
% 
% %set(h,'ylim',[1 length(map)]);
% ylim = get(hCbar,'ylim') ;
% yal = linspace(ylim(1),ylim(2),10) ;
% set(hCbar,'ytick',yal);
% 
% % Create the yticklabels
% if isempty(ip.Results.CLim),
%     zlim = [min(intensity(:)) max(intensity(:))];
% else
%     zlim = ip.Results.CLim;
% end
% ytl=linspace(zlim(1),zlim(2),10);
% %set(h,'ytick',ytl);
% s=char(10,4);
% for i=1:10
%     if min(abs(ytl)) >= 0.001
%         B=sprintf('%-4.3f',ytl(i));
%     else
%         B=sprintf('%-3.1E',ytl(i));
%     end
%     s(i,1:length(B))=B;
% end
% set(hCbar,'yticklabel',s);
% grid on
