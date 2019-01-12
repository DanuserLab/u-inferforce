classdef FlowHistogramDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
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
    properties
        Marker = 'none';
        Linewidth = 2;
        nBins=20;
    end
    methods
        function obj=FlowHistogramDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            % Read data number
            nData=numel(data);
            dataIdx=[1 round(nData/2) nData];
            colors =hsv(numel(dataIdx));
            
            % define small and large fonts     
            tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
            sfont = {'FontName', 'Helvetica', 'FontSize', 18};
            lfont = {'FontName', 'Helvetica', 'FontSize', 22};
            
            % Generate plot
            hold on;
            h=-1*ones(numel(dataIdx),1);
            for i=1:numel(dataIdx),
                [n,x]=hist(data{dataIdx(i)},obj.nBins);
                h(i)=plot(x,n,'Color',colors(i,:),'Linewidth',obj.Linewidth);
            end
            legend(arrayfun(@(x) ['Frame ' num2str(x)],dataIdx,'UniformOutput',false),...
                'Location', 'NorthEast', tfont{:});
            xlabel('Flow (pixels/frame)',lfont{:});
            ylabel('Number',lfont{:});
            set(gca, 'LineWidth', 1.5, sfont{:});
            
            set(h,'Tag',tag);
        end
        function updateDraw(obj,h,data)
            tag = get(h(1),'Tag');
            cla(get(get(h(1),'Parent')))
            obj.initDraw(data,tag);
        end
    end    
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Marker';
            params(1).validator=@ischar;
            params(2).name='Linewidth';
            params(2).validator=@isscalar;
            params(3).name='nBins';
            params(3).validator=@isscalar;
        end
        
        function f=getDataValidator()
            f=@iscell;
        end
    end    
end