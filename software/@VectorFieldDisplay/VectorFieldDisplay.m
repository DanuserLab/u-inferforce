classdef VectorFieldDisplay < MovieDataDisplay
    %Conrete class for displaying flow
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
    properties
        Color='k';  
        vectorScale=1;
        Linewidth = 1;
        Linestyle = '-';
        Colormap=[];
        CLim = [];
    end
    methods
        function obj=VectorFieldDisplay(varargin)
            obj@MovieDataDisplay(varargin{:});
        end
        
        function h=initDraw(obj,data,tag,varargin)
            if isempty(obj.vectorScale),autoscale='on'; else autoscale='off'; end
  
            if isempty(obj.Colormap)
                % Create non-color coded quiverplot to retrieve arrow heads
                % delete NaNs for easy tracking of vectors
                noNans = ~isnan(data(:,3));
                h=quiver(data(noNans,1),data(noNans, 2),obj.vectorScale*(data(noNans,3)),...
                    obj.vectorScale*(data(noNans,4)),'Autoscale',autoscale,...
                    'Linestyle',obj.Linestyle,'Linewidth',obj.Linewidth,...
                    'Color',obj.Color,varargin{:});
%                 h=quiver(data(:,1),data(:, 2),obj.vectorScale*(data(:,3)),...
%                     obj.vectorScale*(data(:,4)),'Autoscale',autoscale,...
%                     'Linestyle',obj.Linestyle,'Linewidth',obj.Linewidth,...
%                     'Color',obj.Color,varargin{:});
            else
                h=quiverColormap(data(:,1),data(:, 2),obj.vectorScale*(data(:,3)),...
                    obj.vectorScale*(data(:,4)),'Autoscale',autoscale,...
                    'Linestyle',obj.Linestyle,'Linewidth',obj.Linewidth,...
                    'Colormap',obj.Colormap,'CLim',obj.vectorScale*obj.CLim,varargin{:});
            end
               
            set(h,'Tag',tag);
        end

        function updateDraw(obj,h,data)
            tag=get(h(1),'Tag');
            delete(h);
            obj.initDraw(data,tag);
        end
    end    
    
    methods (Static)
         function params=getParamValidators()
            params(1).name='Color';
            params(1).validator=@(x)ischar(x) ||isvector(x);
            params(2).name='vectorScale';
            params(2).validator=@isscalar;
            params(3).name='Colormap';
            params(3).validator=@(x) ischar(x) || isnumeric(x);
            params(4).name='CLim';
            params(4).validator=@isvector;
            params(5).name='Linestyle';
            params(5).validator=@ischar;
            params(6).name='Linewidth';
            params(6).validator=@isscalar;
        end

        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end