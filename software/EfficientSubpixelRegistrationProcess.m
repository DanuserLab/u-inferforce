classdef  EfficientSubpixelRegistrationProcess < StageDriftCorrectionProcess
    % Concrete class for a stage drift correction process
    %
    % Andrew R. Jamieson Feb. 2017
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
    
    methods
        function obj = EfficientSubpixelRegistrationProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Constructor for EfficientSubpixelRegistrationProcess

                super_args{1} = owner;
                super_args{2} = EfficientSubpixelRegistrationProcess.getName;
                super_args{3} = @efficientSubPixelRegistration;
                if isempty(funParams)
                    funParams = EfficientSubpixelRegistrationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@StageDriftCorrectionProcess(super_args{:});
        end
        
        
        function output = getDrawableOutput(obj, varargin)
            % Rename default registration output, add bead tracking flow
            output = getDrawableOutput@StageDriftCorrectionProcess(obj);
        end
    end
    methods (Static)
        
        function name = getName()
            name = 'Efficient Subpixel Registration';
        end
        function h = GUI()
            h = @EfficientSubpixelRegistrationProcessGUI;
        end
        
        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'stageDriftCorrection_Efficientsubpixel'];
            funParams.referenceFramePath = '';
            funParams.referenceFrameNum = 1;
            funParams.usfac = 20;
        end
    end
end