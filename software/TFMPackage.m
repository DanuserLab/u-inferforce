classdef TFMPackage < Package
    % Main class of the TFM package
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
    
    % Sebastien Besson, Aug 2011
    
    methods
        function obj = TFMPackage(owner,varargin)
            % Constructor of class TFMPackage
            
            if nargin == 0
                super_args = {};
            else
                % Check input
                ip =inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieObject'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                
                super_args{1} = owner;
                super_args{2} = [outputDir  filesep 'TFMPackage'];
            end
            
            % Call the superclass constructor
            obj = obj@Package(super_args{:});
        end
        
        
        function [status processExceptions] = sanityCheck(obj,varargin) 
            
            % Check that the at least one channel has a valid psf sigma
            hasPsfSigma = any(~cellfun(@isempty, {obj.owner_.channels_.psfSigma_}));
            if ~hasPsfSigma
                error(['Missing standard deviation of the theoretical point-spread function! '...
                    'Please fill the numerical aperture, pixel size and'...
                    ' emission wavelengths for the beads channel!']);            
            end
            
            % Check that the movie has a frame rate
            if isempty(obj.owner_.timeInterval_)
                error('Missing frame rate! Please fill the time interval!');            
            end
            [status, processExceptions] = sanityCheck@Package(obj,varargin{:});
        end

    end
    
    methods (Static)
        function name = getName()
            name = 'Traction Force Microscopy';
        end
        
        function m = getDependencyMatrix(i,j)
  
            m = [0 0 0 0 0;   %1 StageDriftCorrectionProcess [optional]
                 2 0 0 0 0;   %2 DisplacementFieldCalculationProcess
                 0 1 0 0 0;   %3 DisplacementFieldCorrectionProcess [optional]
                 0 1 2 0 0;   %4 ForceFieldCalculationProcess
                 0 1 2 1 0];  %5 StrainEnergyCalculationProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = tfmPackageGUI(varargin{:});
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            TFMProcConstr = {
                @EfficientSubpixelRegistrationProcess,...
                @DisplacementFieldCalculationProcess,...
                @DisplacementFieldCorrectionProcess,...
                @ForceFieldCalculationProcess,...
                @StrainEnergyCalculationProcess};
            
            if nargin==0, index=1:numel(TFMProcConstr); end
            procConstr=TFMProcConstr(index);
        end
        function classes = getProcessClassNames(index)
            TFMClasses = {
                'StageDriftCorrectionProcess',...
                'DisplacementFieldCalculationProcess',...
                'DisplacementFieldCorrectionProcess',...
                'ForceFieldCalculationProcess',...
                'StrainEnergyCalculationProcess'};
            if nargin==0, index=1:numel(TFMClasses); end
            classes=TFMClasses(index);
        end
    end
   
end