classdef DisplacementFieldCalculationProcess < ImageAnalysisProcess
    % Concrete class for a displacement field calculation process
    %
    % Sebastien Besson, Aug 2011
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
    properties (SetAccess = protected)  
        tMapLimits_
    end
    
    methods
        function obj = DisplacementFieldCalculationProcess(owner,varargin)
            
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
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = DisplacementFieldCalculationProcess.getName;
                super_args{3} = @calculateMovieDisplacementField;
                if isempty(funParams)
                    funParams=DisplacementFieldCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4}=funParams;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        
        function status = checkChannelOutput(obj,varargin)
            status = logical(exist(obj.outFilePaths_{1},'file'));
        end
        
        function sanityCheck(obj)
            sanityCheck@ImageAnalysisProcess(obj);
            channelIndex = obj.funParams_.ChannelIndex;
            assert(isscalar(channelIndex), 'lccb:Process:sanityCheck',...
                'A single bead channel must be selected to run this process');
            psfSigma = obj.owner_.channels_(channelIndex).psfSigma_;
            assert(~isempty(psfSigma), 'lccb:Process:sanityCheck',...
                ['The beads channel does not have a valid '...
                'standard deviation of the Gaussian point-spread function.']);
        end
        
        function varargout = loadChannelOutput(obj,varargin)
            
            outputList = {'displField','dMap'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'DisplacementFieldCalculationProcess'));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
            
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            iOut = cellfun(@(x) strcmp(x,output),outputList);
            s = load(obj.outFilePaths_{iOut},output{:});
            
            varargout = cell(numel(output),1);
            if numel(iFrame)>1
                for i=1:numel(output),
                    varargout{i}=s.(output{i});
                end
            else
                for i=1:numel(output),
                    varargout{i}=s.(output{i})(iFrame);
                end
            end
        end
        
        function setTractionMapLimits(obj,tMapLimits)
            obj.tMapLimits_ = tMapLimits;
        end
        
        function h=draw(obj,iFrame,varargin)
            % Function to draw process output (template method)
            
            if ~ismember('getDrawableOutput',methods(obj)), h=[]; return; end
            outputList = obj.getDrawableOutput();

            rendertMap = any(strcmpi('dMap',varargin));
            if rendertMap
                % Input parser
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
%                 ip.addRequired('iChan',@isnumeric);
                ip.addRequired('iFrame',@isnumeric);
                ip.addParamValue('output',outputList(2).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{1:end})
                iFrame=ip.Results.iFrame;
                data=obj.loadChannelOutput('iFrame',iFrame,'output',ip.Results.output);
                if iscell(data), data = data{1}; end
            else                
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addRequired('iFrame',@isnumeric);
                ip.addParamValue('output',outputList(1).var,@(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,iFrame,varargin{:})
                data=obj.loadChannelOutput(iFrame,'output',ip.Results.output);
            end
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData),
                data=outputList(iOutput,1).formatData(data);
            end
            try
                assert(~isempty(obj.displayMethod_{iOutput,1}));
            catch ME
                obj.displayMethod_{iOutput,1}=...
                    outputList(iOutput).defaultDisplayMethod();
            end

            % Delegate to the corresponding method
            tag = ['process' num2str(obj.getIndex) '_output' num2str(iOutput)];
            drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                2*numel(fieldnames(ip.Unmatched)),1);
            h=obj.displayMethod_{iOutput}.draw(data,tag,drawArgs{:});
        end
        
        function output = getDrawableOutput(obj)
            output(1).name='Displacement field';
            output(1).var='displField';
            output(1).formatData=@(x) [x.pos x.vec];
            output(1).type='movieOverlay';
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');

            output(2).name='Displacement map';
            output(2).var='dMap';
            output(2).formatData=[];
            output(2).type='image';
            output(2).defaultDisplayMethod=@(x) ImageDisplay('Colormap','jet','Colorbar','on','Units',obj.getUnits,'CLim',obj.tMapLimits_);
        end
        
    end
    methods (Static)
        function name =getName()
            name = 'Displacement Field Calculation';
        end
        function h = GUI()
            h= @displacementFieldCalculationProcessGUI;
        end
        
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1;
            funParams.OutputDirectory = [outputDir  filesep 'displacementField'];
            funParams.referenceFramePath='';
            funParams.alpha=.05;
            funParams.minCorLength = 21;
            funParams.maxFlowSpeed =20;
            funParams.highRes = true;
            funParams.mode = 'fast';
            funParams.useGrid = false;
            funParams.noFlowOutwardOnBorder = true;
            funParams.lastToFirst=false;
            funParams.addNonLocMaxBeads = false;
            funParams.trackSuccessively = false;
        end
        function units = getUnits(varargin)
            units = 'Displacement (Pix)';
        end
    end
end