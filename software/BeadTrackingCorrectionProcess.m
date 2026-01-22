classdef BeadTrackingCorrectionProcess < StageDriftCorrectionProcess
    % Concrete class for bead-tracking based stage drift correction process
    %
    % Sebastien Besson, Sep 2011
    % Andrew R. Jamieson Feb. 2017
%
% Copyright (C) 2026, Danuser Lab - UTSouthwestern 
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
        function obj = BeadTrackingCorrectionProcess(owner,varargin)
            
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
                
                % Constructor for BeadTrackingCorrectionProcess

                super_args{1} = owner;
                super_args{2} = BeadTrackingCorrectionProcess.getName;
                super_args{3} = @correctMovieStageDrift;
                if isempty(funParams)
                    funParams=BeadTrackingCorrectionProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@StageDriftCorrectionProcess(super_args{:});
        end
        
        function sanityCheck(obj)
            sanityCheck@ImageProcessingProcess(obj);
            channelIndex = obj.funParams_.ChannelIndex;
            psfSigma = obj.owner_.channels_(channelIndex(1)).psfSigma_;
            assert(~isempty(psfSigma), 'MovieData:Process:sanityCheck',...
                ['The beads channel does not have a valid '...
                'standard deviation of the Gaussian point-spread function.']);
        end

        function h=draw(obj,varargin)
            % Function to draw process output
            
            outputList = obj.getDrawableOutput();
            drawGraph = any(strcmp(varargin,'x-flow') | strcmp(varargin,'y-flow') |...
                strcmp(varargin,'refFrame'));
            
            
            if drawGraph
                % Use dedicated draw method for plotting flow histograms
                ip = inputParser;
                ip.addRequired('obj');
                ip.addParamValue('output',outputList(2).var,@(x) all(ismember(x,{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{:})
                
                [~,iOutput] =ismember(ip.Results.output,{outputList.var});
                if regexp(outputList(iOutput).var,'(.+)flow','once')
                    s=load(obj.outFilePaths_{3,1},'flow');
                    data=s.flow;
                else
%                     data=imread(obj.funParams_.referenceFramePath);
                    data=imread(obj.outFilePaths_{2,1});
                end
                
                if ~isempty(outputList(iOutput).formatData),
                    data=outputList(iOutput).formatData(data);
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
            else
                % Call superclass method
                h=draw@StageDriftCorrectionProcess(obj,varargin{1},varargin{2},...
                    varargin{3:end});
            end
        end
        
        function output = getDrawableOutput(obj, varargin)
            % Rename default registration output, add bead tracking flow
            output = getDrawableOutput@StageDriftCorrectionProcess(obj);
            output(4).name='Flow along x-axis';
            output(4).var='x-flow';
            output(4).formatData=@(x)getFlow(x,1);
            output(4).type='movieGraph';
            output(4).defaultDisplayMethod=@FlowHistogramDisplay;
            output(5).name='Flow along y-axis';
            output(5).var='y-flow';
            output(5).formatData=@(x)getFlow(x,2);
            output(5).type='movieGraph';
            output(5).defaultDisplayMethod=@FlowHistogramDisplay;
        end

    end
    methods (Static)
        
        function name = getName()
            name = 'Bead Tracking Drift Correction';
        end
        function h = GUI()
            h= @BeadTrackingCorrectionProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'stageDriftCorrection_BeadTrack'];
            funParams.referenceFramePath = '';
            funParams.minCorLength = 51;
            funParams.maxFlowSpeed =5;
            funParams.alpha = .05;
            funParams.cropROI=[1 1 owner.imSize_(end:-1:1)];
            funParams.doPreReg=1;
        end
    end
end

function data=getFlow(data,i)

data=cellfun(@(x) x(:,i+2)-x(:,i),data,'UniformOutput',false);
end