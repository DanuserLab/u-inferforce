classdef ForceFieldCalculationProcess < DataProcessingProcess
    % Concrete process for calculating a force field
    %
    % Sebastien Besson, Aug 2011
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
    properties (SetAccess = protected)  
        tMapLimits_
        dELimits_
        distBeadMapLimits_
    end
    
    methods
        function obj = ForceFieldCalculationProcess(owner,varargin)
            
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
                super_args{2} = ForceFieldCalculationProcess.getName;
                super_args{3} = @calculateMovieForceField;
                if isempty(funParams)
                    funParams=ForceFieldCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            
            obj = obj@DataProcessingProcess(super_args{:});
            
        end
        
        function status = checkChannelOutput(obj,varargin)
            
            status = logical(exist(obj.outFilePaths_{1},'file'));
            
        end
        
        function varargout = loadChannelOutput(obj,varargin)
            
            outputList = {'forceField','tMap','forceFieldShifted','forceFieldShiftedColor'};
            ip =inputParser;
            ip.addRequired('obj',@(x) isa(x,'ForceFieldCalculationProcess'));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
%             ip.addOptional('iOut',1,@isnumeric);
%             ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) ismember(x,1:obj.owner_.nFrames_));
%             ip.addParamValue('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.addParameter('output',outputList{1},@(x) all(ismember(x,outputList)));
            ip.parse(obj,varargin{:})
            iFrame = ip.Results.iFrame;
%             iOut = ip.Results.iOut;
            
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            if strcmp(output,outputList{1}) 
                iOut=1;
                s = load(obj.outFilePaths_{iOut},output{1});
            elseif strcmp(output,outputList{2})
                iOut=2;
                s = load(obj.outFilePaths_{iOut},output{1});
            elseif strcmp(output,outputList{3})
                iOut=1;
                s = load(obj.outFilePaths_{iOut},output{1});
            elseif strcmp(output,outputList{4})
                iOut=1;
                output = outputList{3};
                s = load(obj.outFilePaths_{iOut},outputList{3});
%             elseif strcmp(output,outputList{5})
%                 [OutputDirectory,tMapFolder] = fileparts(obj.outFilePaths_{2});
%                 % Set up the output directories
%                 outputDir = fullfile(OutputDirectory,tMapFolder);
%                 outFileTMap=@(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
            end
            if ischar(output), output = {output}; end
            
%             if numel(iFrame)>1,
            varargout{1}=s.(output{1})(iFrame);
%             else
%                 varargout{1}=s.(output{1});
%             end
        end
                
        function h=draw(obj,varargin)
            % Function to draw process output
            
            outputList = obj.getDrawableOutput();
            drawLcurve = any(strcmpi('lcurve',varargin));
            rendertMap = any(strcmpi('tMap',varargin) | ...
                strcmpi('dErrMap',varargin) | strcmpi('distBeadMap',varargin) );
            if drawLcurve %Lcurve
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addParameter('output',outputList(1).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{:})
                data=obj.outFilePaths_{4,1};
            elseif rendertMap % forceMap
                % Input parser
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addRequired('iChan',@isnumeric);
                ip.addRequired('iFrame',@isnumeric);
                ip.addParameter('output',outputList(2).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{1:end})
                iFrame=ip.Results.iFrame;
                data=obj.loadChannelOutput('iFrame',iFrame,'output',ip.Results.output);
                if iscell(data), data = data{1}; end
            else % forcefield
                % Input parser
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addRequired('iFrame',@isnumeric);
                ip.addParameter('output',outputList(1).var,...
                    @(x) any(cellfun(@(y) isequal(x,y),{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{1},varargin{2:end})
                iFrame=ip.Results.iFrame;
                
                data=obj.loadChannelOutput('iFrame',iFrame,'output',ip.Results.output);
            end
            iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
            if ~isempty(outputList(iOutput).formatData)
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
        end
        
        function setTractionMapLimits(obj,tMapLimits)
            obj.tMapLimits_ = tMapLimits;
        end
        function setDisplErrMapLimits(obj,dELimits)
            obj.dELimits_ = dELimits;
        end
        function setDistBeadMapLimits(obj,distBeadMapLimits)
            obj.distBeadMapLimits_ = distBeadMapLimits;
        end
        
        function output = getDrawableOutput(obj)
            output(1).name='Force  field';
            output(1).var='forceField';
            output(1).formatData=@(x) [x.pos x.vec(:,1)/nanmean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5) x.vec(:,2)/nanmean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5)];
            output(1).type='movieOverlay';
%             output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color','r');
            output(1).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',[75/255 0/255 130/255]);
            
            output(2).name='Traction map';
            output(2).var='tMap';
            output(2).formatData=[];
            output(2).type='image';
            output(2).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units',obj.getUnits,'CLim',obj.tMapLimits_);

            output(3).name='Force field shifted';
            output(3).var='forceFieldShifted';
            output(3).formatData=@(x) [x.pos x.vec(:,1)/mean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5) x.vec(:,2)/mean((x.vec(:,1).^2+x.vec(:,2).^2).^0.5)];
            output(3).type='movieOverlay';
            output(3).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',[175/255 30/255 230/255]);

            output(4).name='Force field shifted (c)';
            output(4).var='forceFieldShiftedColor';
            output(4).formatData=@(x) [x.pos x.vec(:,1) x.vec(:,2)];
            output(4).type='movieOverlay';
            output(4).defaultDisplayMethod=@(x) VectorFieldDisplay('Colormap',jet,'Linewidth',1);
            
            if ~strcmp(obj.funParams_.method,'FTTC')

                output(5).name='Lcurve';
                output(5).var='lcurve';
                output(5).formatData=[];
                output(5).type='movieGraph';
                output(5).defaultDisplayMethod=@FigFileDisplay;


                %% TODO -- Ensure outputs are generated and available for display
                % output(6).name='Prediction Err map';
                % output(6).var='dErrMap';
                % output(6).formatData=[];
                % output(6).type='image';
                % output(6).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                %     'Colorbar','on','Units',obj.getUnits,'CLim',obj.dELimits_);

                % output(7).name='Map of distance to bead';
                % output(7).var='distBeadMap';
                % output(7).formatData=[];
                % output(7).type='image';
                % output(7).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                %     'Colorbar','on','Units',obj.getUnits,'CLim',obj.distBeadMapLimits_);


            end                
        end
        
        
    end
    methods (Static)
        function name =getName()
            name = 'Force Field Calculation';
        end
        function h = GUI()
            h= @forceFieldCalculationProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'forceField'];
            funParams.YoungModulus = 8000;
            funParams.PoissonRatio = .5;
            funParams.method = 'FastBEM';
            funParams.meshPtsFwdSol = 4096;
            funParams.regParam=1e-4;
            funParams.solMethodBEM='1NormReg';
            funParams.basisClassTblPath='';
            funParams.LcurveFactor=10;
            funParams.thickness=34000;
            funParams.useLcurve=true;
            funParams.lastToFirst=false;
            funParams.lcornerOptimal='optimal';
        end
        function units = getUnits(varargin)
            units = 'Traction (Pa)';
        end
    end
end