classdef StageDriftCorrectionProcess < ImageProcessingProcess
    % Abstract class for a stage drift correction process
    % Feb 2017 - Now Implementing multiple concrete class options for stage shift correction
    %
    % Sebastien Besson, Sep 2011
    % Andrew R. Jamieson Feb. 2017
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
    
    methods (Access = public)
        function obj = StageDriftCorrectionProcess(owner, name, funName, funParams)
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;                
            end
            
            obj = obj@ImageProcessingProcess(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;                              
            end
            if nargin > 3
               obj.funParams_ = funParams;              
            end
        end
        
        function sanityCheck(obj)
            sanityCheck@ImageProcessingProcess(obj);
            % channelIndex = obj.funParams_.ChannelIndex;
            % psfSigma = obj.owner_.channels_(channelIndex(1)).psfSigma_;
            % assert(~isempty(psfSigma), 'MovieData:Process:sanityCheck',...
            %     ['The beads channel does not have a valid '...
            %     'standard deviation of the Gaussian point-spread function.']);
        end
        
        function h = draw(obj, varargin)
            % Function to draw process output
            outputList = obj.getDrawableOutput();

            drawRefFrame = any(strcmp(varargin,'refFrame'));     
                
            if drawRefFrame

                % Use dedicated draw method for reference frame
                ip = inputParser;
                ip.addRequired('obj');
                ip.addParameter('output',outputList(2).var,@(x) all(ismember(x,{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj,varargin{:})
                
%                 [~,iOutput] =ismember(ip.Results.output,{outputList.var});
                iOutput = find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
                data = imread(obj.outFilePaths_{2,1});

                if ~isempty(outputList(iOutput).formatData),
                    data = outputList(iOutput).formatData(data);
                end

                try
                    assert(~isempty(obj.displayMethod_{iOutput,1}));
                catch ME
                    obj.displayMethod_{iOutput,1}=...
                        outputList(iOutput).defaultDisplayMethod();
                end
                % Delegate to the corresponding method
                drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                    2*numel(fieldnames(ip.Unmatched)),1);                
                tag = ['process' num2str(obj.getIndex) '_output' num2str(iOutput)];
                h = obj.displayMethod_{iOutput}.draw(data, tag, ip.Unmatched);
                
            else
                ip = inputParser;
                ip.addRequired('obj',@(x) isa(x,'Process'));
                ip.addRequired('iChan',@isnumeric);
                ip.addOptional('iFrame',[],@isnumeric);
                ip.addParameter('output',outputList(3).var,@(x) all(ismember(x,{outputList.var})));
                ip.KeepUnmatched = true;
                ip.parse(obj, varargin{:});
 
                if strcmp('merged',ip.Results.output)
                    if numel(obj.owner_.channels_) > 1, cdim=3; else cdim=1; end
                        data = zeros([obj.owner_.imSize_ cdim]);

                    iOutput = find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));

                    for iChan = 1:numel(obj.owner_.channels_)
                        imData = obj.loadChannelOutput(iChan, ip.Results.iFrame);
                        data(:,:,iChan) = outputList(iOutput).formatData(imData);
                    end                  

                    try
                        assert(~isempty(obj.displayMethod_{iOutput,1}));
                    catch ME
                        obj.displayMethod_{iOutput,1}=...
                            outputList(iOutput).defaultDisplayMethod();
                    end

                    % Create graphic tag and delegate drawing to the display class
                    tag = ['process' num2str(obj.getIndex()) '_MergedOutput'];
                    h = obj.displayMethod_{3}.draw(data, tag, ip.Unmatched);

                else
                    
                    % Call superclass method
                    h = draw@ImageProcessingProcess(obj,varargin{1},varargin{2},...
                            varargin{3:end});
                end
            end
        end

        function output = getDrawableOutput(obj)
            output = getDrawableOutput@ImageProcessingProcess();
            output(1).name = 'Registered images';
            output(2).name = 'Reference frame';
            output(2).var = 'refFrame';
            output(2).formatData = @mat2gray;
            output(2).defaultDisplayMethod=@ImageDisplay;
            output(2).type = 'movieGraph';
            output(3).name = 'Merged';
            output(3).var = 'merged';
            output(3).formatData = @mat2gray;
            output(3).defaultDisplayMethod=@ImageDisplay;
            output(3).type = 'image';

            %%TODO ? OUTPUT #3 Add new registration display similar to imshowpair(fixed, RegisteredCell{p}, 'Scaling','joint');  
            %%TODO ? OUTPUT #4 (display transforamtions for each frame?)
        end
    end

    methods (Static)
        
        function name =getName()
            name = 'Stage Drift Correction';
        end
        function h = GUI()
            h = @abstractProcessGUI;
        end
        function procClasses = getConcreteClasses()
            procClasses = ...
                {...
                 @BeadTrackingCorrectionProcess;
                 @EfficientSubpixelRegistrationProcess;
                };
            procClasses = cellfun(@func2str, procClasses, 'Unif', 0);
        end
    end
end