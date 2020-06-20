classdef  TiffImagesReader < TiffSeriesReader
    % TiffImagesReader reads tiff images containing a single imFolder.
    % TiffImagesReader is a reader for ImageData AND ImFolder, it inherits properties
    % and methods from its superclasses, TiffSeriesReader and Reader, but
    % overwrites some methods in its superclasses to fit the special needs
    % for ImageData and ImFolder.
    %
    % TiffImageReader's properties - sizeXmax, sizeYmax, nImages, bitDepthMax, sizeZ, and filenames are overwritable.
    %
    % See also Reader, TiffSeriesReader
    %
    % Qiongjing (Jenny) Zou, June 2020
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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
        sizeXmax     % max of the length of the x dimension, horizontal, in one imFolders. If as ImageData.reader, then is a cell array of sizeXmax of all imFolders.
        sizeYmax     % max of the length of the y dimension, vertical, in one imFolders. If as ImageData.reader, then is a cell array of sizeYmax of all imFolders.
        nImages      % number of images in one imFolders. If as ImageData.reader, then is a cell array of nImages of all imFolders. Also see reader.sizeT.
        bitDepthMax     % max of the bitDepth, among all images, in one imFolders. If as ImageData.reader, then is a cell array of bitDepthMax of all imFolders.
    end
    
    methods
        
        % Constructor
        function obj = TiffImagesReader(imFolderPaths,varargin)
            obj = obj@TiffSeriesReader(imFolderPaths,varargin{:});
        end
        
        function getDimensions(obj)
            % Set/update sizeXmax, sizeYmax, nImages, bitDepthMax, sizeZ, and filenames for
            % the reader, and those properties are overwritable.
            
            % QZ sizeZ is always 1, since right now do not consider 3D iamges.
            % QZ Image dimensions and Bit dept do not need to be consistent.
            
            obj.sizeXmax = cell(obj.getSizeC(), 1);
            obj.sizeYmax = cell(obj.getSizeC(), 1);
            obj.nImages = cell(obj.getSizeC(), 1);
            obj.sizeZ = cell(obj.getSizeC(), 1);
            obj.bitDepthMax = cell(obj.getSizeC(), 1);
            for iImFol = 1 : obj.getSizeC()
                fileNames = obj.getImageFileNames(iImFol); % TiffImagesReader.filenames was set here.
                imInfo = cellfun(@(x) imfinfo([obj.paths{iImFol} filesep x]), fileNames, 'unif', 0);
                % sizeX(iChan) = unique(cellfun(@(x)(x.Width), imInfo));
                % sizeY(iChan) = unique(cellfun(@(x)(x.Height), imInfo));
                obj.sizeXmax{iImFol} = max(cellfun(@(x)(x.Width), imInfo));
                obj.sizeYmax{iImFol} = max(cellfun(@(x)(x.Height), imInfo));
                
                
                if length(fileNames)>1
                    obj.sizeZ{iImFol} = unique(cellfun(@numel, imInfo));
                    obj.nImages{iImFol} = length(fileNames);
                else % if single file, assume stack and check for # of files
                    if(obj.force3D)
                        info = imfinfo(fullfile(obj.paths{iImFol}, fileNames{1}));
                        obj.nImages{iImFol} = 1;
                        obj.sizeZ{iImFol} = numel(info);
                    else
                        info = imfinfo(fullfile(obj.paths{iImFol}, fileNames{1}));
                        obj.nImages{iImFol} = numel(info);
                        obj.sizeZ{iImFol} = 1 ;
                    end
                end
                obj.bitDepthMax{iImFol} = max(cellfun(@(x)(x.BitDepth), imInfo));
            end

            if obj.getSizeC() == 1
                obj.sizeXmax = obj.sizeXmax{:};
                obj.sizeYmax = obj.sizeYmax{:};
                obj.nImages = obj.nImages{:};
                obj.sizeZ = obj.sizeZ{:};
                obj.bitDepthMax = obj.bitDepthMax{:};
            end
        end
        
        function nImages = getNumImages(obj)
            % QZ allow adding more images to the ImFolder, and update
            % nImages
            obj.getDimensions();
            nImages = obj.nImages;
        end
        
        function filenames = getImageFileNames(obj, iImFol, iImage)
            % QZ to overwrite TiffSeriesReader.getImageFileNames
            % Allow filenames to be overwritable.
            
            obj.checkPath(iImFol);
            [files, nofExt] = imDir(obj.paths{iImFol}, true);
            assert(nofExt~=0,['No proper image files are detected in:'...
                '\n\n%s\n\nValid image file extension: tif, TIF, STK, bmp, BMP, jpg, JPG.'],obj.paths{iImFol});
            assert(nofExt==1,['More than one type of image files are found in:'...
                '\n\n%s\n\nPlease make sure all images are of same type.'],obj.paths{iImFol});
            
            obj.filenames{iImFol} = arrayfun(@(x) x.name, files, 'unif', 0);
            
            % if index has been supplied & frames are not stored in single stack
            if nargin>2 && ~obj.isSingleMultiPageTiff(iImFol)
                filenames = obj.filenames{iImFol}(iImage);
            else
                filenames = obj.filenames{iImFol};
            end
        end
        
        function I = loadImage(obj, c, t, varargin)
            % QZ to overwrite TiffSeriesReader.loadImage
            
        % loadImage reads a single image plane as a 2D, YX Matrix
        %
        % loadImage(c,t,z) reads the YX plane at (c,t,z)
        %
        % loadImage(c,t) reads the YX plane at (c,t,1)
        %
        % Note: Override if you need to overload or change how the input is
        % checked. Otherwise override loadImage_.
        %
        % Example:
        %   reader = movieData.getReader();
        %   I = reader.loadImage(1,1);
        %   imshow(I,[]);
        %
        
        % Backwards compatability before 2015/01/01:
        % Previously this function was abstract and therefore should be
        % overridden by all subclasses written prior to 2015/01/01
         
            ip = inputParser;
            ip.addRequired('c', ...
                @(c) isscalar(c) && ismember(c, 1 : obj.getSizeC()));
            % TiffSeriesReader allows for multiple t values
            ip.addRequired('t', ... 
                @(t) all(ismember(t, 1 : obj.getSizeT())));
            ip.addOptional('z', 1, ...
                @(z) isscalar(z) && ismember(z, 1 : obj.getSizeZ()) || ...
                    isempty(z));
            ip.parse(c, t, varargin{:});
                     
            z = ip.Results.z;
            if(isempty(z))
                z = 1;
            end
                      
            I = obj.loadImage_(c , t , z);
        end
        
        function sizeT = getSizeT(obj)
            if isempty(obj.nImages), obj.getDimensions(); end
            sizeT = obj.nImages;
        end
    end
    
    methods ( Access = protected )
        
        function I = loadImage_(obj, iChan, iFrame, iZ)
            % QZ to overwrite TiffSeriesReader.loadImage_
            
            if ~obj.isSingleMultiPageTiff(iChan)
                % Read individual files
                fileNames = obj.getImageFileNames(iChan, iFrame);
                
                % Initialize array
                sizeX = obj.getSizeX(); % QZ WIP continue from here!!! need to add sizeX cell for TiffImagesreader, need to write getSizeX(iFrame).
                sizeY = obj.getSizeY();
                bitDepth = obj.getBitDepth();
                I = zeros([sizeY, sizeX, numel(iFrame)], ['uint' num2str(bitDepth)]);
                
                for i=1:numel(iFrame)
                    I(:,:,i) = imread([obj.paths{iChan} filesep fileNames{i}], iZ);
                end
            else % if the channel is stored as a multi-page TIFF
                I = readtiff(fullfile(obj.paths{iChan}, obj.filenames{iChan}{1}), iFrame);
            end
        end
    end
end