classdef CziReader < handle
    %CZIReader Reads czis using libCZI.
    %   Tiles, positions in XY, subsets in XY, complex, and RGB images are
    %   not supported at the moment.
        
    properties (SetAccess = private)
        FileName
        ImageSize
    end
    properties
        ZenOrientation = false
    end
    properties (Access = private)
        DimensionOrder = 'ZHTSC'
        MetadataXml
        PixelType
        HigherDimensionOrder
        HigherDimensionIndices
        CoordinateIndex
        ZoomIndex
        SubblockBounds
    end
    
    methods
        function obj = CziReader(filename, varargin)
            % Construct an instance of CziReader
            %   filename: czi file name to read
            %   optional: desired dimension order
            obj.FileName = filename;
            narginchk(1,2);
            if nargin == 2
                obj.DimensionOrder = erase(upper(varargin{1}),{'X','Y'});
            end
            
            % initial open to get bounds and metadata
            czih = MEXlibCZI('Open', filename);
            try
                cziinfo = MEXlibCZI('GetInfo', czih);
                if (cziinfo.subblockcount == 0)
                    MEXlibCZI('Close', czih);
                    error(['File' filename ' is empty']);
                end
                obj.SubblockBounds = cell(1, cziinfo.subblockcount);

                bounds.X = cziinfo.boundingBox(3);
                bounds.Y = cziinfo.boundingBox(4);
                for dim = obj.DimensionOrder
                    if isfield(cziinfo.dimBounds, dim)
                        bounds.(dim) = cziinfo.dimBounds.(dim)(2);
                    end
                end
                
                obj.ImageSize = bounds;
                obj.MetadataXml = MEXlibCZI('GetMetadataXml', czih);

                % data type and subblock info structure
                sh = MEXlibCZI('GetSubBlock', czih, 0);
                sinfo = MEXlibCZI('GetInfoFromSubBlock', czih, sh);
                for i=1:length(sinfo)
                    % if strcmp(sinfo(i).property, 'Pixeltype')
                    %     pixelMap = struct('gray8','uint8','gray16','uint16','gray32float','single'); % looks like complex is not implemented by MEXlibCZI
                    %     obj.PixelType = pixelMap.(sinfo(i).value);
                        obj.PixelType = sinfo.Pixeltype;
                    % elseif strcmp(sinfo(i).property, 'Coordinate')
                    %     obj.CoordinateIndex = i;
                    %     cs = sinfo(i).value;
                    %     alldims = cs(regexp(cs,'\D'));
                    % 
                    %     for dim = obj.DimensionOrder
                    %         di = strfind(alldims, dim);
                    %         if ~isempty(di) && isfield(bounds, dim) && bounds.(dim) > 1
                    %             obj.HigherDimensionIndices.(dim) = di;
                    %         end
                    %     end
                    % elseif strcmp(sinfo(i).property, 'Zoom')
                    %     obj.ZoomIndex = i;
                    % end
                end
                MEXlibCZI('ReleaseSubBlock', czih, sh);

                if isstruct(obj.HigherDimensionIndices)
                    higherdims = fieldnames(obj.HigherDimensionIndices);
                    obj.HigherDimensionOrder = [higherdims{:}];
                end

                MEXlibCZI('Close', czih);
            catch me
                MEXlibCZI('Close', czih);
                rethrow(me);
            end
        end

        function data = Read(obj, varargin)
            % Reads data from a czi.
            %   optional argument: slice in higher dimensions
            filter = struct([]);
            if nargin == 2
                filter = varargin{1};
            end
            
            dataSize = [];
            for dim = obj.HigherDimensionOrder
                if isfield(filter, dim)
                    dimSize = 1;
                else
                    dimSize = obj.ImageSize.(dim);
                end
                dataSize = [dataSize dimSize]; %#ok<AGROW>
            end
            
            wh = [obj.ImageSize.X obj.ImageSize.Y];
            if (obj.ZenOrientation)
                wh = fliplr(wh);
            end
            sizeXY = int64(prod(wh));
            % data = zeros([wh dataSize], obj.PixelType);
            data = zeros([wh dataSize], 'single'); % changed
            
            czih = MEXlibCZI('Open', obj.FileName);
            try
                cziinfo = MEXlibCZI('GetInfo', czih);

                for i = 1:cziinfo.subblockcount
                    sh = -1;
                    if isempty(obj.SubblockBounds{i})
                        sh = MEXlibCZI('GetSubBlock', czih, i-1);
                        sinfo = MEXlibCZI('GetInfoFromSubBlock', czih, sh);
                        % cs = sinfo(obj.CoordinateIndex).value;
                        % zoom = sinfo(obj.ZoomIndex).value;
                        cs = sinfo.Coordinate;
                        zoom = sinfo.Zoom;
                        if zoom < 1
                            obj.SubblockBounds{i} = 'Pyramid';
                            MEXlibCZI('ReleaseSubBlock', czih, sh);
                        else
                            obj.SubblockBounds{i} = cs;
                        end
                    end
                    
                    if strcmp(obj.SubblockBounds{i}, 'Pyramid')
                        continue;
                    end
                    
                    start = str2double(regexp(obj.SubblockBounds{i},'\d*','match'));
                    offset = int64(0);
                    factor = sizeXY;
                    inRange = 1;
                    for dim = obj.HigherDimensionOrder
                        if isfield(filter, dim)
                            if start(obj.HigherDimensionIndices.(dim)) + 1 ~= filter.(dim)
                                inRange = 0;
                            end
                        else
                            offset = offset + factor * start(obj.HigherDimensionIndices.(dim));
                            factor = factor * int64(obj.ImageSize.(dim));
                        end
                    end

                    if (inRange)
                        if sh < 0
                            sh = MEXlibCZI('GetSubBlock', czih, i-1);
                        end

                        if (obj.ZenOrientation)
                            b = MEXlibCZI('GetBitmapFromSubBlock', czih, sh);
                        else
                            b = MEXlibCZI('GetBitmapFromSubBlock', czih, sh)';
                        end

                        data(offset+1:offset+sizeXY)=b(:);
                    end

                    if sh > 0
                        MEXlibCZI('ReleaseSubBlock', czih, sh);
                    end
                end

                MEXlibCZI('Close', czih);
            catch me
                MEXlibCZI('Close', czih);
                rethrow(me);
            end
        end
        
        function metadata = ReadMetadata(obj, metadataPath)
            metadata = searchinxml(obj.MetadataXml, metadataPath);
        end

    end
    
    methods(Static)
        function data = ReadAll(filename, varargin)
            if nargin == 2
                cr = CziReader(filename, varargin{1});
            else
                cr = CziReader(filename);
            end
            data = cr.Read;
        end
    end
end

function metadata = searchinxml(xmlstring, xmlpath)

    NET.addAssembly('System.Runtime');
    NET.addAssembly('System.Xml');
    import System.IO.*
    import System.Xml.*

    stringReader = StringReader(xmlstring);
    xmlReader = XmlReader.Create(stringReader);

    xmlReader.MoveToContent();
    metadataFound = true;
    
    for subpath = regexp(xmlpath,'[.]','split')
        elementNameAndIndex = regexp(subpath{1},'[[]]','split');
        elementName = elementNameAndIndex{1};
        metadataFound = metadataFound && xmlReader.ReadToDescendant(elementName);
            
        if length(elementNameAndIndex) > 1
            elementIndex = str2double(elementNameAndIndex{2});
            elementName = elementName(1:end-1); % e.g. Channels -> Channel TODO is this valid in general?
            metadataFound = metadataFound && xmlReader.ReadToDescendant(elementName);
            for i=1:elementIndex
                xmlReader.Skip;
                metadataFound = metadataFound && xmlReader.ReadToFollowing(elementName);
            end
        end
    end
    
    if metadataFound
        metadata = char(xmlReader.ReadElementContentAsString);
        metadataAsNumber = str2double(metadata);
        if ~isnan(metadataAsNumber)
            metadata = metadataAsNumber;
        end
    else
        metadata = [];
    end
    
    xmlReader.Close();
    xmlReader.Dispose();
    stringReader.Close();
    stringReader.Dispose();
end

