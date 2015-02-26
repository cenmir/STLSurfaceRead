classdef STLSurfaceRead < matlab.mixin.Copyable
    %STLSurfaceRead Read a STL surface to matlab
    %   Detailed explanation goes here
    
    properties
        Points
        Connectivity
        VertexNormals
        FaceNormals
    end
    
    properties (Access = private)
        Triangulated = 0
        Welded = 0
        hp
    end
    
    properties (Dependent)
        EdgeColor
        FaceColor
    end
    
    
    methods
        function S = STLSurfaceRead(filename)
            fv = S.stlread(filename);
            warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId')
            T = triangulation(fv.faces,fv.vertices);
            warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId')
            S.Points = fv.vertices;
            S.Connectivity = fv.faces;
            S.VertexNormals = vertexNormal(T);
            S.FaceNormals = faceNormal(T);
        end
        
        function S = TriangulateSurface(S)
            if  S.Triangulated || S.Welded
                return
            end
            
            X = S.Points;
            NF = S.FaceNormals;
            
            nTri = size(X,1)/3;
            tri = zeros(nTri, 3);
            itri = 1;
            for iP = 1:3:size(X,1)
                if itri == 1
                    tri(1,1:3) = [1,2,3];
                    xe = X(tri(itri,:),:);
                    n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
                    n = n/norm(n);
                    
                    if dot(NF(itri,:),n) < 0
                        tri(1,1:3) = [1,3,2];
                    end
                    %         tri(itri,:)
                    
                    itri = itri+1;
                else
                    % the local triangle node index is set to 1
                    tind = 1;
                    for iPnt = iP:iP+2
                        PrevPoints = [1:iPnt-1]';
                        SearchInds = PrevPoints;
                        
                        searchSpace = X(SearchInds,:);
                        
                        % Chose a point in the current set of trinagle points to measure the distance from
                        Xp = X(iPnt,:);
                        
                        % Choose a set of points to measure the distance to. The set of
                        % points is defined by the searchSpace.
                        P = ones(size(searchSpace,1),1)*Xp;
                        
                        % The distance vectors between the chosen point and the searchSpace
                        % points.
                        D = searchSpace-P;
                        
                        % Eucledian distance
                        NDist = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);
                        
                        % Find the indices to the point that are identical (with room
                        % for numerical error)
                        indt = SearchInds(find(NDist<eps*100,1 ));
                        
                        % if duplicate points exist; set the local triangle index to
                        % the index of the searchSpace and increment the triangle
                        % index.
                        % If no duplicate points are found; set the local triangle
                        % index to the index of the current triangle set point
                        if ~isempty(indt)
                            tri(itri,tind) = indt;
                            tind = tind +1;
                        else
                            tri(itri,tind) = iPnt;
                            tind = tind +1;
                        end
                    end
                    xe = X(tri(itri,:),:);
                    
                    n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
                    n = n/norm(n);
                    
                    if dot(NF(itri,:),n) < 0
                        tri(1,1:3) = [1,3,2];
                    end
                    %         tri(itri,:)
                    
                    itri = itri+1;
                    
                end
            end
            
            S.Connectivity = tri;
            
            warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId')
            T = triangulation(tri,X);
            warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId')
            S.VertexNormals = vertexNormal(T);
            S.FaceNormals = faceNormal(T);
            
            S.Triangulated = 1;
        end
        
        function S = WeldPoints(S)
            
            if ~S.Welded
                if ~S.Triangulated 
                    S.TriangulateSurface;
                end
            else
                return
            end
            
            X = S.Points;
            tri = S.Connectivity;
            
            mt = 0;
            map = zeros(size(X,1),2);
            c = 0;
            for iel = 1:size(tri)
                itri = tri(iel,:);
                
                if iel == 1
                    for it = itri
                        mt=mt+1;
                    end
                    map(1:3,1)=[1:3]';
                    map(1:3,2)=[1:3]';
                    c = 3;
                    mt = 4;
                else
                    for indt=1:3
                        it = itri(indt);
                        if it >= mt && ~any(it==map(:,1))
                            map(c+1,:) = [it,mt];
                            mt=mt+1;
                            c=c+1;
                        else
                        end
                        
                    end
                    
                end
            end
            map = map(all(map~=0,2),:);
            
            for im = 1:size(map,1)
                tri(tri==map(im,1)) = map(im,2);
            end
            X =  X(map(:,1),:);
            S.Points = X
            S.Connectivity = tri;
            
            T = triangulation(tri, X);
            S.VertexNormals = vertexNormal(T);
            S.FaceNormals = faceNormal(T);
            
            S.Welded = 1;
            
        end
        
        function VizSurface(S, varargin)
            if ~exist('xfigure','file') == 2
                S.RequiredFileMissing('xfigure', 'https://raw.githubusercontent.com/cenmir/xfigure/master/xfigure.m')
            end
            if ~exist('xfigure_KPF','file') == 2
                S.RequiredFileMissing('xfigure_KPF', 'https://raw.githubusercontent.com/cenmir/xfigure/master/xfigure_KPF.m')
            end
            
            xfigure;
            axis equal;
            S.hp = patch('Faces',S.Connectivity,'Vertices',S.Points, 'FaceColor','c', 'Facelighting','gouraud');
            S.EdgeColor = 'k';
            S.FaceColor = 'c';
            light
            view(3)
            if S.isenabled(varargin, 'FaceNormals')
                T = triangulation(S.Connectivity, S.Points);
                XC = incenter(T);
                FN = S.FaceNormals;
                hold on;
                quiver3(XC(:,1),XC(:,2),XC(:,3),FN(:,1),FN(:,2),FN(:,3),'Color','b')
                hold off;
            end
            
            if S.isenabled(varargin, 'VertexNormals')
%                 T = triangulation(S.Connectivity, S.Points);
%                 XC = incenter(T);
                X = S.Points;
                FV = S.VertexNormals;
                hold on;
                quiver3(X(:,1),X(:,2),X(:,3),FV(:,1),FV(:,2),FV(:,3),'Color','r')
                hold off;
            end
            
        end
        
        function EdgeColor = get.EdgeColor(S)
            EdgeColor = S.hp.EdgeColor;
        end
        function set.EdgeColor(S,C)
           S.hp.EdgeColor = C;
        end  
        
        function FaceColor = get.FaceColor(S)
            FaceColor = S.hp.FaceColor;
        end
        function set.FaceColor(S,C)
           S.hp.FaceColor = C;
        end  
        
    end
    
    methods (Access = private)
        function varargout = stlread(S,file)
            % STLREAD imports geometry from an STL file into MATLAB.
            %    FV = STLREAD(FILENAME) imports triangular faces from the binary
            %    STL file idicated by FILENAME, and returns the patch struct FV, with fields
            %    'faces' and 'vertices'.
            %
            %    [F,V] = STLREAD(FILENAME) returns the faces F and vertices V separately.
            %
            %    [F,V,N] = STLREAD(FILENAME) also returns the face normal vectors.
            %
            %    The faces and vertices are arranged in the format used by the PATCH plot
            %    object.
            
            % Copyright 2011 The MathWorks, Inc.
            
            if ~exist(file,'file')
                error(['File ''%s'' not found. If the file is not on MATLAB''s path' ...
                    ', be sure to specify the full path to the file.'], file);
            end
            
            fid = fopen(file,'r');
            if ~isempty(ferror(fid))
                error(lasterror); %#ok
            end
            
            M = fread(fid,inf,'uint8=>uint8');
            fclose(fid);
            
            [f,v,n] = S.stlbinary(M);
            
            varargout = cell(1,nargout);
            switch nargout
                case 2
                    varargout{1} = f;
                    varargout{2} = v;
                case 3
                    varargout{1} = f;
                    varargout{2} = v;
                    varargout{3} = n;
                otherwise
                    varargout{1} = struct('faces',f,'vertices',v);
            end
            
        end
        
        function [F,V,N] = stlbinary(S,M)
            
            F = [];
            V = [];
            N = [];
            
            if length(M) < 84
                error('MATLAB:stlread:incorrectFormat', ...
                    'Incomplete header information in binary STL file.');
            end
            
            % Bytes 81-84 are an unsigned 32-bit integer specifying the number of faces
            % that follow.
            numFaces = typecast(M(81:84),'uint32');
            %numFaces = double(numFaces);
            if numFaces == 0
                warning('MATLAB:stlread:nodata','No data in STL file.');
                return
            end
            
            T = M(85:end);
            F = NaN(numFaces,3);
            V = NaN(3*numFaces,3);
            N = NaN(numFaces,3);
            
            numRead = 0;
            while numRead < numFaces
                % Each facet is 50 bytes
                %  - Three single precision values specifying the face normal vector
                %  - Three single precision values specifying the first vertex (XYZ)
                %  - Three single precision values specifying the second vertex (XYZ)
                %  - Three single precision values specifying the third vertex (XYZ)
                %  - Two unused bytes
                i1    = 50 * numRead + 1;
                i2    = i1 + 50 - 1;
                facet = T(i1:i2)';
                
                n  = typecast(facet(1:12),'single');
                v1 = typecast(facet(13:24),'single');
                v2 = typecast(facet(25:36),'single');
                v3 = typecast(facet(37:48),'single');
                
                n = double(n);
                v = double([v1; v2; v3]);
                
                % Figure out where to fit these new vertices, and the face, in the
                % larger F and V collections.
                fInd  = numRead + 1;
                vInd1 = 3 * (fInd - 1) + 1;
                vInd2 = vInd1 + 3 - 1;
                
                V(vInd1:vInd2,:) = v;
                F(fInd,:)        = vInd1:vInd2;
                N(fInd,:)        = n;
                
                numRead = numRead + 1;
            end
            
        end
        
        function rv = isenabled(S,mode, varargin)
            % ISENABLED  Checks if mode exists in the cell-array varargin.
            % isenabled(mode,varargin{:}) return true or false.
            if nargin < 1
                error('No arguments')
            end
            varargin = varargin{:};
            
            ind = find(strcmpi(varargin,mode), 1);
            if ~isempty(ind)
                rv = 1;
            else
                rv = 0;
            end
        end
        
        function RequiredFileMissing(S,filename, RemoteDestination)
            % If we're going to download a whole bunch of files, it is better to set
            % RequiredFilesDir to be a global and not have to ask the user to
            % specify a destination folder for every file...
            GetABunchOfFiles = 1;
            
            if GetABunchOfFiles
                global RequiredFilesDir
            else
                RequiredFilesDir = [];
            end
            
            warning([filename,' is missing!'])
            disp([filename,' is missing!'])
            disp(['Trying to download ',filename,' from ',RemoteDestination])
            
            
            if isempty(RequiredFilesDir)
                scriptPath = mfilename('class');
                [ScriptDir,~,~] = fileparts(scriptPath);
                DestDir = uigetdir(ScriptDir,['Select where to save ',filename,'. Make sure its either the script directory or a directory on the Path.']);
                if DestDir == 0
                    error(['Failed to select folder, failed to install ',filename])
                end
                
                RequiredFilesDir = DestDir;
            end
            DestFile= [RequiredFilesDir,'/',filename];
            
            % Download the RemoteFile and save it to DestFile
            websave(DestFile,RemoteDestination);
            
            % Give up to 10 seconds for the file to show up, otherwise send error
            % message.
            tic
            while 1
                if exist('xfigure','file') == 2
                    break
                end
                pause(0.1)
                t1 = toc;
                if t1 > 10
                    error(['Failed to download ',filename,'! Timeout.'])
                end
            end
        end
        
        
    end
    
    %Hide some of the inherited methods from handle
    methods(Hidden)
        function lh = addlistener(varargin)
            lh = addlistener@handle(varargin{:});
        end
        function notify(varargin)
            notify@handle(varargin{:});
        end
        function delete(varargin)
            delete@handle(varargin{:});
        end
        function Hmatch = findobj(varargin)
            Hmatch = findobj@handle(varargin{:});
        end
        function p = findprop(varargin)
            p = findprop@handle(varargin{:});
        end
        function TF = eq(varargin)
            TF = eq@handle(varargin{:});
        end
        function TF = ne(varargin)
            TF = ne@handle(varargin{:});
        end
        function TF = lt(varargin)
            TF = lt@handle(varargin{:});
        end
        function TF = le(varargin)
            TF = le@handle(varargin{:});
        end
        function TF = gt(varargin)
            TF = gt@handle(varargin{:});
        end
        function TF = ge(varargin)
            TF = ge@handle(varargin{:});
        end
    end
end

