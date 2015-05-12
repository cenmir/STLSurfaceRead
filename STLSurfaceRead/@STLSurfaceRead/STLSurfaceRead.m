classdef STLSurfaceRead < matlab.mixin.Copyable
    %STLSurfaceRead Read an STL surface to matlab
    %   Reads an STL and creates an object.
    %   S = STLSurfaceRead(filename)
    
    properties
        Points % Coordinates m-by-3
        Connectivity % Triangle connectivity matrix m-by-3
        VertexNormals 
        FaceNormals
    end
    
    properties (Access = private)
        Triangulated = 0
        Welded = 0
        hp
        GD
    end
    
    properties (Dependent)
        EdgeColor %Get or set edge color
        FaceColor %Get or set edge color
    end
    
    
    methods
        function S = STLSurfaceRead(filename) % Constructor
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
            % Triangulates the surface.
            if  S.Triangulated || S.Welded
                return
            end
            
            X = S.Points;
            NF = S.FaceNormals;

            tri = PrivateTriangulate(X,NF,0);
            [X,tri] = weldPoints(X,tri);
          
            %%
            S.Connectivity = tri;
            S.Points = X;
            
            warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId')
            T = triangulation(tri,X);
            warning('on', 'MATLAB:triangulation:PtsNotInTriWarnId')
            S.VertexNormals = vertexNormal(T);
            S.FaceNormals = faceNormal(T);
            
            S.Triangulated = 1;
            S.Welded = 1;
            
        end
        
        function VizSurface(S, varargin)
            if ~exist('xfigure','file') == 2
                S.RequiredFileMissing('xfigure', 'https://raw.githubusercontent.com/cenmir/xfigure/master/xfigure.m')
            end
            if ~exist('xfigure_KPF','file') == 2
                S.RequiredFileMissing('xfigure_KPF', 'https://raw.githubusercontent.com/cenmir/xfigure/master/xfigure_KPF.m')
            end
            
            [h1.fig, xfigure_This] = xfigure;
            axis equal;
            S.hp = patch('Faces',S.Connectivity,'Vertices',S.Points, 'FaceColor','c', 'Facelighting','gouraud');
            S.EdgeColor = 'k';
            S.FaceColor = 'c';
            S.GD.lights(1) = light;
            view(3)
            if S.isenabled(varargin, 'FaceNormals')
                T = triangulation(S.Connectivity, S.Points);
                XC = incenter(T);
                FN = S.FaceNormals;
                hold on;
                S.GD.FNquiver = quiver3(XC(:,1),XC(:,2),XC(:,3),FN(:,1),FN(:,2),FN(:,3),'Color','b');
                hold off;
            end
            
            if S.isenabled(varargin, 'VertexNormals')
%                 T = triangulation(S.Connectivity, S.Points);
%                 XC = incenter(T);
                X = S.Points;
                FV = S.VertexNormals;
                hold on;
                S.GD.VNquiver = quiver3(X(:,1),X(:,2),X(:,3),FV(:,1),FV(:,2),FV(:,3),'Color','r');
                hold off;
            end
            
            %Leave this be, we need the standard keypress function!
            
            h1.fig.KeyPressFcn = {@GKPF2,h1,xfigure_This,S};
            
            
            
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

function GKPF2(src,evnt,h,xfigure_This,S)
    xfigure_KPF(src, evnt, xfigure_This); %Do not touch

    if strcmp(evnt.Key,'s')
        if strcmpi(S.hp.FaceLighting,'gouraud')
            S.hp.FaceLighting = 'flat';
            set(xfigure_This.StatusBox, 'String', 'FaceLighting: flat')
        else
            S.hp.FaceLighting = 'gouraud';
            set(xfigure_This.StatusBox, 'String', 'FaceLighting: gouraud')
        end
    end
    
    if strcmpi(evnt.Character,'l')
        for ilight = S.GD.lights
            if strcmpi(ilight.Visible,'on')
                ilight.Visible = 'off';
                set(xfigure_This.StatusBox, 'String', 'Lights: off')
            else
                ilight.Visible = 'on';
                set(xfigure_This.StatusBox, 'String', 'Lights: on')
            end
        end
    end
    
    if strcmpi(evnt.Character,'v')
        if isfield(S.GD,'VNquiver')
            if strcmpi(S.GD.VNquiver.Visible,'on')
                S.GD.VNquiver.Visible = 'off';
                set(xfigure_This.StatusBox, 'String', 'Vertex normals: off')
            else
                S.GD.VNquiver.Visible = 'on';
                set(xfigure_This.StatusBox, 'String', 'Vertex normals: on')
            end
        end
    end
    
    if strcmpi(evnt.Character,'f')
        if isfield(S.GD,'FNquiver')
            if strcmpi(S.GD.FNquiver.Visible,'on')
                S.GD.FNquiver.Visible = 'off';
                set(xfigure_This.StatusBox, 'String', 'Face normals: off')
            else
                S.GD.FNquiver.Visible = 'on';
                set(xfigure_This.StatusBox, 'String', 'Face normals: on')
            end
        end
    end
    
    if strcmpi(evnt.Character,'e')
        if strcmpi(S.hp.EdgeColor, 'none')
            S.hp.EdgeColor = 'k';
            set(xfigure_This.StatusBox, 'String', 'EdgeColor: k')
        else
            S.hp.EdgeColor = 'none';
            set(xfigure_This.StatusBox, 'String', 'EdgeColor: none')
        end
    end
    
%     if strcmpi(evnt.Key,'numpad4')
%        [az,el] = view();
%        az = az+5;
%        view(az,el)
%        set(xfigure_This.StatusBox, 'String', ['Az: ',num2str(az), ' El: ', num2str(el) ])
%     end
%     if strcmpi(evnt.Key,'numpad6')
%        [az,el] = view();
%        az = az-5;
%        view(az,el)
%        set(xfigure_This.StatusBox, 'String', ['Az: ',num2str(az), ' El: ', num2str(el) ])
%     end
%     if strcmpi(evnt.Key,'numpad8')
%        [az,el] = view();
%        el = el-5;
%        view(az,el)
%        set(xfigure_This.StatusBox, 'String', ['Az: ',num2str(az), ' El: ', num2str(el) ])
%     end
%     if strcmpi(evnt.Key,'numpad2')
%        [az,el] = view();
%        el = el+5;
%        view(az,el)
%        set(xfigure_This.StatusBox, 'String', ['Az: ',num2str(az), ' El: ', num2str(el) ])
%     end
%     if strcmpi(evnt.Key,'numpad5')
%         [az,el] = view();
%         v = [az,el];
%         roundTargets = [-360 -270 -180 -90 0 90 180 270 360];
%         vRounded = interp1(roundTargets,roundTargets,v,'nearest');
%         az = vRounded(1);
%         el = vRounded(2);
%         view([az, el]);
%         set(xfigure_This.StatusBox, 'String', ['Az: ',num2str(az), ' El: ', num2str(el) ])
%     end
    
    xfigure_This.uiTextHelp.String = {'Press "G" to toggle grid';
                                      'Press "CTRL+R" to reset axis';
                                      'Press "Middlemousebuton and drag" to pan axis';
                                      'Press "CTRL+Click and drag" to rotate axis';
                                      'Press "SHIFT+S" snap the view to the view-box';
                                      'Press "SHIFT+P" to save the current window position';
                                      'Press numpad keys to rotate';
                                      'Press "e" to toggle edge off and on';
                                      'Press "l" to toggle lights';
                                      'Press "v" to toggle vertex normals';
                                      'Press "f" to toggle face normals';
                                      'Press "s" to toggle surface shading'};
    
    xfigure_This.uiTextHelp.Position = [5    20   261   176];
end

function [P,T] = weldPoints(X,tri)
    mt = 0;
    map = zeros(size(X,1),2);
    % Renumbering nodes by creating map
    % Map which nodes to replace with other nodes
    % Second column idices replaces first column
    c = 0;
    for iel = 1:size(tri)
        itri = tri(iel,:);
        if iel == 1
            map(1:3,:)=[(1:3)',(1:3)'];
            c = 3;
            mt = 4;
        else
            for indt=1:3
                it = itri(indt);
                if it >= mt && ~any(it==map(:,1))
                    map(c+1,:) = [it,mt];
                    mt=mt+1;
                    c=c+1;
                end
                
            end
            
        end
    end
    map = map(all(map~=0,2),:);
    
    for im = 1:size(map,1)
        tri(tri==map(im,1)) = map(im,2);
    end
    T = tri;
    P =  X(map(:,1),:);
end

function bandInds = computeBoundingBox(X,iP,ShowBox)
    ipnts = X(iP:iP+2,:);
    xc0 = min(ipnts(:,1));
    xc1 = max(ipnts(:,1));
    yc0 = min(ipnts(:,2));
    yc1 = max(ipnts(:,2));
    zc0 = min(ipnts(:,3));
    zc1 = max(ipnts(:,3));
    bandIndsx = X(:,1)>=xc0 & X(:,1)<=xc1;
    bandIndsy = X(:,2)>=yc0 & X(:,2)<=yc1;
    bandIndsz = X(:,3)>=zc0 & X(:,3)<=zc1;
    bandInds = all([bandIndsx,bandIndsy,bandIndsz],2);
    if ShowBox
        drawBBox(X,bandInds,xc1,xc0,yc1,yc0,zc1,zc0,iP)
    end
end

function tri = PrivateTriangulate(X,NF,ShowBox)
% Strategy:
% The overal strategy is to loop over each set of triangle points
% and compare each point to the previous added list of points. If
% the current point exists in the list, then it was previously added
% and its position in the list is added to the list. Otherwise the
% number of nodes is incremented and that number is added to the
% triangulation list.
%   Example:
%   First set is initiated to [1,2,3]
%   In the next iteration a point in the set is already in the list,
%   the point number is determined to be 2, so two new nodes are
%   added:
%   [1,2,3;
%    4,5,3]
%   And so on...
%
% To improve performance a bounding box is sweeping over every triangle
% (every third point of the list). The bounding box ecapsulates a number of
% points which are the search space. So instead of for each point
% compare the distance to every other point, we compute the distance to
% a constant set of other points.
    nnod = size(X,1);
    nele = nnod/3;
    tri = zeros(nele, 3);
    itri = 1;
    for iP = 1:3:size(X,1)
        if itri == 1
            tri(1,1:3) = [1,2,3];
            
            xe = X(tri(itri,:),:);
            n = cross( xe(2,:)-xe(1,:) , xe(3,:)-xe(1,:) );
            n = n/norm(n);
            a = NF(itri,:)*n(:);
            if a < 0
                tri(1,:) = [1,3,2];
            end
            
            itri = itri +1;
        else
            %% Bounding Box
            bandInds = computeBoundingBox(X,iP,ShowBox);
            PrevPoints = [1:(iP+2-1)]';
            ninds = sum(bandInds);
            if ( length(PrevPoints) >= ninds )
                searchSpaceBand = X(bandInds,:);
                SearchIndsBand = find(bandInds); 
            end
            %% Loop over the three next points in the set of triangles
            % the local triangle node index is set to 1
            tind = 1;
            for iPnt = iP:iP+2
                % If the neighborhood point space is smaller than the space of
                % points we have found so far, then define the searchSpace as
                % the neighborhood points. Otherwise define the searchSpace as
                % the points we have found so far.
                % This is done because in the first few iterations the points
                % found so far are fewer than all the neighborhood points, but
                % later the neighborhgood is pretty mush constant where as the
                % points found so far would only increase.

                PrevPoints = [1:iPnt-1]';


                ninds = sum(bandInds);
                if ( length(PrevPoints) > ninds )
                    searchSpace = searchSpaceBand;
                    SearchInds = SearchIndsBand;
                else
                    searchSpace = X(PrevPoints,:);
                    SearchInds = PrevPoints;
                end
                %

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
                indt = SearchInds(NDist<eps*100);


                % if duplicate points exist; set the local triangle index to
                % the index of the searchSpace and increment the triangle
                % index.
                % If no duplicate points are found; set the local triangle
                % index to the index of the current triangle set point
                if ~isempty(indt)
                    indt = indt(1);
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
            
            itri = itri +1;
        end
    end
end

function drawBBox(X,bandInds,xc1,xc0,yc1,yc0,zc1,zc0,iP)
    %% Draw bounding box
    BX = [xc0,yc0,zc0;
        xc1,yc0,zc0;
        xc1,yc1,zc0;
        xc0,yc1,zc0;
        xc0,yc0,zc1;
        xc1,yc0,zc1;
        xc1,yc1,zc1;
        xc0,yc1,zc1];
    faces = [1,2,3,4;
        5,6,7,8;
        1,2,6,5;
        2,3,7,6;
        3,4,8,7;
        4,1,5,8];
    xfigure(23);
    hold off
    plot3(BX(:,1),BX(:,2),BX(:,3),'sk');hold on;
    patch('Faces',faces,'Vertices',BX,'FaceColor','none');
    axis equal;

    ssb = X(bandInds,:);
    plot3(ssb(:,1),ssb(:,2),ssb(:,3),'ok');
    
    PrevPoints = [1:(iP+2)];
    pnts = iP:iP+2;
    ssp = X(PrevPoints,:);
    plot3(ssp(:,1),ssp(:,2),ssp (:,3),'*b')
    plot3(X(pnts,1),X(pnts,2),X (pnts,3),'+r','MarkerSize',12)
%     dx = (xc1-xc0);
%     dy = (yc1-yc0);
%     dz = zc1-zc0;
%     f = 2;
    X0 = min(X(:,1));
    X1 = max(X(:,1));
    Y0 = min(X(:,2));
    Y1 = max(X(:,2));
    Z0 = min(X(:,3));
    Z1 = max(X(:,3));
    axis([X0,X1,Y0,Y1,Z0,Z1])
%     axis([xc0-dx*f,xc1+dx*f,yc0-dy*f,yc1+dy*f,zc0-dz*f,zc1+dz*f])
    legend('','','Local Searchspace','Previous Points','Current Point')
    title([num2str(size(ssb,1)),' points in box. ',num2str(length(PrevPoints)),' previous points'])
% 	commandwindow
%     waitfor(gcf,'CurrentCharacter');
    disp('Press any key or click on the figure to continue')
    w = waitforbuttonpress;
end
