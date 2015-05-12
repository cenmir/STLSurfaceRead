clear
clc
close all

S = STLSurfaceRead('STLSurfaceRead/Examples/teapot.stl')
VizSurface(S,'FaceNormals','VertexNormals')
title('Stl surface without welded vertices')

%% Triangulate surface
disp('Triangulating vertices...')
tic
S.TriangulateSurface();
toc
%% Viz surface
% S2 = S.copy;
VizSurface(S,'FaceNormals','VertexNormals')
title('Stl surface with welded vertices')



%%