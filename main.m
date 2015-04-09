clear
clc
close all

S = STLSurfaceRead('STLSurfaceRead/Examples/teapot.stl')
VizSurface(S,'FaceNormals','VertexNormals')
title('Stl surface without welded vertices')

%% Weld Vertecies
disp('Welding vertices...')
tic
S.WeldPoints();
toc
%% Viz surface
S2 = S.copy;
VizSurface(S2,'FaceNormals','VertexNormals')
title('Stl surface with welded vertices')



