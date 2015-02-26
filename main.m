clear
clc
close all

S = STLSurfaceRead('STLSurfaceRead/Examples/teapot.stl')
VizSurface(S,'FaceNormals','VertexNormals')

%% Weld Vertecies
disp('Welding vertices...')
tic
S.WeldPoints();
toc
%% Viz surface
VizSurface(S,'FaceNormals','VertexNormals')