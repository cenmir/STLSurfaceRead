# STLSurfaceRead
![surface.png](https://raw.githubusercontent.com/cenmir/STLSurfaceRead/master/STLSurfaceRead/Examples/surface.png)
## Contents
**STLSurfaceRead class**
Constructor.
	S = STLSurfaceRead(filname)

**VizSurface**
Visualize the surface

	VizSurface(S,'FaceNormals','VertexNormals')

**WeldPoints**
Removes duplicate points and creates surface triangulation.

	WeldPoints(S)

or

	S.WeldPoints


## Example
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