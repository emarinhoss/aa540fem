// Flow past a cylinder in a channel: the Schaefer-Turek 2D benchmark geometry.
// Channel [0, 2.2] x [0, 0.41], cylinder of radius 0.05 centred at (0.2, 0.2).
// Physical curves: inlet, outlet, walls, cylinder.  Meshed by make_meshes.py.
lc = 0.05;
lcc = 0.006;

Point(1) = {0, 0, 0, lc};
Point(2) = {2.2, 0, 0, lc};
Point(3) = {2.2, 0.41, 0, lc};
Point(4) = {0, 0.41, 0, lc};
Point(5) = {0.2, 0.2, 0, lcc};
Point(6) = {0.25, 0.2, 0, lcc};
Point(7) = {0.2, 0.25, 0, lcc};
Point(8) = {0.15, 0.2, 0, lcc};
Point(9) = {0.2, 0.15, 0, lcc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};

Physical Curve("inlet") = {4};
Physical Curve("outlet") = {2};
Physical Curve("walls") = {1, 3};
Physical Curve("cylinder") = {5, 6, 7, 8};
Physical Surface("fluid") = {1};

// Refine towards the cylinder and its wake
Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling = 200;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lcc;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.03;
Field[2].DistMax = 0.5;
Field[3] = Box;
Field[3].VIn = 0.015;
Field[3].VOut = lc;
Field[3].XMin = 0.1;
Field[3].XMax = 0.9;
Field[3].YMin = 0.1;
Field[3].YMax = 0.3;
Field[4] = Min;
Field[4].FieldsList = {2, 3};
Background Field = 4;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
