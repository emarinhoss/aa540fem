// Annulus r1 <= r <= r2, used as a curved-boundary verification case:
// with T = 0 on the inner circle and T = 1 on the outer circle and
// isotropic conductivity the exact solution is T = ln(r/r1) / ln(r2/r1).
//
// Meshed by make_meshes.py (linear/quadratic, triangles/quads); to mesh
// by hand:  gmsh -2 annulus.geo -o annulus.msh
r1 = 1.0;
r2 = 2.0;
lc = 0.2;

Point(1) = {0, 0, 0, lc};
Point(2) = { r1,   0, 0, lc};
Point(3) = {  0,  r1, 0, lc};
Point(4) = {-r1,   0, 0, lc};
Point(5) = {  0, -r1, 0, lc};
Point(6) = { r2,   0, 0, lc};
Point(7) = {  0,  r2, 0, lc};
Point(8) = {-r2,   0, 0, lc};
Point(9) = {  0, -r2, 0, lc};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Curve Loop(1) = {5, 6, 7, 8};
Curve Loop(2) = {1, 2, 3, 4};
Plane Surface(1) = {1, 2};

Physical Curve("inner") = {1, 2, 3, 4};
Physical Curve("outer") = {5, 6, 7, 8};
Physical Surface("domain") = {1};
