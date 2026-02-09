// Simple 10x10 uniform quad mesh for freestream test
lc = 0.1;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Force structured quad mesh
Transfinite Curve {1, 3} = 11 Using Progression 1;
Transfinite Curve {2, 4} = 11 Using Progression 1;
Transfinite Surface {1};
Recombine Surface {1};

// Physical groups
Physical Curve("farfield", 1) = {1, 2, 3, 4};
Physical Surface("domain", 2) = {1};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 8;  // Frontal-Delaunay for quads
