Mesh.MshFileVersion = 2.2;

// Geometric parameters
H = 10;  //  Height of the Bi-Material Beam
L = H/2;  //  Width of the Bi-Material Beam

// Discretization Parameter
lc1 = 1; // Element Size on the Edges

// Domain construction

Point(1) = {L,0.0,0.0,lc1};
Point(2) = {L,H,0.0,lc1};
Point(3) = {0.0,H,0.0,lc1};
Point(4) = {0.0,0.0,0.0,lc1}; 

Point(5) = {0.0,H/2,0.0,lc1};
Point(6) = {L,H/2,0.0,lc1};

Line(1) = {1,6};
Line(2) = {6,5};
Line(3) = {5,4};
Line(4) = {4,1};

Line(5) = {2,6};
Line(7) = {5,3};
Line(6) = {3,2};

Line Loop(11) = {1,2,3,4};
Plane Surface(1) = {11};

Line Loop(12) = {5,2,7,6};
Plane Surface(2) = {12};

// Return 2D Element in .msh File
Physical Surface(1) = {1};
Physical Surface(2) = {2};
