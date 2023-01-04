Mesh.MshFileVersion = 2.2;

// Geometric parameters
H = 10;  //  semiheight of plate
L = H/2;  //  semiwidth of plate
R = L/4;

// Discretization parameters
lc1 = 1; // element size at the border
lc2 = 1; // element size at the crack tip

// Domain construction
Point(1) = {L,0.0,0.0,lc1};
Point(2) = {L,H,0.0,lc1};
Point(3) = {0.0,H,0.0,lc1};
Point(4) = {0.0,0.0,0.0,lc1};
Point(5) = {(L/2)+R,H/2,0.0,lc1};
Point(6) = {L/2,H/2,0.0,lc1};
Point(7) = {(L/2)-R,H/2,0.0,lc1};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Circle(5) = {7,6,5};
Circle(6) = {5,6,7};



Line Loop(11) = {1,2,3,4};
Line Loop(12) = {5,6};

Plane Surface(1) = {11,12};
Plane Surface(2) = {12};



// To return only 2D element in msh file
Physical Surface(1) = {1};
Physical Surface(2) = {2};
