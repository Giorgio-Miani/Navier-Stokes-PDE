H = 0.41;   // Square side
L = 2.5;    // Prism lenght
R = 0.05;   // Cylinder radius
N = 2;      // Cylinder Layers
a = 0.5;    // Cylinder center to Inflow Plane
b = 0.2;    // Cylinder center to Bottom

lc1 = 1; // Mesh size 1
lc2 = 1; // Mesh size 2


// Points. //////////////////////////////////////////////////////////////////////

// Inflow Plane //

Point(1) = {0, 0, 0, lc1};
Point(2) = {0, 0, H, lc1};
Point(3) = {0, H, H, lc1};
Point(4) = {0, H, 0, lc1};

// Outflow Plane //

Point(5) = {L, 0, 0, lc1};
Point(6) = {L, 0, H, lc1};
Point(7) = {L, H, H, lc1};
Point(8) = {L, H, 0, lc1};

// Circles //

Point(9)  = {a, b, 0, lc1};

Point(10) = {a + R, b, 0, lc2};
Point(11) = {a, b + R, 0, lc2};
Point(12) = {a - R, b, 0, lc2};
Point(13) = {a, b - R, 0, lc2};

Point(14)  = {a, b, H, lc1};

Point(15) = {a + R, b, H, lc2};
Point(16) = {a, b + R, H, lc2};
Point(17) = {a - R, b, H, lc2};
Point(18) = {a, b - R, H, lc2};


// Lines ////////////////////////////////////////////////////////////////////////

// Squares //

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Sides //

Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

// Circle //

Circle(13) = {10, 9, 11};
Circle(14) = {11, 9, 12};
Circle(15) = {12, 9, 13};
Circle(16) = {13, 9, 10};

Circle(17) = {15, 14, 16};
Circle(18) = {16, 14, 17};
Circle(19) = {17, 14, 18};
Circle(20) = {18, 14, 15};


// Squares Loops //

Line Loop(1) = {1, 2, 3, 4};        // Inflow Plane
Line Loop(2) = {5, 6, 7, 8};        // Outflow Plane

// Sides Loops //

Line Loop(3) = {1, 10, -5, -9};     // Bottom
Line Loop(4) = {2, 11, -6, -10};    // Right
Line Loop(5) = {3, 12, -7, -11};    // Top
Line Loop(6) = {4, 9, -8, -12};     // Left

// Circle Loops //

Line Loop(7) = {13, 14, 15, 16};
Line Loop(8) = {17, 18, 19, 20};


// Surfaces /////////////////////////////////////////////////////////////////////

Plane Surface(1) = {1};             // Inflow Plane
Plane Surface(2) = {2};             // Outflow Plane

Plane Surface(3) = {3};             // Bottom
Plane Surface(4) = {4, 8};          // Right
Plane Surface(5) = {5};             // Top
Plane Surface(6) = {6, 7};          // Left

Surface(7) = {7};                   // Circle

Extrude {0, 0, H} {                 //Cylinder
    Surface{7}; Layers{N};
}

Delete{ Volume{1}; Surface{7}; Surface{42}; Point{9}; Point{14};}

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Surface Loop(2) = {29:41:4};


// Volumes //////////////////////////////////////////////////////////////////////

Volume(1) = {1, 2};


// Set tags to the boundaries. //////////////////////////////////////////////////

Physical Surface(0) = {1};                              // Inflow Plane
Physical Surface(1) = {2};                              // Outflow Plane
Physical Surface(2) = {3, 4, 5, 6, 29, 33, 37, 41};     // Sides and Cylinder

Physical Volume(10) = {1};


// Generate mesh ////////////////////////////////////////////////////////////////

Mesh 3;
Save "mesh.msh";