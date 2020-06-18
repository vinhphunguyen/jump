h=40.;
r1 = 0.5e3;

Point(1) = {0,0, 0, h};
Point(2) = {r1,0, 0, h};
Point(3) = {0, r1, 0, h};
Point(4) = {-r1,0, 0, h};
Point(5) = {0, -r1, 0, h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(5) = {1, 2, 3, 4};

Plane Surface(60) = {5};
Physical Surface("All")={60};
Physical Line("boundary")={1,2,3,4};

