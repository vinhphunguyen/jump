h=2.;
L0 = 2.5;
r1 = 37;

n = 6;

Point(1) = {0,0, 0, h/4};
Point(2) = {r1,0, 0, h};
Point(3) = {0, 0,r1, h};
Point(4) = {-r1,0, 0, h};
Point(5) = {0, 0, -r1,h};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(5) = {1, 2, 3, 4};

Plane Surface(60) = {5};


Extrude {0, L0,0} {
  Surface{60}; Layers{n}; Recombine;
}

Physical Volume("All") = {1}; // force
Physical Surface("fix") = {77,73,81,69}; // force