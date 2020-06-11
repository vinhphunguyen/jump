h=0.3;
D0 = 7.6;
L0 = 25.4;
r1 = D0/2;

n = 40;

Point(1) = {0,0, 0, h};
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
  Surface{60}; Layers{n};
}

Physical Volume("All") = {1}; // force