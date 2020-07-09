h=1.2;
h1=0.9;
h2=0.1;
r1 = 8.0;
n = 1;

Point(1) = {0,0, 0, h};
Point(2) = {r1,0, 0, h1};
Point(3) = {0, -r1, 0, h2};
Point(4) = {-r1,0, 0, h1};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Line(3)   = {2,4};
Line Loop(5) = {-1, -2,3};

Plane Surface(60) = {5};


//+
Extrude {0, 0, 1} {
  Surface{60};  Layers{n}; Recombine;
}

Physical Surface("boundary")={76,68};
Physical Surface("force")={72};
Physical Volume("All")={1};