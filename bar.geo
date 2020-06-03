l  = 1000.;
w1 = 1000.;
w  = 1000.;


h = 100.;

n = 10;

Point(1) = {0,0,0,h};
Point(2) = {l,0,0,h};
Point(3) = {l,w1,0,h};
Point(4) = {0,w1,0,h};


Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};


Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

Transfinite Line{1,2,3,4} = n+1;
Transfinite Surface{1} = {1,2,3,4};


Recombine Surface{1};

Physical Point(111) = {4}; // force

//+
Extrude {0, 0, 1000} {
  Surface{1}; Layers{10}; Recombine;
}

Physical Surface("TopSurface") = {128}; // force
Physical Volume("All") = {1}; // force
