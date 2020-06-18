l  = 1625;
w1 = 100;
w  = 100;


h = 5.;

m = 100;
n = 2;

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

Transfinite Line{1,3} = m+1;
Transfinite Line{2,4} = n+1;
Transfinite Surface{1} = {1,2,3,4};


Recombine Surface{1};


Physical Line("boundary") = {1,2,3,4}; // force
Physical Line("BotSurface") = {1}; // force
Physical Surface("All") = {1}; // force
