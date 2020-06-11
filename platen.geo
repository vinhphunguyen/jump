l  = 200.;
w1 = 10.;
w  = 1000.;


h = 100.;

n = 30;

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

Transfinite Line{1,3} = n+1;
Transfinite Line{2,4} = 2;
Transfinite Surface{1} = {1,2,3,4};


Recombine Surface{1};

Physical Point(111) = {4}; // force

Physical Surface("All") = {1}; // force

