//bar
ri = 80;
ro = 150;
rm = 0.5*(ri+ro);
r1  = rm;
r2  = ro;

h = 1.0;//for ref
n1 = 20; n2 =8;


Point(1) = {0, 0, 0, h};
Point(2) = {r1, 0, 0, h}; Point(3) = {-r1, 0, 0, h};
Point(4) = {r2, 0, 0, h}; Point(5) = {-r2, 0, 0, h};


Circle(1) = {2,1,3}; Circle(2) = {3,1,2};
Circle(3) = {4,1,5}; Circle(4) = {5,1,4};

Line(5) = {2,4}; 
Line(6) = {3,5};


Line Loop(1) = {-1,5,3,-6}; Plane Surface(1) = {1};
Line Loop(2) = {2,5,-4,-6}; Plane Surface(2) = {2}; 

// transfinite 
Transfinite Line {1,2,3,4} = n1 Using Progression 1;
Transfinite Line {5,6} = n2 Using Progression 1;
Transfinite Surface {1,2};


//for Q4/Q8 elements
Mesh.RandomFactor = 1e-9;//default = 1e-9
Recombine Surface{1,2};//T3->T4
//Mesh.ElementOrder = 2;//T4->Q9 //Mesh.SecondOrderLinear = 0;
//Mesh.SecondOrderIncomplete=1;//Q9->Q8


//+
Extrude {0, 0, 10} {
  Surface{1}; Surface{2}; Layers{1}; Recombine;
}


Physical Surface("fix")={19,27};
Physical Volume("All")={1};
