// Nicolo Grilli
// University of Oxford
// AWE project 2020

// Square geometry, 20 um edge
// 0.25 um element size

Point(1) = {0, 0, 0, 0.25};
Point(2) = {20, 0, 0, 0.25};
Point(3) = {20, 20, 0, 0.25};
Point(4) = {0, 20, 0, 0.25};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Line loops

Line Loop(11) = {1, 2, 3, 4};

Plane Surface(21) = {11};

// Transfinite Line

Transfinite Line {1, 2, 3, 4} = 81 Using Progression 1;

// Recombine Surface

Recombine Surface {21};

Extrude {0, 0, 0.25} {
  Surface{21};
  Layers{1};
  Recombine;
}

Physical Volume("Part1") = {1};

Physical Surface("Y0") = {30};
Physical Surface("Ymax") = {38};
Physical Surface("X0") = {42};
Physical Surface("Xmax") = {34};
Physical Surface("Z0") = {21};
Physical Surface("Zmax") = {43};


