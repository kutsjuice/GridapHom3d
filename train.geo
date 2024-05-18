// Gmsh project created on Sat May 11 17:01:24 2024
SetFactory("OpenCASCADE");

//+
Rectangle(8) = {-1., -1, -0.5, 2, 2, 0};
//+
Rectangle(9) = {-1, -1, 0.5, 2, 2, 0};

//+
Sphere(2) = {-0.5, -0.5, -0.5, 0.1, -Pi/2, Pi/2, 2*Pi};
//+
Rectangle(11) = {-1, -1, 0.5, 2, 2, 0};
//+
Rectangle(12) = {-1, -1, 0.5, 2, 2, 0};
//+
Rectangle(13) = {-1, -1, 0.5, 2, 2, 0};
//+
Rectangle(14) = {-1, -1, 0.5, 2, 2, 0};
