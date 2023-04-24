input = "rectangle.dat";

rectangle =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "PlaneStrain";
    E    = 1.e5;
    nu   = 0.25;
  };
};

solver =
{
  type     = "LinearSolver";
};

outputModules = ["vtk"];

vtk =
{
  type = "MeshWriter";
};

graph =
{
  type = "GraphWriter";
  onScreen = true;

  columns = ["disp", "load"];

  disp =
  {
    type = "state";
    node = "left";
    dof  = "u";
  };

  load =
  {
    type = "fint";
    node = "right";
    dof  = "u";
  };
};

output =
{
  type = "OutputWriter";

  onScreen = true;
};
