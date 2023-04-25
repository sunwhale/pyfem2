input = "rectangle.dat";

rectangle =
{
  type = "SmallStrainContinuum";

  material =
  {
    type   = "IsotropicKinematicHardening";
    E      = 1.e6;
    nu     = 0.25;
    syield = 10000.;
    hard   = 1.e4;
  };
};

solver =
{
  type = "NonlinearSolver";

  maxCycle = 20;
};

outputModules = ["vtk", "GraphWriter"];

vtk =
{
  type = "MeshWriter";
};

GraphWriter =
{
  type = "GraphWriter";
  onScreen = true;

  columns = ["disp", "load"];

  disp =
  {
    type = "state";
    node = "right";
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
