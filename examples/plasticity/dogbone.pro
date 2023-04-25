input = "dogbone.dat";

Continuum =
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

outputModules = ["vtk" , "GraphWriter"];

vtk =
{
  type = "MeshWriter";
};

GraphWriter =
{
  onScreen = true;

  columns = [ "disp" , "load" ];

  disp =
  {
    type = "state";
    node = 17;
    dof  = 'u';
  };
 
  load =
  {
    type = "fint";
    node = load_nodes1;
    dof  = 'u';
  };
};


