
############################################################################
#  Description: The general PyFEM input file of the example presented in   #
#               section 2.6 of the book, pages 53--62.                     #
#                                                                          #
#  Use:         pyfem PatchTest4.pro                                       #
############################################################################

input = "PatchTest8_3D.dat";

ContElem =
{
  type = "SmallStrainContinuum";

  material =
  {
    type = "Isotropic";
    E    = 1.e6;
    nu   = 0.25;
  };
};

solver =
{
  type = "LinearSolver";
};

outputModules = ["vtk","output"];

vtk =
{
  type = "MeshWriter";
};

output =
{
  type = "OutputWriter";

  onScreen = true;
};
