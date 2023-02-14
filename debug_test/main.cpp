
#include <string>
#include <iostream>

#include "MappedOrganism.h"

const double PI = 3.14159265358979323846;

int main()
{

  auto plant = std::make_shared<CPlantBox::MappedPlant>();
  const int leaf_resolution = 20;
  const int time = 28;
  plant->readParameters("C:/work/CPlantBox/tutorial/examples/python/results/P0_plant.xml");
  for (auto p : plant->getOrganRandomParameter(CPlantBox::Organism::ot_leaf))
  {
    auto leaf = std::dynamic_pointer_cast<CPlantBox::LeafRandomParameter>(p);
    leaf->lb = 0.0;
    leaf->la = 38.41053981;
    leaf->lmax = 38.41053981;
    leaf->areaMax = 54.45388021;
    leaf->leafGeometryPhi = {-PI/2.0, -80.0/180.0*PI, -PI/4.0, 0.0, PI/4.0, PI/2.0};
    leaf->leafGeometryX = { 38.41053981, 1, 1, 0.3, 1, 38.41053981 };
    leaf->tropismT = 1;
    leaf->tropismS = 0.05;
    leaf->tropismAge = 5;
    leaf->createLeafRadialGeometry(leaf->leafGeometryPhi, leaf->leafGeometryX, leaf_resolution);
  }

  for(auto p : plant->getOrganRandomParameter(CPlantBox::Organism::ot_stem))
  {
    auto stem = std::dynamic_pointer_cast<CPlantBox::StemRandomParameter>(p);
    auto r = 0.758517633;
    stem->r = r;
    stem->lmax = static_cast<double>(time - 7) * r;
  }
  for (auto p : plant->getOrganRandomParameter(CPlantBox::Organism::ot_root))
  {
    auto root = std::dynamic_pointer_cast<CPlantBox::RootRandomParameter>(p);
    auto r = 0.758517633;
    root->r = r;
    root->lmax = static_cast<double>(time - 7) * r;
  }

  plant->initialize();
  plant->simulate(time, true);
  std::cout << "Computing Geometry" << std::endl;
  plant->ComputeGeometry();
  auto points = plant->GetGeometry();

  std::cout << "Generated points are "  << points.size() << std::endl;
}
