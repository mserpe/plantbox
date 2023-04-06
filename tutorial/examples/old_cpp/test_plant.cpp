#include "../../../src/structural/MappedOrganism.h"
#include <iostream>
#include <string> 

int main(int argc, char **argv) 
{
	std::string path = "../../../modelparameter/structural/plant/";
    std::string name = "Heliantus_Pagès_2013";
    auto rs = std::make_shared<CPlantBox::MappedPlant>(2);//pb.RootSystem()  # the original
	rs->readParameters(path + name + ".xml");
	
	rs->initialize(false);
	rs->simulate(76,false);//db
	std::cout <<path + name + ".xml"<<std::endl;
    return 0;
}