// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "leafparameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <assert.h>

namespace CPlantBox {

/**
 * @return Mean maximal leaf length of this leaf type
 */
double LeafSpecificParameter::getK() const {
	double l = std::accumulate(ln.begin(), ln.end(), 0.);
	return l+la+lb;
}

/**
 * @copydoc OrganParameter::toString()
 */
std::string LeafSpecificParameter::toString() const
{
	std::stringstream str;
	str << "subType\t" << subType << std::endl;
	str << "lb\t" << lb << std::endl << "la\t" << la << std::endl;
	str << "r\t" << r << std::endl << "a\t" << a << std::endl;
	str << "theta\t" << theta << std::endl << "rlt\t" << rlt << std::endl;
	str << "ln\t";
	for (int i=0; i<ln.size(); i++) {
		str << ln[i] << " ";
	}
	str << std::endl;
	return str.str();
}



/**
 * Default constructor sets up hashmaps for class introspection
 */
LeafRandomParameter::LeafRandomParameter(std::shared_ptr<Organism> plant) :OrganRandomParameter(plant)
{
	// base class default values
	name = "undefined";
	organType = Organism::ot_leaf;
	subType = -1;
	f_tf = std::make_shared<Tropism>(plant);
	bindParameters();
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
std::shared_ptr<OrganRandomParameter> LeafRandomParameter::copy(std::shared_ptr<Organism> plant)
{
    // std::cout << "LeafRandomParameter::copy\n"<< std::flush;
	auto r = std::make_shared<LeafRandomParameter>(*this); // copy constructor breaks class introspection
	r->plant = plant;
	r->bindParameters(); // fix class introspection
	r->f_tf = f_tf->copy(plant); // copy call back classes
	r->f_gf = f_gf->copy();
	r->f_se = f_se->copy();
	r->f_sa = f_sa->copy();
	r->f_sbp = f_sbp->copy();
	return r;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific leaf from the leaf type parameters.
 * @return Specific leaf parameters derived from the leaf type parameters
 */
std::shared_ptr<OrganSpecificParameter> LeafRandomParameter::realize()
{
    auto p = plant.lock();
	// type does not change
	double lb_ = std::max(lb + p->randn()*lbs,double(0)); // length of basal zone
	double la_ = std::max(la + p->randn()*las,double(0)); // length of apical zone
	std::vector<double> ln_; // stores the inter-distances

	// stores the inter-distances
	int nob_ = std::max(round(nob() + p->randn()*nobs()),double(0)); // maximal number of leafs
	switch(lnf) {
	case 0: // homogeneously distributed stem nodes
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d = std::max(ln + p->randn()*lns,1e-9); //Normal function of equal internode distance
			ln_.push_back(d);
		}
		break;
	case 1: //nodes distance increase linearly TODO
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln*(1+i) + p->randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
			ln_.push_back(0);
		}
		break;
	case 2: //nodes distance decrease linearly TODO
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d =  std::max(ln*(1+i) + p->randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
		}
		break;
	case 3: //nodes distance increase exponential TODO
		for (int i = 0; i<nob_-1; i++) { // create inter-stem distances
			double d =  std::max(ln + p->randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
		}
		break;
	case 4://nodes distance decrease exponential TODO
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln/(1+i) + p->randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
			ln_.push_back(0);
		}
		break;
	case 5://nodes distance decrease exponential
		for (int i = 0; i<nob_*2-1; i++) { // create inter-stem distances
			double d =  std::max(ln/(1+i) + p->randn()*lns,1e-9); //std::max(  );//ln + randn()*lns,1e-9);
			ln_.push_back(d);
		}
		break;
	default:
		throw 1; // TODO make a nice one
	}
	double r_ = std::max(r + p->randn()*rs,double(0)); // initial elongation
	double a_ = std::max(a + p->randn()*as,double(0)); // radius
	double theta_ = std::max(theta + p->randn()*thetas,double(0)); // initial elongation
	double rlt_ = std::max(rlt + p->randn()*rlts,double(0)); // leaf life time
	return std::make_shared<LeafSpecificParameter>(subType,lb_,la_,ln_,r_,a_,theta_,rlt_);
}

/**
 * Choose (dice) lateral type based on leaf parameters successor and successorP
 *
 * @param pos       spatial position (for coupling to a soil model)
 * @return          leaf sub type of the lateral leaf
 */
int LeafRandomParameter::getLateralType(const Vector3d& pos)
{
	assert(successor.size()==successorP.size()
			&& "LeafTypeParameter::getLateralType: Successor sub type and probability vector does not have the same size");
	if (successorP.size()>0) { // at least 1 successor type
		double d = plant.lock()->rand(); // in [0,1]
		int i=0;
		double p=successorP.at(i);
		i++;
		while ((p<d) && (i<successorP.size())) {
			p+=successorP.at(i);
			i++;
		}
		if (p>=d) { // success
			// std::cout << "lateral type " << successor.at(i-1) << "\n" << std::flush;
			return successor.at(i-1);
		} else { // no successors
			// std::cout << "no lateral type " << std::flush;
			return -1;
		}
	} else {
		return -1; // no successors
	}
}

/**
 * todo docme
 *
 * todo I have no idea why this holds...
 */
double LeafRandomParameter::nobs() const
{
    double nobs = (lmaxs/lmax - lns/ln)*lmax/ln; // error propagation
    if (la>0) {
        nobs -= (las/la - lns/ln)*la/ln;
    }
    if (lb>0) {
        nobs -= (lbs/lb - lns/ln)*lb/ln;
    }
    return std::max(nobs,0.);
}

/**
 * @copydoc OrganTypeParameter::toString()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
std::string LeafRandomParameter::toString(bool verbose) const {

	if (verbose) {
		std::string s = OrganRandomParameter::toString(true);
		std::stringstream str;
		str << "successor\t";
		for (int i=0; i<successor.size(); i++) {
			str << successor[i] << " ";
		}
		str << "\t" << description.at("successor") << std::endl;
		str << "successorP\t";
		for (int i=0; i<successorP.size(); i++) {
			str << successorP[i] << " ";
		}
		str << "\t" << description.at("successorP") << std::endl;
		return s.insert(s.length()-4, str.str());
	} else {
		return OrganRandomParameter::toString(false);
	}

}

/**
 * @copydoc OrganTypeParameter::readXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 *
 * If the parameter successor or successorP are not in the element, they are set to zero size.
 */
void LeafRandomParameter::readXML(tinyxml2::XMLElement* element)
{
	OrganRandomParameter::readXML(element);
	tinyxml2::XMLElement* p = element->FirstChildElement("parameter");
	successor.resize(0);
	successorP.resize(0);
	while(p) {
		std::string key = p->Attribute("name");
		if (key.compare("successor")==0)  {
			successor.push_back(p->IntAttribute("type"));
			successorP.push_back(p->DoubleAttribute("percentage"));
		}
		p = p->NextSiblingElement("parameter");
	}
	double p_ = std::accumulate(successorP.begin(), successorP.end(), 0.);
	if  ((p_<1) && (p_!=0))  {
		std::cout << "LeafRandomParameter::readXML: Warning! percentages to not add up to 1. \n";
	}
	assert(successor.size()==successorP.size() &&
			"LeafTypeParameter::readXML: Successor sub type and probability vector does not have the same size" );
}

/**
 * @copydoc OrganTypeParameter::writeXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
tinyxml2::XMLElement* LeafRandomParameter::writeXML(tinyxml2::XMLDocument& doc, bool comments) const
{
	assert(successor.size()==successorP.size() &&
			"LeafTypeParameter::writeXML: Successor sub type and probability vector does not have the same size" );
	tinyxml2::XMLElement* element = OrganRandomParameter::writeXML(doc, comments);
	for (int i = 0; i<successor.size(); i++) {
		tinyxml2::XMLElement* p = doc.NewElement("parameter");
		p->SetAttribute("name", "successor");
		p->SetAttribute("number", i);
		p->SetAttribute("type", successor[i]);
		p->SetAttribute("percentage", float(successorP[i]));
		element->InsertEndChild(p);
		if (comments) {
			std::string str = description.at("successor");
			tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
			element->InsertEndChild(c);
		}

	}
	double p_ = std::accumulate(successorP.begin(), successorP.end(), 0.);
	if ((p_<1) && (p_!=0)) {
		std::cout << "LeafRandomParameter::writeXML: Warning! percentages do not add up to 1. = " << p_ << "\n";
	}
	return element;
}

/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void LeafRandomParameter::bindParameters()
{
    OrganRandomParameter::bindParameters();
	bindParameter("lb", &lb, "Basal zone [cm]", &lbs);
	bindParameter("la", &la, "Apical zone [cm]", &las);
	bindParameter("ln", &ln, "Inter-lateral distance [cm]", &lns);
    bindParameter("lmax", &lmax, "Maximal stem length [cm]", &lmaxs);
	bindParameter("r", &r, "Initial growth rate [cm day-1]", &rs);
	bindParameter("a", &a, "Leaf width [cm]", &as);
    bindParameter("RotBeta", &rotBeta, "RevRotation of the leaf"); /// todo improve description, start lower letter
    bindParameter("BetaDev", &betaDev, "RevRotation deviation"); /// todo improve description, start lower letter
    bindParameter("InitBeta", &initBeta, "Initial RevRotation"); /// todo improve description, start lower letter
	bindParameter("tropismT", &tropismT, "Type of leaf tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
	bindParameter("tropismN", &tropismN, "Number of trials of leaf tropism");
	bindParameter("tropismS", &tropismS, "Mean value of expected change of leaf tropism [1/cm]");
	bindParameter("dx", &dx, "Axial resolution [cm] (maximal segment size)");
    bindParameter("dxMin", &dxMin, "Axial resolution [cm] (minimal segment size)");
	bindParameter("theta", &theta, "Angle between leaf and parent leaf [rad]", &thetas);
	bindParameter("rlt", &rlt, "Leaf life time [day]", &rlts);
	bindParameter("gf", &gf, "Growth function number [1]", &rlts);
	bindParameter("lnf", &lnf, "Type of inter-branching distance (0 homogeneous, 1 linear inc, 2 linear dec, 3 exp inc, 4 exp dec)");
	// other parameters (descriptions only)
	description["successor"] = "Sub type of lateral leaf veins";
	description["successorP"] = "Probability of each sub type to occur";
}

} // end namespace CPlantBox
