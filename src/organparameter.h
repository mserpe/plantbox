// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANPARAMETER_H_
#define ORGANPARAMETER_H_

#include <string>
#include <map>
#include <memory>
#include <iostream>
#include <vector>
#include <string>

#include "external/tinyxml2/tinyxml2.h"

/**
 * This file describes the classes OrganSpecificParameter and OrganRandomParameter.
 * OrganSpecificParameter are drawn from the OrganRandomParameter class
 */

namespace CPlantBox {

class Organism; // forward declaration
class GrowthFunction;
class ExponentialGrowth;

/**
 * Parameters for a specific organ
 */
class OrganSpecificParameter {
public:

    OrganSpecificParameter(int t, double a, double lb, double la): subType(t), lb(lb), la(la) , a(a) { }

    virtual ~OrganSpecificParameter() { }

    int subType = -1; ///< sub type of the organ
    double lb = 0.; 		///< Basal zone of leaf (leaf-stem) [cm]
	double la = 0.;			///< Apical zone of leaf vein [cm];
	double a = 0.; ///< radius of the organ [cm]
    virtual std::string toString() const; ///< quick info for debugging

};

/**
 * Contains a parameter set describing as single sub type of an organ,
 * specific parameters are then created with realize().
 *
 * Organizes parameters in hash-maps for scalar double and scalar int values.
 * For this reason derived classes getParameter(), toString(), readXML(), and writeXML() should work out of the box.
 * For other parameter types the methods must be overwritten, see e.g. RootRandomParameter.
 *
 * The factory function copy() has to be overwritten for each specialization.
 */
class OrganRandomParameter
{
public:

    OrganRandomParameter(std::shared_ptr<Organism> plant); ///< default constructor
    virtual ~OrganRandomParameter() { };

    virtual std::shared_ptr<OrganRandomParameter> copy(std::shared_ptr<Organism> plant); ///< copies the root type parameter into a new plant

    virtual std::shared_ptr<OrganSpecificParameter> realize(); ///< creates a specific organ from the root parameter set

    virtual double getParameter(std::string name) const; // get a scalar parameter

    virtual std::string toString(bool verbose = true) const; ///< info for debugging

    virtual void readXML(tinyxml2::XMLElement* element); ///< reads a single sub type organ parameter set
    void readXML(std::string name); ///< reads a single sub type organ parameter set
    virtual tinyxml2::XMLElement* writeXML(tinyxml2::XMLDocument& doc, bool comments = true) const; ///< writes a organ root parameter set
    void writeXML(std::string name) const; ///< writes a organ root parameter set

    virtual void bindParameters(); ///<sets up class introspection

    void bindParameter(std::string name, int* i, std::string descr = "", double* dev = nullptr); ///< binds integer to parameter name
    void bindParameter(std::string name, double* d, std::string descr = "", double* dev = nullptr); ///< binds double to parameter name

    std::string name = "organ";
    int organType = 0;
    int subType = 0;
    double tropismN = 1.;	///< Leaf tropism parameter (number of trials)
	double tropismS = 0.2;	///< Leaf tropism parameter (mean value of expected changeg) [1/cm]
	double lb = 0.; 	///< Basal zone [cm]
	double lbs = 0.;  	///< Standard deviation basal zone [cm]
	double la = 10.;	///< Apical zone [cm];
	double las = 0.;	///< Standard deviation apical zone [cm];
	double a = 0.1; 		///< Root radius [cm]
    double as = 0.; 		///< Standard deviation root radius [cm]
    double dx = 0.25; 		///< Maximal segment size [cm]
	double dxMin = 1e-6; 	///< threshold value, smaller segments will be skipped (otherwise stem tip direction can become NaN)

    std::weak_ptr<Organism> plant;
    std::shared_ptr<GrowthFunction> f_gf;
	std::vector<double> string2vector(std::string xmlInput);///<convert string to vector<double>, to simplifiy xml input

protected:

    /* class introspection */
    std::map<std::string, double*> dparam; ///< Parameters with type double that can be read and written
    std::map<std::string, int*> iparam; ///< Parameters with type double that can be read and written
    std::map<std::string, double*> param_sd; ///< Deviations of parameters
    std::map<std::string, std::string> description; ///< Parameter descriptions
};

} // namespace

#endif
