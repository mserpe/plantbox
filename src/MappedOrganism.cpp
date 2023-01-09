// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "MappedOrganism.h"

#include "SegmentAnalyser.h"
#include "growth.h"
#include <algorithm>
#include <functional>
#include <cmath>

namespace CPlantBox {

  /** 
   * Hidden Template method to write a set of vectors into a buffer
   * @param buffer    the buffer
   * @param offset    the offset in the buffer
   * @param v...     the vectors
   */
  inline unsigned int vec2Buf(std::vector<double>& buffer, unsigned int offset, const Vector3d& v) {
    buffer[offset] = v.x;
    buffer[offset + 1] = v.y;
    buffer[offset + 2] = v.z;
    return offset + 3;
  }
  template<typename... Args>
  inline unsigned int vec2Buf(std::vector<double>& buffer, unsigned int offset, const Vector3d& v, Args... args) {
    offset = vec2Buf(buffer, offset, v);
    return vec2Buf(buffer, offset, args...);
  }

/**
 * A static plant, as needed for flux computations, represented as
 *
 * @param nodes     	coordinates [cm]
 * @param nodeCTs   	node creation times [d]
 * @param segs      	describes a segment with two node indices segs.x, and segs.y [1]
 * @param radii     	segment radius [cm]
 * @param subTypes     	root type or order of the segment [1]
 * @param organTypes    organ type [1]
 */
MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
		std::vector<double> radii, std::vector<int> subTypes, std::vector<int> organTypes) :
									nodes(nodes), nodeCTs(nodeCTs), segments(segs), radii(radii), subTypes(subTypes), organTypes(organTypes)
{
	assert((nodes.size()==nodeCTs.size()) && "MappedSegments::MappedSegments: Unequal vector sizes nodes and nodeCTs");
	assert((segments.size()==radii.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and radii");
	assert((segments.size()==subTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and subTypes");
	assert((segments.size()==organTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and organ types");
}

/**
 * A static root system, as needed for flux computations, represented as
 *
 * @param nodes     	coordinates [cm]
 * @param nodeCTs   	node creation times [d]
 * @param segs      	describes a segment with two node indices segs.x, and segs.y [1]
 * @param radii     	segment radius [cm]
 * @param subTypes     	root type or order of the segment [1]
 */
MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<double> nodeCTs, std::vector<Vector2i> segs,
		std::vector<double> radii, std::vector<int> subTypes) :
									nodes(nodes), nodeCTs(nodeCTs), segments(segs), radii(radii), subTypes(subTypes)
{
	organTypes.resize(segments.size());
	std::fill(organTypes.begin(), organTypes.end(), Organism::ot_root);
	assert((nodes.size()==nodeCTs.size()) && "MappedSegments::MappedSegments: Unequal vector sizes nodes and nodeCTs");
	assert((segments.size()==radii.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and radii");
	assert((segments.size()==subTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and subTypes");
	assert((segments.size()==organTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and organ types");
}
/**
 *  A static root system, as needed for flux computations.
 *
 * @param nodes     [cm]
 * @param segs      describes a segment with two node indices segs.x, and segs.y [1]
 * @param radii     [cm] segment radius
 *
 * nodeCTs is set to 0. for all segments
 * subTypes are set to type 0 for all segments
 */

MappedSegments::MappedSegments(std::vector<Vector3d> nodes, std::vector<Vector2i> segs, std::vector<double> radii)
:nodes(nodes), segments(segs), radii(radii) {
	nodeCTs.resize(nodes.size());
	std::fill(nodeCTs.begin(), nodeCTs.end(), 0.);
	organTypes.resize(segments.size());
	std::fill(organTypes.begin(), organTypes.end(), Organism::ot_root);
	setSubTypes(0);
	assert((segments.size()==radii.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and radii");
	assert((segments.size()==subTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and subTypes");
	assert((segments.size()==organTypes.size()) && "MappedSegments::MappedSegments: Unequal vector sizes segments and organTypes");
}

/**
 * Sets the radius @param a [cm] for all segments.
 */
void MappedSegments::setRadius(double a) {
	radii.resize(segments.size());
	std::fill(radii.begin(), radii.end(), a);
}

/**
 * Sets the sub type of all segments to @param t.
 */
void MappedSegments::setSubTypes(int t) {
	subTypes.resize(segments.size());
	std::fill(subTypes.begin(), subTypes.end(), t);
}

/**
 * Sets the soil cell index call back function, and defines a rectangular grid for cutting the segments.
 * First cuts all segments at the grid boundaries, then resets and updates the mappers.
 *
 * @param s 		the callback function picks a cell with spatial coordinate [cm] and returns the index of the cell [1]
 * @param min 		minimum of the soil domain [cm]
 * @param max		maximum of the soil domain [cm]
 * @param res	 	resolution, how many cells in each dimension [1]
 */
void MappedSegments::setSoilGrid(const std::function<int(double,double,double)>& s, Vector3d min, Vector3d max, Vector3d res, bool cut) {
	soil_index = s;
	this->setRectangularGrid(min,max,res, cut);
	this->setSoilGrid(s);
}

/**
 * Sets the soil cell index call back function, resets and updates the mappers.
 *
 * @param s 		the callback function picks a cell with spatial coordinate [cm] and returns the index of the cell [1]
 */
void MappedSegments::setSoilGrid(const std::function<int(double,double,double)>& s) {
	soil_index = s;
	seg2cell.clear(); // re-map all segments
	cell2seg.clear();
	mapSegments(segments);
}

/**
 * Sets a rectangular grid, and cuts all segments along the grid cells
 *
 * @param min 		minimum of the soil domain [cm]
 * @param max		maximum of the soil domain [cm]
 * @param res	 	resolution, how many cells in each dimension [1]
 * @param cut 		determines if the segments are cut at the rectangular grid faces
 */
void MappedSegments::setRectangularGrid(Vector3d min, Vector3d max, Vector3d res, bool cut)
{
	minBound = min;
	maxBound = max;
	resolution = res;
	cutAtGrid = cut;
	// std::cout << "setRectangularGrid: cutSegments \n" << std::flush;
	if (cutAtGrid) {
		cutSegments(); // re-add (for cutting)
	}
	// std::cout << "setRectangularGrid: sort \n" << std::flush;
	sort(); // todo should not be necessary, or only in case of cutting?
	seg2cell.clear(); // re-map all segments
	cell2seg.clear();
	// std::cout << "setRectangularGrid: map \n" << std::flush;
	mapSegments(segments);
}


/**
 * Update the mappers root2cell, which maps root segment index to soil cell index, and
 * cell2seg which maps soil cell index to multiple root segments.
 *
 * @param segs      the (new) segments that need to be mapped
 */
void MappedSegments::mapSegments(const std::vector<Vector2i>& segs) {
	for (auto& ns : segs) {
		Vector3d mid = (nodes[ns.x].plus(nodes[ns.y])).times(0.5);
		int cellIdx = soil_index(mid.x,mid.y,mid.z);
		int segIdx = ns.y-1; // this is unique in a tree like structured
		seg2cell[segIdx] = cellIdx;
		if (cell2seg.count(cellIdx)>0) {
			cell2seg[cellIdx].push_back(segIdx);
		} else {
			cell2seg[cellIdx] = std::vector<int>({segIdx});
		}
	}
}

/**
 * Cuts segments @param segs at a rectangular grid (@see MappedSegments::setSoilGrid)
 */
void MappedSegments::cutSegments() {
	assert(segments.size()==radii.size() && "MappedSegments::addSegments: number of segments and radii disagree!");
	assert(segments.size()==subTypes.size() && "MappedSegments::addSegments: number of segments and subTypes disagree!");
	assert(segments.size()==organTypes.size() && "MappedSegments::addSegments: number of segments and organTypes disagree!");
	int n = segments.size(); // segs.size() will change within the loop (recursive implementation)
	for (int i=0; i<n; i++ ) {
		addSegment(segments[i], radii[i], subTypes[i], organTypes[i], i);
	}
}

/**
 * Adds and cuts a single segment at index @param ii. If the segment is cut, appends the remaining segments.
 * Used by cutSegments
 * This approach may run into problems if a segment is located exactly along a face.
 *
 * @param ns 		the segment to add and cut
 * @param r 		segment radius [cm]
 * @param st 		segment sub type
 * @param ot 		segment organ type
 * @param ii		index to insert the segment, -1 to append the segment
 */
void MappedSegments::addSegment(Vector2i ns, double r,  int st, int ot, int ii) {
	Vector3d n1 = nodes[ns.x];
	Vector3d n2 = nodes[ns.y];
	Vector3d mid = (n1.plus(n2)).times(0.5);
	int im = soil_index(mid.x,mid.y,mid.z); // cell indices
	int in1 = soil_index(n1.x,n1.y,n1.z);
	int in2 = soil_index(n2.x,n2.y,n2.z);
	if ((im!=in1) or (im!=in2)) { // cut
		// build SDF
		auto width = maxBound.minus(minBound); // construct sdf
		Vector3d dx(width.x/resolution.x, width.y/resolution.y, width.z/resolution.z);
		auto mid0 = mid.minus(minBound);
		int x = std::floor(mid0.x/dx.x);
		int y = std::floor(mid0.y/dx.y);
		int z = std::floor(mid0.z/dx.z);
		SDF_Cuboid sdf; // create a signed distance function for cutting
		Vector3d minB(x*dx.x, y*dx.y, z*dx.z);
		minB = minB.plus(minBound);
		Vector3d maxB((x+1)*dx.x, (y+1)*dx.y, (z+1)*dx.z);
		maxB = maxB.plus(minBound);
		sdf.min = minB;
		sdf.max = maxB;
		// std::cout << minB.toString() << ", " << maxB.toString() << ", width " << width.toString() << ", " << sdf.min.toString() << ", " << sdf.max.toString() << "\n";

		im = sdf.getDist(mid)>0; // redo indices, since accuracy of pickking may differ
		in1 = sdf.getDist(n1)>0;
		in2 = sdf.getDist(n2)>0;
		if ((im!=in1) or (im!=in2)) {
			Vector3d cPoint;
			if (im==in1) { // is one node at mid (sort accordingly)
				// std::cout << "n1 " << sdf.getDist(n1) << " mid " << sdf.getDist(mid) << " n2 " << sdf.getDist(n2) << ",indices "<< in1 << ", " << im << ", " << in2 << "\n";
				if (sdf.getDist(n2)<0) {
					cPoint = SegmentAnalyser::cut(n2, mid, std::make_shared<SDF_Cuboid>(sdf), eps);
				} else {
					cPoint = SegmentAnalyser::cut(mid, n2, std::make_shared<SDF_Cuboid>(sdf), eps);
				}
				nodeCTs.push_back(nodeCTs[ns.y]); // todo: we might linearly interpolate
			} else if (im==in2) {
				// std::cout << "n2 " << sdf.getDist(n2) << " mid " << sdf.getDist(mid) << " n1 " << sdf.getDist(n1) << ", " << in2 << ", " << im << ", " << in1 << "\n";
				if (sdf.getDist(n1)<0) {
					cPoint = SegmentAnalyser::cut(n1, mid, std::make_shared<SDF_Cuboid>(sdf), eps);
				} else {
					cPoint = SegmentAnalyser::cut(mid, n1, std::make_shared<SDF_Cuboid>(sdf), eps);
				}
				nodeCTs.push_back(nodeCTs[ns.x]); // todo: we might linearly interpolate
			} else { // otherwise split in mid, use cutSegments on those
				cPoint = mid;
				nodeCTs.push_back(0.5*(nodeCTs[ns.x]+nodeCTs[ns.y]));
			}
			// std::cout << "[" << n1.toString() << n2.toString() << "] -> [" << nodes[ns.x].toString() << ", " << nodes.back().toString() << "], ["<< nodes.back().toString() << ", " << n2.toString() << "], " << "\n";
			nodes.push_back(cPoint);
			Vector2i s1(ns.x, nodes.size()-1);
			Vector2i s2(nodes.size()-1, ns.y);
			if ((length(s1)<eps) ||  (length(s2)<eps)) { // if the cut segments are too small, just give up
				add(ns, r, st, ot, ii);
				nodes.pop_back(); // remove cPoint
				nodeCTs.pop_back();
			} else {
				addSegment(s1, r, st , ot, ii); // first segment replaces at index ii
				addSegment(s2, r, st , ot, -1); // append second segment
			}
		} else { // im==in1==in2, dont't cut
			// std::cout << "ok " << ii <<": (" << ns.x <<", " << ns.y << ") [" << n1.toString() <<", "<< n2.toString() <<"\n";
			add(ns, r, st, ot, ii);
		}
	} else { // im==in1==in2, dont't cut
		// std::cout << "ok " << ii <<": (" << ns.x <<", " << ns.y << ") [" << n1.toString() <<", "<< n2.toString() <<"\n";
		add(ns, r, st, ot, ii);
	}
}

/**
 * Adds the segment at index i, or appends it, if i = -1
 *
 * @param s 		the segment to append or insert
 * @param r 		segment radius [cm]
 * @param st 		segment sub type
 * @param ot 		segment organ type
 * @param i			index to insert the segment, -1 to append the segment
 */
void MappedSegments::add(Vector2i s, double r,  int st, int ot,  int i) {
	if (i>=0) {
		segments[i] = s;
		radii[i] = r;
		subTypes[i] = st;
		organTypes[i] = ot;
	} else {
		segments.push_back(s);
		radii.push_back(r);
		subTypes.push_back(st);
		organTypes.push_back(ot);
	}
}

/**
 * Length of the segment @param s
 */
double MappedSegments::length(const Vector2i& s) const {
	return (nodes.at(s.y).minus(nodes.at(s.x))).length();
}

/**
 * Removes segments @param segs from the mappers
 */
void MappedSegments::unmapSegments(const std::vector<Vector2i>& segs) {
	for (auto& ns : segs) {
		int cellIdx = -1;
		int segIdx = ns.y-1;
		if (seg2cell.count(segIdx)>0) { // remove from seg2cell
			cellIdx = seg2cell[segIdx];
			auto it = seg2cell.find(segIdx);
			seg2cell.erase(it);
		} else {
			throw std::invalid_argument("MappedSegments::removeSegments: warning segment index "+ std::to_string(segIdx)+ " was not found in the seg2cell mapper");
		}
		if (cell2seg.count(cellIdx)>0) {
			auto& csegs= cell2seg[cellIdx];
			int c = 0;
			for (int i=0; i<csegs.size(); i++) {
				if (csegs[i] == segIdx) {
					csegs.erase(csegs.begin() + c, csegs.begin() + c);
					break; // inner for
				}
				c++;
			}
		} else {
			throw std::invalid_argument("MappedSegments::removeSegments: warning cell index "+ std::to_string(cellIdx)+ " was not found in the cell2seg mapper");
		}
	}
}

/**
 * Maps a point into a cell and return the cells linear index (for a equidistant rectangular domain)
 */
int MappedSegments::soil_index_(double x, double y, double z) {
	Vector3d p(x,y,z);
	std::array<double,3>  r = { resolution.x, resolution.y, resolution.z};
	auto w = maxBound.minus(minBound);
	auto p0 = p.minus(minBound);
	std::array<double,3> i = { p0.x/w.x*r[0], p0.y/w.y*r[1], p0.z/w.z*r[2] };
	for (int k=0; k<3; k++) {
		if ((i[k] < 0) or (i[k] >= r[k])) {
			return -1; // point is out of domain
		}
	}
	return std::floor(i[2]) * r[0] * r[1] + std::floor(i[1]) * r[0] + std::floor(i[0]); // a linear index not periodic
}

/**
 * Sorts the segments, so that the segment index == second node index -1 (unique mapping in a tree)
 */
void MappedSegments::sort() {
	auto newSegs = segments;
	auto newRadii = radii;
	auto newSubTypes = subTypes;
	auto newTypesorgan = organTypes;
	for (int i=0; i<newSegs.size(); i++) {
		int ind = segments[i].y-1;
		newSegs[ind] = segments[i];
		newRadii[ind] = radii[i];
		newSubTypes[ind] = subTypes[i];
		newTypesorgan[ind] = organTypes[i];
	}
	segments = newSegs;
	radii = newRadii;
	subTypes = newSubTypes;
	organTypes = newTypesorgan;
}

/**
 * Calculates outer segment radii [cm], so that the summed segment volumes per cell equals the cell volume
 * @param type 			prescribed cylinder volume proportional to 0: segment volume, 1: segment surface, 2: segment length
 * @param vols 			(optional) in case of non-equidistant grids, volumes per cell must be defined
 */
std::vector<double> MappedSegments::segOuterRadii(int type, const std::vector<double>& vols) const {
	double cellVolume;
	auto lengths =  this->segLength();
	auto width = maxBound.minus(minBound);
	std::vector<double> outer_radii = std::vector<double>(segments.size());
	std::fill(outer_radii.begin(), outer_radii.end(), 0.);
	for(auto iter = cell2seg.begin(); iter != cell2seg.end(); ++iter) {
		int cellId =  iter->first;
		if (vols.size()==0) {
			cellVolume = width.x*width.y*width.z/resolution.x/resolution.y/resolution.z;
		} else {
			cellVolume = vols.at(cellId);
		}
		auto segs = cell2seg.at(cellId);
		double v = 0.;  // calculate sum of root volumes or surfaces over cell
		for (int i : segs) {
			if (type==0) { // volume
				v += M_PI*(radii[i]*radii[i])*lengths[i];
			} else if (type==1) { // surface
				v += 2*M_PI*radii[i]*lengths[i];
			} else if (type==2) { // length
				v += lengths[i];
			}
		}
		for (int i : segs) { // calculate outer radius
			double l = lengths[i];
			double t =0.; // proportionality factor (must sum up to == 1 over cell)
			if (type==0) { // volume
				t = M_PI*(radii[i]*radii[i])*l/v;
			} else if (type==1) { // surface
				t = 2*M_PI*radii[i]*l/v;
			} else if (type==2) { // length
				t = l/v;
			}
			double targetV = t * cellVolume;  // target volume
			outer_radii[i] = std::sqrt(targetV/(M_PI*l)+radii[i]*radii[i]);
		}
	}
	return outer_radii;
}

/**
 * Calculates segment lengths [cm]
 */
std::vector<double> MappedSegments::segLength() const {
	std::vector<double> lengths = std::vector<double>(segments.size());
	for(int i=0; i<lengths.size(); i++) {
		auto n1 = nodes[segments[i].x];
		auto n2 = nodes[segments[i].y];
		lengths[i] = (n2.minus(n1)).length();
	}
	return lengths;
}


/**
 * Calculates the minimum of node coordinates
 * (e.g. minimum corner of bounding box)
 * value not cached
 */
Vector3d MappedSegments::getMinBounds() {
    Vector3d min_ = Vector3d(nodes[0].x, nodes[0].y, nodes[0].z); 
    for (const auto& n : nodes) {
        if (n.x < min_.x) {
            min_.x = n.x;
        }
        if (n.y < min_.y) {
            min_.y = n.y;
        }
        if (n.z < min_.z) {
            min_.z = n.z;
        }
    }
    return min_;
}



/**
 * Overridden, to map initial shoot segments (@see RootSystem::initialize).
 *
 * Shoot segments have per default radii = 0.1 cm, types = 0, orgtype = 2
 * This can be changed by directly accessing the member variables.
 * @param basaltype			subtype of basal roots 	(default = 4)
 * @param shootbornetype	subtype of shootborn roots (default = 5)
 * @param LB		 		implement length-based waiting time before growth (true) of laterals or delay-based (false)? (default = true)
 */
void MappedRootSystem::initialize_(int basaltype, int shootbornetype, bool verbose, bool LB) {
	std::cout << "MappedRootSystem::initialize \n" << std::flush;
	if(LB){
		RootSystem::initializeLB( basaltype, shootbornetype, verbose);
	}else{RootSystem::initializeDB( basaltype, shootbornetype, verbose);}
	
	segments = this->getShootSegments();
	nodes = this->getNodes();
	nodeCTs = this->getNodeCTs();
	radii.resize(segments.size());
	std::fill(radii.begin(), radii.end(), 0.1);
	subTypes.resize(segments.size());
	std::fill(subTypes.begin(), subTypes.end(), 0);
	organTypes.resize(segments.size());
	std::fill(organTypes.begin(), organTypes.end(), Organism::ot_root); //root organ type = 2
	mapSegments(segments);
}



/**
 * Simulates the development of the organism in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void MappedRootSystem::simulate(double dt, bool verbose)
{
	if (soil_index==nullptr) {
		throw std::invalid_argument("MappedRootSystem::simulate():soil was not set, use MappedRootSystem::simulate::setSoilGrid" );
	}

	RootSystem::simulate(dt,verbose);

	auto uni = this->getUpdatedNodeIndices(); // move nodes
	auto unodes = this->getUpdatedNodes();
	assert(uni.size()==unodes.size() && "updated node indices and number of nodes must be equal");
	int c = 0;
	for (int i : uni) {
		nodes.at(i) = unodes[c];
		c++;
	}
	if (verbose) {
		std::cout << "nodes moved "<< uni.size() << "\n" << std::flush;
	}
	auto newnodes = this->getNewNodes(); // add nodes
	nodes.reserve(nodes.size()+newnodes.size());
	for (auto& nn : newnodes) {
		nodes.push_back(nn);
	}
	auto newnode_cts = this->getNewNodeCTs(); // add node cts
	nodeCTs.reserve(nodeCTs.size()+newnode_cts.size());
	for (auto& nct : newnode_cts) {
		nodeCTs.push_back(nct);
	}
	if (verbose) {
		std::cout << "new nodes added " << newnodes.size() << "\n" << std::flush;
	}
	auto newsegs = this->getNewSegments(); // add segments (TODO cutting)
	segments.resize(segments.size()+newsegs.size());
	for (auto& ns : newsegs) {
		segments[ns.y-1] = ns;
	}
	if (verbose) {
		std::cout << "segments added "<< newsegs.size() << "\n" << std::flush;
	}
	auto newsegO = this->getNewSegmentOrigins(); // to add radius and type (TODO cutting)
	radii.resize(radii.size()+newsegO.size());
	subTypes.resize(subTypes.size()+newsegO.size());
	organTypes.resize(organTypes.size()+newsegO.size());
	c = 0;
	if (verbose) {
		std::cout << "Number of segments " << radii.size() << ", including " << newsegO.size() << " new \n"<< std::flush;
	}
	for (auto& so : newsegO) {
		int segIdx = newsegs[c].y-1;
		c++;
		radii[segIdx] = so->getParam()->a;
		subTypes[segIdx] = so->getParam()->subType;
		organTypes[segIdx] = so->organType();
	}
	// map new segments
	this->mapSegments(newsegs);

	// update segments of moved nodes
	std::vector<Vector2i> rSegs;
	for (int i : uni) {
		int segIdx = i -1;
		int cellIdx = seg2cell[segIdx];
		auto s = segments[segIdx];
		Vector3d mid = (nodes[s.x].plus(nodes[s.y])).times(0.5);
		int newCellIdx = soil_index(mid.x,mid.y,mid.z);
		// 1. check if mid is still in same cell (otherwise, remove, and add again)
		// 2. if cut is on, check if end point is in same cell than mid point (otherwise remove and add again)
		bool remove = false;
		if (cellIdx==newCellIdx) {
			if (cutAtGrid) {
				auto endPoint = nodes[s.y];
				newCellIdx = soil_index(endPoint.x,endPoint.y,endPoint.z);
				remove = (newCellIdx!=cellIdx);
			}
		} else {
			remove = true;
		}
		if (remove) {
			rSegs.push_back(s);
		}
	}
	MappedSegments::unmapSegments(rSegs);
	MappedSegments::mapSegments(rSegs);
}


/**
 * initialization of mappedplant
 * @param verbose 		indicates if status is written to the console (cout) (default = false)
 * @param stochastic 	keep stochasticity in simulation? (default = true)
 * @param LB		 	implement length-based waiting time before growth (true) of laterals or delay-based (false)? (default = true)
 */
void MappedPlant::initialize_(bool verbose, bool stochastic, bool LB)
{
  reset(); // just in case
	std::cout << "MappedPlant::initialize \n" << std::flush;
	this->stochastic = stochastic;
	if(LB){	Plant::initializeLB(verbose);
	}else{Plant::initializeDB(verbose);}
	Plant::setStochastic(stochastic);
	nodes = this->getNodes();
	nodeCTs = this->getNodeCTs();
	mapSegments(segments);
	mapSubTypes();
	plantParam = this->organParam;
}

/**
 * creates a map to reorder sub types, so that
 * the N subtypes of one organ type go from 0 to N-1
 */
void MappedPlant::mapSubTypes(){
	for(int ot = 0; ot < organParam.size();ot++ )
	{
		//std::cout<<"MappedPlant::mapSubTypes for organtype "<<ot<<" with "<<organParam[ot].size()<<" subtypes "<<std::endl;
		int stNew = 0;
		for(int stOld_ = 1; stOld_ < organParam[ot].size();stOld_++)//skipe stOld ==0, not realted to any organ st
		{
			if(organParam[ot][stOld_] != NULL) {
				int stOld = organParam[ot][stOld_]->subType;
				st2newst[std::make_tuple(ot, stOld)] = stNew;
				//std::cout<<"old st: "<<stOld<<", new st: "<< stNew <<std::endl;
				stNew ++;
			}// else {std::cout<<"subType n#"<<stOld_<<" does not exist, skip "<<std::endl;}
		}
	}
}

/**
 * Simulates the development of the organism in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void MappedPlant::simulate(double dt, bool verbose)
{
	if (soil_index==nullptr) {
		throw std::invalid_argument("MappedPlant::simulate():soil was not set, use MappedPlant::simulate::setSoilGrid" );
	}
	Plant::simulate( dt,  verbose);
	auto uni = this->getUpdatedNodeIndices(); // move nodes
	auto unodes = this->getUpdatedNodes();
	auto uncts = this->getUpdatedNodeCTs();
	assert(uni.size()==unodes.size() && "updated node indices and number of nodes must be equal");
	int c = 0;
	for (int i : uni) {
		nodes.at(i) = unodes[c];
		nodeCTs.at(i) = uncts[c];
		c++;
	}

	if (verbose) {
		std::cout << "nodes moved "<< uni.size() << "\n" << std::flush;
	}
	auto newnodes = this->getNewNodes(); // add nodes
	nodes.reserve(nodes.size()+newnodes.size());
	for (auto& nn : newnodes) {
		nodes.push_back(nn);
	}
	auto newnode_cts = this->getNewNodeCTs(); // add node cts
	nodeCTs.reserve(nodeCTs.size()+newnode_cts.size());
	for (auto& nct : newnode_cts) {
		nodeCTs.push_back(nct);
	}
	if (verbose) {
		std::cout << "new nodes added " << newnodes.size() << "\n" << std::flush;
	}
	auto newsegs = this->getNewSegments(); // add segments (TODO cutting)
	segments.resize(segments.size()+newsegs.size());
	for (auto& ns : newsegs) {
		segments[ns.y-1] = ns;
	}
	if (verbose) {
		std::cout << "segments added "<< newsegs.size() << "\n" << std::flush;
	}
	auto newsegO = this->getNewSegmentOrigins(); // to add radius and type (TODO cutting)
	radii.resize(radii.size()+newsegO.size());
	subTypes.resize(subTypes.size()+newsegO.size());
	organTypes.resize(organTypes.size()+newsegO.size());
	segVol.resize(segVol.size()+newsegO.size());
	bladeLength.resize(bladeLength.size()+newsegO.size());
	leafBladeSurface.resize(leafBladeSurface.size()+newsegO.size()); 
	c = 0;
	if (verbose) {
		std::cout << "Number of segments " << radii.size() << ", including " << newsegO.size() << " new \n"<< std::flush;
	}
	std::vector<int> vsegIdx;
	for (auto& so : newsegO) {
		int segIdx = newsegs[c].y-1;
		vsegIdx.push_back(segIdx);
		radii[segIdx] = so->getParam()->a;
		organTypes.at(segIdx) = so->organType();
		subTypes.at(segIdx) = st2newst[std::make_tuple(organTypes[segIdx],so->getParam()->subType)];//new st 
		
		if(organTypes[segIdx] == Organism::ot_leaf) //leaves can be cylinder, cuboid or characterized by user-defined 2D shape
		{
			int index;
			auto nodeIds = so->getNodeIds();
			auto it = find(nodeIds.begin(), nodeIds.end(), newsegs[c].y);
			if (it != nodeIds.end()){ index = it - nodeIds.begin() -1;
			}else { 
				throw std::runtime_error("MappedPlant::simulate: global segment index not found in organ");
			}
			int localSegId = index;
			bool realized = true; bool withPetiole = false;
			segVol.at(segIdx) = -1; 
			bladeLength.at(segIdx) = std::static_pointer_cast<Leaf>(so)->leafLengthAtSeg(localSegId, withPetiole);
			leafBladeSurface.at(segIdx) =  std::static_pointer_cast<Leaf>(so)->leafAreaAtSeg(localSegId,realized, withPetiole);
			withPetiole = true;
			segVol.at(segIdx) = std::static_pointer_cast<Leaf>(so)->leafVolAtSeg(localSegId, realized, withPetiole);//* thickness;
			assert((segVol.at(segIdx) >= 0)&&"MappedPlant::simulate: computation of leaf volume failed");
			
		}else{ //stems and roots are cylinder
			auto s = segments.at(segIdx);
			double length_seg = (nodes.at(s.x).minus(nodes.at(s.y))).length();
			segVol.at(segIdx) = radii.at(segIdx) * radii.at(segIdx) * M_PI * length_seg;
			bladeLength.at(segIdx) = 0;
			leafBladeSurface.at(segIdx) = 0;
		}
		c++;
	}

	// map new segments
	this->mapSegments(newsegs);

	// update segments of moved nodes
	std::vector<Vector2i> rSegs;
	for (int i : uni) {
		int segIdx = i -1;
		int cellIdx = seg2cell[segIdx];
		auto s = segments[segIdx];
		Vector3d mid = (nodes[s.x].plus(nodes[s.y])).times(0.5);
		int newCellIdx = soil_index(mid.x,mid.y,mid.z);
		// 1. check if mid is still in same cell (otherwise, remove, and add again)
		// 2. if cut is on, check if end point is in same cell than mid point (otherwise remove and add again)
		bool remove = false;
		if (cellIdx==newCellIdx) {
			if (cutAtGrid) {
				auto endPoint = nodes[s.y];
				newCellIdx = soil_index(endPoint.x,endPoint.y,endPoint.z);
				remove = (newCellIdx!=cellIdx);
			}
		} else {
			remove = true;
		}
		if (remove) {
			rSegs.push_back(s);
		}
	}
	MappedSegments::unmapSegments(rSegs);
	MappedSegments::mapSegments(rSegs);
	if(kr_length > 0.){calcExchangeZoneCoefs();}

}



/**
 * computes coeficients for kr
 * when root kr > 0 up to kr_length cm from the root tip
 * see @XylemFlux::kr_RootExchangeZonePerType()
 **/
void MappedPlant::calcExchangeZoneCoefs() { //
	exchangeZoneCoefs.resize(segments.size(), -1.0);
	auto orgs = getOrgans(-1);
	for(auto org: orgs)
	{
		for(int localIdx = 1; localIdx < org->getNumberOfNodes();localIdx++)
		{
			int globalIdx_x = org->getNodeId(localIdx -1 );
			int globalIdx_y = org->getNodeId(localIdx);
			if(org->organType() != Organism::ot_root){exchangeZoneCoefs.at(globalIdx_y-1) = 1;
			}else{
				auto n1 = nodes.at(globalIdx_x);
				auto n2 = nodes.at(globalIdx_y);
				auto v = n2.minus(n1);
				double l = v.length();
				double distance2RootTip_y = org->getLength(true) - org->getLength(localIdx);
				double length_in_exchangeZone = std::min(l,std::max(kr_length - std::max(distance2RootTip_y,0.),0.));
				exchangeZoneCoefs.at(globalIdx_y-1) = length_in_exchangeZone/l;
			}
		}
	}
	const int notFound = std::count(exchangeZoneCoefs.cbegin(), exchangeZoneCoefs.cend(), -1.0);
	if(notFound != 0)
	{
		std::stringstream errMsg;
		errMsg <<"MappedPlant::calcExchangeZoneCoefs(): "<<notFound<<" elements not initalized";
		throw std::runtime_error(errMsg.str().c_str());
		std::cout<<"notFound "<<notFound<<std::endl;
	}
}


/**
 *Gives an overview of the mappedplant object (for debugging)
 *
 **/
void MappedPlant::printNodes() {

	std::cout << "\n MappedPlant::printnodes \n"<< std::flush;
	std::cout << "\n nodes \n"<< std::flush;
	nodes = this->getNodes();
	for (auto nd : nodes) {
		std::cout <<nd.toString()<< std::flush;
	}
	std::cout << "\n nodes size \n" <<  nodes.size() << std::flush;
	std::cout << "\n organ types \n" << std::flush;
	for (auto ot : organTypes) {
		std::cout << ot << std::flush;
	}
	std::cout << "\n organTypes size\n"<<  organTypes.size() << std::flush;
	std::cout << "\n subtypes \n"<< std::flush;
	for (auto st : subTypes) {
		std::cout << st << std::flush;
	}
	std::cout << "\n subtypes size \n"<< subTypes.size() << std::flush;
	std::cout << "\n cts \n"<< std::flush;
	for (auto to : nodeCTs) {
		std::cout << to << std::flush;
	}
	std::cout << "\n cts size \n"<< nodeCTs.size() << std::flush;
	std::cout << "\n segments \n"<< std::flush;
	for (auto to : segments) {
		std::cout << to.toString() << std::flush;
	}
	std::cout << "\n segments size \n"<< segments.size() << std::flush;
}

/**
 *	index of node of organtype ot
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          Id of each segment
 */
std::vector<int> MappedPlant::getSegmentIds(int ot) const
{
	std::vector<int> segId;// = std::vector<int>(segments.size());
    for (int i=0; i<segments.size(); i++) {
		if((ot == -1)||(organTypes[segments[i].y-1]== ot)){
			segId.push_back(segments[i].y-1);
		}
    }
	return segId;
}

/**
 * index of segment of organ type ot
 * @param ot        the expected organ type, where -1 denotes all organ types (default)
 * @return          Id of each segment
 */
std::vector<int> MappedPlant::getNodeIds(int ot) const
{
	std::vector<int> nodeId;// = std::vector<int>(segments.size());
    for (int i=0; i<segments.size(); i++) {
		if((ot == -1)||(organTypes[segments[i].y-1]== ot)){
			nodeId.push_back(segments[i].y);
		}
    }
	return nodeId;
}

void MappedPlant::ComputeGeometryForOrgan(int organId)
{
  auto result = this->getOrgans(-1, false);
  auto organ_it = std::find_if(result.begin(), result.end(), [organId](const auto& o) {
    return o->getId() == organId;
  });
  if(organ_it == result.end())
  {
    std::stringstream errMsg;
    errMsg << "MappedPlant::ComputeGeometryForOrgan: organ not found: " << organId;
    throw std::runtime_error(errMsg.str().c_str());
  }
  auto organ = *organ_it;
  
  if(organ->organType() == 4)
  {
    // 4 SHOULD mean leaf, so we do not check for successful cast
    GenerateLeafGeometry(std::dynamic_pointer_cast<Leaf>(organ));
  }
  else
  {
    GenerateStemGeometry(organ);
  }
}

void MappedPlant::ComputeGeometryForOrganType(int organType)
{
  auto organ_list = this->getOrgans(organType, false);
  // First we check if we have enough memory to support the geometry
  unsigned int point_space = 0;
  unsigned int cell_space = 0;
  for(auto organ : organ_list)
  {
    // Check against object, because organType can be -1
    if(organ->organType() == 4)
    {
      // 4 SHOULD mean leaf, so we do not check for successful cast
      point_space += organ->getNumberOfNodes() * 6 * 3;
      cell_space += organ->getNumberOfNodes() * 6;
    }
    else
    {
      point_space += organ->getNumberOfNodes() * 3 * geometry_resolution;
      cell_space += (organ->getNumberOfNodes() - 1) * 2 * geometry_resolution;
    }
  }
  geometry.reserve(point_space);
  geometryNormals.reserve(point_space);
  geometryColors.reserve(point_space);
  geometryIndices.reserve(cell_space);
  geometryTextureCoordinates.reserve(point_space);

  point_space = 0;
  cell_space = 0;
  for(auto organ : organ_list)
  {
    if(organType >= 0 && organ->organType() != organType)
    {
      continue;
    }
    if(organ->organType() == 4)
    {
      // 4 SHOULD mean leaf, so we do not check for successful cast
      GenerateLeafGeometry(std::dynamic_pointer_cast<Leaf>(organ), point_space, cell_space);
      point_space += organ->getNumberOfNodes() * 6 * 3;
      cell_space += organ->getNumberOfNodes() * 6;
    }
    else
    {
      GenerateStemGeometry(organ, point_space, cell_space);
      point_space += organ->getNumberOfNodes() * 3 * geometry_resolution;
      cell_space += (organ->getNumberOfNodes() - 1) * 2 * geometry_resolution;
    }
  }
}

void MappedPlant::ComputeGeometry()
{
  this->ComputeGeometryForOrganType(-1);
}

void MappedPlant::MapPropertyToColors(std::string property, std::vector<double> minMax)
{
  
}

void MappedPlant::GenerateLeafGeometry(std::shared_ptr<Leaf> leaf, unsigned int p_o, unsigned int c_o)
{
  // std::vector::reserve should be idempotent.
  geometry.reserve(p_o + leaf->getNumberOfNodes() * 2);
  geometryNormals.reserve(p_o + leaf->getNumberOfNodes() * 2);
  geometryIndices.reserve(c_o + (leaf->getNumberOfNodes() - 1) * 12);
  geometryColors.reserve(p_o + leaf->getNumberOfNodes() * 2 );
  geometryTextureCoordinates.reserve(p_o + leaf->getNumberOfNodes() * 2);
  geometryNodeIds.reserve(p_o + leaf->getNumberOfNodes() * 2);
  unsigned int points_index = p_o;
  unsigned int cell_index = c_o;
  Quaternion heading = Quaternion::FromMatrix3d(leaf->iHeading);
  Quaternion rot = Quaternion::geodesicRotation({1,0,0},{0,0,1});
  Vector3d lastPosition = leaf->getNode(0);
  // TODO: check if getLength(true) is the correct length
  double totalLenght = leaf->getLength(true);
  double currentLength = 0;
  for(int i = 0; i < leaf->getNumberOfNodes(); ++i)
  {
    // we presume that the nodes follow the center line of the plant
    // and that the leaf is oriented along the x axis
    auto position = leaf->getNode(i);
    auto id = leaf->getNodeId(i);
    Vector3d dist;
    // This part is for the rotation of the leaf segment
    if(i + 1 < leaf->getNumberOfNodes())
    {
      dist = leaf->getNode(i + 1) - position;
    }
    else
    {
      dist = leaf->getNode(i - 1) - position;
    }
    // This, in contrast, is for the texture mapping
    currentLength += (position - lastPosition).length() / totalLenght;
    rot *= Quaternion::geodesicRotation(rot.Forward(), dist);
    // TODO: Check with mona on what the Vector3d coordinates of
    // this function are, and if we need to change them
    auto vis = leaf->getLeafVis(id);
    // We don't normally split normals, but in this case we have a flat surface
    // and we want to have a smooth shading
    geometryTextureCoordinates.insert(geometryTextureCoordinates.begin() + points_index,
                                            {currentLength, 0.0, currentLength, 1.0});
    geometryNodeIds.insert(geometryNodeIds.begin() + points_index, {id, id});
    // TODO: it not obvious that points_index can be changed by the insert here
    vec2Buf(geometry, points_index, vis[0], vis[1]);
    points_index = vec2Buf(geometryNormals, points_index, rot.Up(), rot.Up());
    vec2Buf(geometry, points_index, vis[0], vis[1]);
    points_index = vec2Buf(geometryNormals, points_index ,-rot.Up(), -rot.Up());
    // The triangles are defined clockwise for the front face and counter clockwise for the back face
    if(i > 0)
    {
      geometryIndices.insert(geometryIndices.begin() + cell_index, 
      {points_index-6, points_index-7, points_index-3, points_index-6, points_index-3, points_index-2,
       points_index-4, points_index-5, points_index-1, points_index-4, points_index-1, points_index-0});
      cell_index += 12;
    }
  }
}

void MappedPlant::GenerateStemGeometry(std::shared_ptr<Organ> stem, unsigned int p_o, unsigned int c_o)
{
  const unsigned int geometry_resolution = 16;
  
  geometry.reserve(p_o + (stem->getNumberOfNodes() * geometry_resolution * 3));
  geometryNormals.reserve(p_o + stem->getNumberOfNodes() * geometry_resolution * 3);
  geometryIndices.reserve(c_o + (stem->getNumberOfNodes() -1) * geometry_resolution * 3);
  geometryColors.reserve(p_o + stem->getNumberOfNodes() * geometry_resolution);
  geometryTextureCoordinates.reserve(p_o + stem->getNumberOfNodes() * geometry_resolution * 2);
  geometryNodeIds.reserve(p_o + stem->getNumberOfNodes() * geometry_resolution);
  Quaternion heading = Quaternion::FromMatrix3d(stem->iHeading);
  Quaternion lastRotation = Quaternion::geodesicRotation({1,0,0},{0,0,1});
  for(auto i = 0; i < stem->getNumberOfNodes(); ++i)
  {
    double diameter = stem->getParameter("radius");
    const auto& node = stem->getNode(i);
    Vector3d dist;
    if(i + 1 < stem->getNumberOfNodes())
    {
      dist = stem->getNode(i + 1) - node;
    }
    else
    {
      dist = stem->getNode(i - 1) - node;
    }
    lastRotation *= Quaternion::geodesicRotation(lastRotation.Forward(), dist);
    auto deltaphi = 2 * M_PI / geometry_resolution;
    for (auto j = 0; j < geometry_resolution; ++j)
    {
      auto phi = j * deltaphi;
      Vector3d outer = {0.0, std::cos(phi) * diameter, std::sin(phi) * diameter};
      outer = lastRotation.Rotate(outer);
      // consequtive points are stored in the buffer
      // the plus 0 is not necessary, but it makes it easier to read
      geometry[3 * (i * geometry_resolution + j) + 0 + p_o] = node.x + outer.x;
      geometry[3 * (i * geometry_resolution + j) + 1 + p_o] = node.y + outer.y;
      geometry[3 * (i * geometry_resolution + j) + 2 + p_o] = node.z + outer.z;
      // the indices are stored in the buffer, and are all front facing
      geometryIndices[6 * (i * geometry_resolution + j) + 0 + c_o] = i * geometry_resolution + j;
      geometryIndices[6 * (i * geometry_resolution + j) + 1 + c_o] = i * geometry_resolution + (j + 1) % geometry_resolution;
      geometryIndices[6 * (i * geometry_resolution + j) + 2 + c_o] = (i + 1) * geometry_resolution + (j + 1) % geometry_resolution;
      geometryIndices[6 * (i * geometry_resolution + j) + 3 + c_o] = i * geometry_resolution + j;
      geometryIndices[6 * (i * geometry_resolution + j) + 4 + c_o] = (i + 1) * geometry_resolution + (j + 1) % geometry_resolution;
      geometryIndices[6 * (i * geometry_resolution + j) + 5 + c_o] = (i + 1) * geometry_resolution + j;
      // the normals are stored in the buffer
      geometryNormals[3 * (i * geometry_resolution + j) + 0 + p_o] = outer.x;
      geometryNormals[3 * (i * geometry_resolution + j) + 1 + p_o] = outer.y;
      geometryNormals[3 * (i * geometry_resolution + j) + 2 + p_o] = outer.z;

      geometryNodeIds[i * geometry_resolution + j + p_o] = stem->getNodeId(i);

      geometryTextureCoordinates[2 * (i * geometry_resolution + j) + 0 + p_o] = i / (double)stem->getNumberOfNodes();
      geometryTextureCoordinates[2 * (i * geometry_resolution + j) + 1 + p_o] = phi / (2 * M_PI);
    }
  }
}

} // namespace
