// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "MappedOrganism.h"

#include "SegmentAnalyser.h"
#include "growth.h"
#include <algorithm>
#include <functional>
#include <iterator>
#include <cmath>
#include "mymath.h"

namespace CPlantBox {

	// an iterator that can give me the vector mirrored around 0
	class MirrorIterator : public std::iterator<std::input_iterator_tag, std::pair<int,double> > {
	public:
		MirrorIterator(const std::vector<double>* v) : v(v) {
			//std::cout << "MirrorIterator was created " << v->size() << std::endl;
			// output all vector elements
      std::copy(v->begin(), v->end(), std::ostream_iterator<double>(std::cout, " ")); std::cout << std::endl;
		}
		std::pair<int,double> operator*() { return std::make_pair(idx(), v->at(i)); }
		MirrorIterator& operator++() { inc(); return *this; }
		bool operator!=(const MirrorIterator& other) { return i != other.i && r != other.r; }
    bool operator==(const MirrorIterator& other) { return i == other.i && r == other.r; }
		std::size_t size() { return (v->size() > 0) ? (v->size() * 2 - (v->back() < std::numeric_limits<float>::epsilon() ? 1 : 0)) : 0; }
		double operator[](int i) {
			if(i < v->size())
				return v->at(i);
			else
			{
				return -v->at(v->size() - (i - v->size() + 1) - ((v->back() < std::numeric_limits<float>::epsilon()) ? 1 : 0));
			}
		 }
		// begin
		MirrorIterator begin() { return MirrorIterator(v); }
		// end
		MirrorIterator end() { return MirrorIterator(v, true, 0); }

    // computes the texture coordinate within the unit interval based on the index
    double texcoord(int i)
		{
		  double d = static_cast<double>(i) / static_cast<double>(size() - 1);
			return d;
		}
		void inc()
		{
			if(r)
			{
				--i;
			}
			else
			{
				++i;
				if(i == v->size())
				{
					r = true;
					i = v->size() - (v->back() < std::numeric_limits<float>::epsilon() ? 2 : 1);
				}
			}
		}
		// checks whether we are in the mirrored part for a given index
		bool isMirrored(int i)
		{
		  return i >= v->size();
    }
		int idx()
		{
			if(r)
				return v->size() + i - (v->back() < std::numeric_limits<float>::epsilon() ? 1 : 0);
			else
				return i;
		}
	private:
		MirrorIterator(const std::vector<double>* v, bool end, int i) : v(v), r(end), i(i) { }
		bool r = false;
		const std::vector<double>* v;
		int i = 0;
	};

  /** 
   * Hidden Template method to write a set of vectors into a buffer
   * @param buffer    the buffer
   * @param offset    the offset in the buffer
   * @param v...     the vectors
   */
  inline unsigned int vec2Buf(std::vector<double>& buffer, unsigned int offset, const Vector3d& v) {
    buffer[offset + 0] = v.x;
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
	// //std::cout << "setRectangularGrid: cutSegments \n" << std::flush;
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
	if ((im!=in1) || (im!=in2)) { // cut
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
		if ((im!=in1) || (im!=in2)) {
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
		if ((i[k] < 0) || (i[k] >= r[k])) {
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
	//std::cout << "MappedRootSystem::initialize \n" << std::flush;
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
		//std::cout << "nodes moved "<< uni.size() << "\n" << std::flush;
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
		//std::cout << "new nodes added " << newnodes.size() << "\n" << std::flush;
	}
	auto newsegs = this->getNewSegments(); // add segments (TODO cutting)
	segments.resize(segments.size()+newsegs.size());
	for (auto& ns : newsegs) {
		segments[ns.y-1] = ns;
	}
	if (verbose) {
		//std::cout << "segments added "<< newsegs.size() << "\n" << std::flush;
	}
	auto newsegO = this->getNewSegmentOrigins(); // to add radius and type (TODO cutting)
	radii.resize(radii.size()+newsegO.size());
	subTypes.resize(subTypes.size()+newsegO.size());
	organTypes.resize(organTypes.size()+newsegO.size());
	c = 0;
	if (verbose) {
		//std::cout << "Number of segments " << radii.size() << ", including " << newsegO.size() << " new \n"<< std::flush;
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
	//std::cout << "MappedPlant::initialize \n" << std::flush;
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
		unsigned int point_space = 0, cell_space = 0;
		
		GenerateStemGeometry(organ, point_space, cell_space);
		point_space += organ->getNumberOfNodes() * 3 * geometry_resolution;
		cell_space += (organ->getNumberOfNodes() - 1) * 2 * geometry_resolution;
		auto leaf = std::dynamic_pointer_cast<Leaf>(organ);
		if(leaf->getLeafRandomParameter()->parametrisationType == 1)
		{
			GenerateRadialLeafGeometry(leaf, point_space, cell_space);
			point_space += leaf->getNumberOfNodes() * 6 * 3;
			cell_space += leaf->getNumberOfNodes() * 6;
		}
		else
		{
			int petiole_zone = 0;
			for(int i = 0; i < leaf->getNumberOfNodes(); i++)
			{
				if(leaf->nodeLeafVis(leaf->getLength(i)))
				{
					petiole_zone = i;
					break;
				}
			}
			if(petiole_zone + 1 < leaf->getNumberOfNodes())
			{
				GenerateLeafGeometry(leaf, petiole_zone, point_space, cell_space);
				point_space += (organ->getNumberOfNodes() - petiole_zone) * 6 * 3;
				cell_space += (organ->getNumberOfNodes() - petiole_zone) * 6;
			}
		}
  }
  else
  {
    GenerateStemGeometry(organ);
  }
}

void MappedPlant::ComputeGeometryForOrganType(int organType)
{
  auto organ_list = this->getOrgans(-1, false);
		
  // First we check if we have enough memory to support the geometry
  unsigned int point_space = 0;
  unsigned int cell_space = 0;
	unsigned int num_organs = 0;

  for(auto organ : organ_list)
  {
    // Check against object, because organType can be -1
		if(organ->organType() == organType || organType < 0)
		{
			if(organ->organType() == 4)
			{
				// 4 SHOULD mean leaf, so we do not check for successful cast
				if(bIncludeMidlineInLeaf)
				{
				  point_space += organ->getNumberOfNodes() * 3 * geometry_resolution;
				  cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution;
				}
			  point_space += (organ->getNumberOfNodes()) * 4 * 3;
				cell_space += ((organ->getNumberOfNodes()) - 1) * 4 * 3 + 40;
			}
			else
			{
				point_space += organ->getNumberOfNodes() * 3 * geometry_resolution;
				cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution;
			}
		}
  }
  //std::cout << "Going to allocate " << point_space << " points and " << cell_space << " cells" << std::endl;
  geometry.clear();
  geometryNormals.clear();
  geometryColors.clear();
  geometryIndices.clear();
  geometryTextureCoordinates.clear();
  geometryNodeIds.clear();
  geometry.reserve(point_space);
  geometryNormals.reserve(point_space);
  geometryColors.reserve(point_space / 3 * 2);
  geometryIndices.reserve(cell_space);
  geometryTextureCoordinates.reserve(point_space);
  geometryNodeIds.reserve(point_space / 3 + 1);


  point_space = 0;
  cell_space = 0;
  unsigned int checked_organs = 0;
  for(auto organ : organ_list)
  {
    checked_organs++;
		//std::cout << "Going through organ " << organ->getId() << std::endl;

    if((organType >= 1 && organ->organType() != organType) || organ->getNumberOfNodes() <= 1)
    {
      continue;
    }
    if(organ->organType() == 4)
    {
			auto leaf = std::dynamic_pointer_cast<Leaf>(organ);
			// render petiole
      //std::cout << "Generating geometry for leaf " << organ->getId() << " with " << organ->getNumberOfNodes() << " nodes." << std::endl;
      //std::cout << "Stem part for petiole and rest" << std::endl;
			if(bIncludeMidlineInLeaf)
			{
        GenerateStemGeometry(organ, point_space, cell_space);
			}
      //std::cout << "Updating buffer positions because the leaf is a two-parter" << std::endl;
      //point_space += organ->getNumberOfNodes() * 3 * geometry_resolution;
      //cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution;
      point_space = geometry.size();
      cell_space = geometryIndices.size();
			
			if(leaf->getLeafRandomParameter()->parametrisationType == 1)
			{
				std::cout << "Generating radial leaf geometry" << std::endl;
				GenerateRadialLeafGeometry(leaf, point_space, cell_space);
				point_space = geometry.size();
				cell_space = geometryIndices.size();
			}
			else
			{
				int petiole_zone = 0;
				for(int i = 0; i < leaf->getNumberOfNodes(); i++)
				{
					if(!leaf->nodeLeafVis(leaf->getLength(i)))
					{
						petiole_zone = i;
					}
					else break;
				}
				if(petiole_zone + 1 < leaf->getNumberOfNodes())
				{
					GenerateLeafGeometry(leaf, petiole_zone, point_space, cell_space);
					point_space = geometry.size();
					cell_space = geometryIndices.size();
				}
			}
    }
    else
    {
			//std::cout << "Organ is a stem" << std::endl;
      //std::cout << "Generating geometry for stem " << organ->getId() << " with " << organ->getNumberOfNodes() << " nodes." << std::endl;
      auto prev_size = geometryIndices.size();
      GenerateStemGeometry(organ, point_space, cell_space);
      //std::cout << "Organ " << organ->getId() << " pushed the size from " << prev_size << " to " << geometryIndices.size() << std::endl;
      //point_space += organ->getNumberOfNodes() * 3 * geometry_resolution;
      //cell_space += (organ->getNumberOfNodes() - 1) * 6 * geometry_resolution;
      point_space = geometry.size();
      cell_space = geometryIndices.size();
    }
    //std::cout << "Done generating geometry for organ " << organ->getId() << std::endl;

  }
  assert(geometry.size() % 3 == 0);
  assert(geometryIndices.size() % 3 == 0);
}

void MappedPlant::ComputeGeometry()
{
  this->ComputeGeometryForOrganType(-1);
	std::cout << "Sanity Check for C++ " << std::endl;
	std::cout << "Geometry size: " << geometry.size() << std::endl;
	std::cout << "Geometry Indices size: " << geometryIndices.size() << std::endl;
	std::cout << "Geometry Normals size: " << geometryNormals.size() << std::endl;
	std::cout << "Geometry Colors size: " << geometryColors.size() << std::endl;
	std::cout << "Geometry Texture Coordinates size: " << geometryTextureCoordinates.size() << std::endl;
	std::cout << "Geometry Node Ids size: " << geometryNodeIds.size() << std::endl;
	std::cout << "This would mean we have " << geometry.size() / 3 << " points and " << geometryIndices.size() / 3 << " cells." << std::endl;
}

void MappedPlant::MapPropertyToColors(std::string property, std::vector<double> minMax)
{
  
}

void MappedPlant::GenerateLeafGeometry(std::shared_ptr<Leaf> leaf, unsigned int petiole_zone, unsigned int p_o, unsigned int c_o)
{
  // std::vector::reserve should be idempotent.
  //std::cout << "Resizing geometry buffers for a leaf with n=" << leaf->getNumberOfNodes() << ", pet=" << petiole_zone << std::endl;
	int first_surface_id = petiole_zone + 1;
	int total_points = leaf->getNumberOfNodes() - first_surface_id;
	geometry.resize(std::max(static_cast<std::size_t>(p_o + total_points * 4 * 3), geometry.size()), -1.0);
	geometryNormals.resize(std::max(static_cast<std::size_t>(p_o + total_points * 4 * 3), geometryNormals.size()), -1.0);
	geometryIndices.resize(std::max(static_cast<std::size_t>(c_o + (total_points - 1) * 12), geometryIndices.size()), static_cast<unsigned int>(-1));
	geometryColors.resize(std::max(static_cast<std::size_t>((p_o/3) + total_points * 4 * 3), geometryColors.size()), static_cast<unsigned short>(-1));
	geometryTextureCoordinates.resize(std::max(static_cast<std::size_t>((p_o/3*2) + total_points * 4 * 2), geometryTextureCoordinates.size()), -1.0);
	geometryNodeIds.resize(std::max(static_cast<std::size_t>(p_o/3 + total_points * 4), geometryNodeIds.size()), -1);

  //std::cout << "Orientation generation" << std::endl;
  unsigned int points_index = p_o;
  unsigned int cell_index = c_o;
  Quaternion heading = Quaternion::FromMatrix3d(leaf->iHeading);
  Quaternion rot = Quaternion::geodesicRotation({1,0,0},{0,0,1}).normalized();
  Vector3d lastPosition = leaf->getNode(first_surface_id);
  // TODO: check if getLength(true) is the correct length
  double totalLenght = leaf->getLength(true);
  double currentLength = leaf->getLength(first_surface_id);
  for(int i = first_surface_id; i < leaf->getNumberOfNodes(); ++i)
  {
    // we presume that the nodes follow the center line of the plant
    //std::cout << " going through node id " << i << " on this leaf." << std::endl;
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
		//std::cout << "[Leaf] Rotating " << rot.toString() << " to get " << dist.toString() << std::endl;
    // TODO: Check with mona on what the Vector3d coordinates of
    // this function are, and if we need to change them
    auto vis = leaf->getLeafVis(i);
    // We don't normally split normals, but in this case we have a flat surface
    // and we want to have a smooth shading
    //std::cout << "Inserting some geometry" << std::endl;
    //std::cout << "Vis has length " << vis.size() << std::endl;
    geometryTextureCoordinates.insert(geometryTextureCoordinates.begin() + (points_index/3*2),
                                            {currentLength, 0.0, currentLength, 1.0});
    //std::cout << "geometryNodeIds[" << points_index << "] = " << id << std::endl;
    //std::cout << "Vector Data: Size=" << geometryNodeIds.size() << ", capacity=" << geometryNodeIds.capacity() << std::endl;
    geometryNodeIds[points_index] = id;
    geometryNodeIds[points_index + 1] = id;
    // TODO: it not obvious that points_index can be changed by the insert here
    vec2Buf(geometry, points_index, vis[0], vis[1]);
    points_index = vec2Buf(geometryNormals, points_index, rot.Up(), rot.Up());
    geometryNodeIds[points_index/3] = id;
    geometryNodeIds[(points_index/3) + 1] = id;
    vec2Buf(geometry, points_index, vis[0], vis[1]);
    points_index = vec2Buf(geometryNormals, points_index ,-rot.Up(), -rot.Up());
    // The triangles are defined clockwise for the front face and counter clockwise for the back face
		
		unsigned int point_index_offset = points_index / 3;
      //std::cout << "Inserting some indices: " << geometryIndices.size() << " + 6 < " << geometryIndices.capacity() << std::endl;
		if(i > first_surface_id)
		{
			geometryIndices[cell_index +  0] = point_index_offset-7;
			geometryIndices[cell_index +  1] = point_index_offset-8;
			geometryIndices[cell_index +  2] = point_index_offset-4;
			geometryIndices[cell_index +  3] = point_index_offset-7;
			geometryIndices[cell_index +  4] = point_index_offset-4;
			geometryIndices[cell_index +  5] = point_index_offset-3;
			geometryIndices[cell_index +  6] = point_index_offset-5;
			geometryIndices[cell_index +  7] = point_index_offset-6;
			geometryIndices[cell_index +  8] = point_index_offset-2;
			geometryIndices[cell_index +  9] = point_index_offset-5;
			geometryIndices[cell_index + 10] = point_index_offset-2;
			geometryIndices[cell_index + 11] = point_index_offset-1;
			cell_index += 12;
		}
		//std::cout << "Inserted" << std::endl;
  }
      //std::cout << "Done" << std::endl;
}

void MappedPlant::GenerateStemGeometry(std::shared_ptr<Organ> stem, unsigned int p_o, unsigned int c_o)
{
  //std::cout << "Generating Stem for " << stem->getId() << " and reserving buffers" << std::endl;
	geometry.resize(std::max(static_cast<std::size_t>(p_o + (stem->getNumberOfNodes() * geometry_resolution * 3)), geometry.size()),-1.0);
	geometryNormals.resize(std::max(static_cast<std::size_t>(p_o + (stem->getNumberOfNodes() * geometry_resolution * 3)), geometryNormals.size()),-1.0);
  //std::cout << "Wanting to resize " << geometryIndices.size() << "/" << geometryIndices.capacity() << " indices to " << static_cast<std::size_t>(c_o + (stem->getNumberOfNodes() -1) * geometry_resolution * 6) << std::endl;
  //std::cout << "c_o=" << c_o << ", stem->getNumberOfNodes()=" << stem->getNumberOfNodes() << ", geometry_resolution=" << geometry_resolution << ", total is " << (stem->getNumberOfNodes() -1) * geometry_resolution * 6 << std::endl;
  geometryIndices.resize(std::max(static_cast<std::size_t>(c_o + (stem->getNumberOfNodes() -1) * geometry_resolution * 6), geometryIndices.size()),static_cast<unsigned int>(-1));
  geometryColors.resize(std::max(static_cast<std::size_t>(p_o + stem->getNumberOfNodes() * geometry_resolution), geometryColors.size()),static_cast<unsigned short>(-1));
  geometryTextureCoordinates.resize(std::max(static_cast<std::size_t>((p_o/3*2) + stem->getNumberOfNodes() * geometry_resolution * 2), geometryTextureCoordinates.size()),-1.0);
  geometryNodeIds.resize(std::max(static_cast<std::size_t>((p_o/3) + stem->getNumberOfNodes() * geometry_resolution), geometryNodeIds.size()),-1);
  Quaternion heading = Quaternion::FromMatrix3d(stem->iHeading);
  Quaternion lastRotation = Quaternion::geodesicRotation({1,0,0},{0,0,1});
  for(auto i = 0; i < stem->getNumberOfNodes(); ++i)
  {
    double diameter = stem->getParameter("radius");
    const auto& node = stem->getNode(i);

		// if the current i is in the last 10% of the nodes
		if(static_cast<float>(i)/static_cast<float>(stem->getNumberOfNodes()) > 0.9f)
		{
		  // reduce the diameter to form a tip
			diameter *= (1.0 - ((static_cast<float>(i)/static_cast<float>(stem->getNumberOfNodes()) - 0.9f) * 10.0f));
		}
		if(i == stem->getNumberOfNodes() - 1)
		{
			diameter = 0.01;
		}

    Vector3d dist;
    if(i + 1 < stem->getNumberOfNodes())
    {
      dist = stem->getNode(i + 1) - node;
    }
    else
    {
      dist = node - stem->getNode(i - 1);
    }
    
		lastRotation = Quaternion::FromForwardAndUp(dist, lastRotation.Up());
    auto deltaphi = 2 * M_PI / geometry_resolution;
		unsigned int point_index_offset = p_o / 3;
    for (auto j = 0; j < geometry_resolution; ++j)
    {
      auto phi = j * deltaphi;
      Vector3d outer = {0.0, std::cos(phi) * diameter, std::sin(phi) * diameter};
      outer = lastRotation.Rotate(outer);
			//std::cout << "[Leaf] Rotating " << lastRotation.toString() << " to get " << outer.toString() << std::endl;
      // consequtive points are stored in the buffer
      // the plus 0 is not necessary, but it makes it easier to read
      geometry[3 * (i * geometry_resolution + j) + 0 + p_o] = node.x + outer.x;
      geometry[3 * (i * geometry_resolution + j) + 1 + p_o] = node.y + outer.y;
      geometry[3 * (i * geometry_resolution + j) + 2 + p_o] = node.z + outer.z;
			// calculating the point index offset from the buffer offset
      // the indices are stored in the buffer, and are all front facing
			if (i > 0)
			{
				geometryIndices[6 * ((i-1) * geometry_resolution + j) + 2 + c_o] = point_index_offset + ((i-1) + 1) * geometry_resolution + (j + 1) % geometry_resolution;
				geometryIndices[6 * ((i-1) * geometry_resolution + j) + 1 + c_o] = point_index_offset +  (i-1) * geometry_resolution + (j + 1) % geometry_resolution;
				geometryIndices[6 * ((i-1) * geometry_resolution + j) + 0 + c_o] = point_index_offset +  (i-1) * geometry_resolution + j;
				geometryIndices[6 * ((i-1) * geometry_resolution + j) + 3 + c_o] = point_index_offset +  (i-1) * geometry_resolution + j;
				geometryIndices[6 * ((i-1) * geometry_resolution + j) + 4 + c_o] = point_index_offset + ((i-1) + 1) * geometry_resolution + (j + 1) % geometry_resolution;
				geometryIndices[6 * ((i-1) * geometry_resolution + j) + 5 + c_o] = point_index_offset + ((i-1) + 1) * geometry_resolution + j;
			}
      // the normals are stored in the buffer
      geometryNormals[3 * (i * geometry_resolution + j) + 0 + p_o] = outer.x;
      geometryNormals[3 * (i * geometry_resolution + j) + 1 + p_o] = outer.y;
      geometryNormals[3 * (i * geometry_resolution + j) + 2 + p_o] = outer.z;

      geometryNodeIds[i * geometry_resolution + j + point_index_offset] = stem->getNodeId(i);

      geometryTextureCoordinates[2 * (i * geometry_resolution + j) + 0 + point_index_offset*2] = i / (double)stem->getNumberOfNodes();
      geometryTextureCoordinates[2 * (i * geometry_resolution + j) + 1 + point_index_offset*2] = phi / (2.0 * M_PI);
    }
  }
}

void MappedPlant::GenerateRadialLeafGeometry(std::shared_ptr<Leaf> leaf, unsigned int p_o, unsigned int c_o)
{
	// Fetch the phi array
	double scaling_factor = leaf->getParameter("areaMax") * leaf->getLength(false) / leaf->getParameter("k");

	// resolution
	int resolution = leaf_resolution;
	// Compute the mid vein of the leaf
	CatmullRomSplineManager midVein = leaf->getNodes();
	// Compute the leaf length
	auto length = leaf->getLength(false);
	// save c_o
	unsigned int start_c_o = c_o;
  // std::cout << "Accessing leaf random parameter for leaf " << leaf->getId() << std::endl;
	// get leaf random parameter
	auto lrp = leaf->getLeafRandomParameter();
	auto stem = leaf->getParent();
	auto min_radius = stem->getParameter("radius");

  // std::cout << "Invoking create leaf radial geometry for leaf " << leaf->getId() << std::endl;

	// create leaf radial geometry
	// greate points for the leaf outer
	// double check that this was not done in python
	if(lrp->leafGeometry.size() == 0)
	{
		lrp->createLeafRadialGeometry(lrp->leafGeometryPhi,lrp->leafGeometryX,resolution);
	}
	// retrieve the leaf geometry
	const auto& outer_geometry_points = lrp->leafGeometry;
  // set buffer sizes
  int point_buffer = 0;
  int index_buffer = 0;
  int last_amount = -1;
  int last_non_petiole = -1;
  double r_max = std::numeric_limits<float>::lowest();
  std::cout << "Counting how much space we need for the leaf geometry" << std::endl;
  for (auto i = 0; i < outer_geometry_points.size(); ++i)
  {
		MirrorIterator helper(&(outer_geometry_points[i]));
		auto current_amount = helper.size();
		r_max = std::max(r_max, *std::max_element(outer_geometry_points[i].begin(), outer_geometry_points[i].end()));
    if(current_amount < 2)
    {
      //std::cout << "Skipping petiole at " << i << " because it has size " << current_amount << std::endl;
			last_non_petiole = -1;
      continue;
    }
    //std::cout << "NOT Skipping petiole at " << i << " because it has size " << current_amount << std::endl;
    point_buffer += current_amount;
    if(last_non_petiole > -1 && i > last_non_petiole)
    {
				//std::cout << "Adding Ts : " << last_amount << "/" << current_amount << std::endl;
			if(last_amount != current_amount)
			{
				index_buffer += (std::min(last_amount, (int)current_amount) - 1) * 6;
				index_buffer += (std::abs(last_amount - (int)current_amount) - 1) * 3;
				point_buffer ++;
			}
			else
			{
				index_buffer += (current_amount - 1) * 6;
			}
    }
    if(current_amount > 1)
    {
      last_non_petiole = i;
			last_amount = current_amount;
    }
    if(i > last_non_petiole)
    {
      last_amount = current_amount;
    }
  }
  std::cout << "Resizing geometry buffers by " << point_buffer << " points and " << index_buffer << " triangle values." << std::endl;
  // increase geometry buffers
  this->geometry.resize(std::max(static_cast<std::size_t>(p_o + point_buffer * 3), this->geometry.size()),-1.0);
	this->geometryIndices.resize(std::max(static_cast<std::size_t>(c_o + index_buffer), this->geometryIndices.size()),static_cast<unsigned int>(-1));
	this->geometryNormals.resize(std::max(static_cast<std::size_t>(p_o + point_buffer * 3), this->geometryNormals.size()),-1.0);
	this->geometryTextureCoordinates.resize(std::max(static_cast<std::size_t>((p_o/3*2) + point_buffer * 2), this->geometryTextureCoordinates.size()),-1.0);
	this->geometryNodeIds.resize(std::max(static_cast<std::size_t>(p_o / 3 + point_buffer), this->geometryNodeIds.size()),-1);
	// get the number of points
  // std::cout << "Iterating through the line intersections and generating the geometry" << std::endl;
  last_amount = -1;
  last_non_petiole = -1;
	Quaternion last_orientation;
	Vector3d last_position;
	int last_index{-1};
  // create two random factors between 0 and 1
  float random_factor_1 = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
  float random_factor_2 = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);



	for(auto i = 0; i < outer_geometry_points.size(); ++i)
	{
    const std::vector<double>& current = outer_geometry_points[i];
		MirrorIterator helper(&current);
		auto current_amount = helper.size();
    if(current_amount < 2)
    {
      last_non_petiole = -1;
      continue;
    }
		// compute the size of the current array, which is double its size unless one point is near zero which is only counted once
		int current_size = current_amount;
		// get the current point
    double t = static_cast<double>(i) / static_cast<double>(resolution);
		double l = t * length;
    auto midpoint = (i == 0) ? leaf->getNode(0) : midVein(t);
    // get the current point
		// get the best spline for the current point
		auto select_spline = midVein.selectSpline(t);
		//auto lowest_possible_spline = std::find(midVein.getSplines().begin(), midVein.getSplines().end(), [l](auto spline) {return spline.getT0() < l; });
		
	  // input points, normaly, ids, texture coordinates
    // iterate through the points
		Quaternion local_q = select_spline.computeOrientation(t);
		auto up = local_q.Up();
    // iterate through the points
    //std::cout << "Iterating through the points of the current line intersection " << i << std::endl;

		float petiole_distance = leaf->getLeafRandomParameter()->lb;

		const Vector3d first_node = leaf->getNode(0);
		const Vector3d first_estimated = midVein(0);
		const auto first_distance = (first_node - first_estimated).length();
		std::cout << "First distance is " << first_distance << std::endl;
    
    for(int p = 0; p < helper.size(); ++p)
    {
      //std::cout << p_o << "/" << geometry.size() << " ";

      auto r = helper[p];
      // get the point
			// get the wave effect which is a sine function along the length of the leaf

			float z_offset =  0.33 * (helper.isMirrored(p) ? -1.0 : 1.0);
			// make two different sine waves for each side
			if(helper.isMirrored(p))
			{
			  z_offset *= std::sin((2.0*random_factor_1 + 2.0) * l / M_PI);
      }
      else
      {
        z_offset *= std::sin((2.0 * random_factor_2 + 2.0) * l / M_PI + M_PI);
      }
			const Vector3d base_direction = r * Vector3d(0.0, 1.0, 0.0);
			//std::cout << base_direction.toString() << std::endl;
			Vector3d updated_direction = local_q.Rotate(base_direction);
			updated_direction = (updated_direction.length() > min_radius) ? updated_direction : min_radius * updated_direction.normalized();
			const Vector3d point = midpoint + updated_direction + up * z_offset;
      //std::cout << "V: " << point.toString() << "; ";
      // set the point
      //std::cout << "p" << " ";
      geometry[p_o + 0] = point.x;
      geometry[p_o + 1] = point.y;
      geometry[p_o + 2] = point.z;
      // set the normal
      //std::cout << "n" << " ";
      geometryNormals[p_o + 0] = up.x;
      geometryNormals[p_o + 1] = up.y;
      geometryNormals[p_o + 2] = up.z;
      // set the texture coordinates
      //std::cout << "t" << " ";
      geometryTextureCoordinates[(p_o/3*2)] = l;
      geometryTextureCoordinates[(p_o/3*2) + 1] = helper.texcoord(p);
			// set the node id
      //std::cout << "i" << " ";
			geometryNodeIds[p_o/3] = 1;
			// increase buffer
			p_o += 3;
    }
    std::cout << std::endl;
    //std::cout << std::endl << "Generating the triangles for the current line intersection " << i << "(" << current.size() << ")" << std::endl;
    if(i > last_non_petiole && last_non_petiole >= 0)
    {
      // use the case distinction between number of intersections
      // for the triangulation between connected sections of the surface
      if(current_amount == last_amount)
      {
        // std::cout << " which is equal to the last one" << std::endl;
        // we construct pairwise triangles for the two sections
        for(auto j = 0; j < current_amount - 1; j += 2)
        {
					// first triangle
					geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount;
					geometryIndices[c_o++] = (p_o/3) - current_amount + j + 1;
					geometryIndices[c_o++] = (p_o/3) - current_amount + j;
					// second triangle
					geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount;
					geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount + 1;
					geometryIndices[c_o++] = (p_o/3) - current_amount + j + 1;
        }
      }
      else
      {
        // set the normal
        geometryNormals[p_o + 0] = up.x;
        geometryNormals[p_o + 1] = up.y;
        geometryNormals[p_o + 2] = up.z;
        // set the texture coordinates
        geometryTextureCoordinates[(p_o/3*2)] = t;
        geometryTextureCoordinates[(p_o/3*2) + 0] = 0.0;
        // set the node id
        geometryNodeIds[p_o/3] = 1.0;
        if(current_amount > last_amount)
        {
					// since we have more points in one of the sections
					// we have to construct triangles with the midvein in mind
					// we construct pairwise triangles for the two sections
					geometry[p_o + 0] = last_position.x;
					geometry[p_o + 1] = last_position.y;
					geometry[p_o + 2] = last_position.z;
          // std::cout << " which is larger to the last one" << std::endl;
          // set the triangles before we increase the buffer to keep the indices correct
          // outer triangles are the only ones connected to the last outline points
          // compute difference between the number of points
          auto diff = current_amount - last_amount;
          // iterate through the top points until we reach the midvein plus the difference
          // std::cout << "iterate through the top points until we reach the midvein plus the difference" << std::endl;
          for(auto j = 0; j < current_amount/2 -diff - 1; ++j)
          {
						// first triangle
						geometryIndices[c_o++] = (p_o/3) - current_amount + j;
						geometryIndices[c_o++] = (p_o/3) - current_amount + j + 1;
						geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount;
						// second triangle
						geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount;
						geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount + 1;
						geometryIndices[c_o++] = (p_o/3) - current_amount + j + 1;
          }
          // iterate through the bottom points starting from the midvein plus the difference
          // std::cout << "iterate through the bottom points starting from the midvein plus the difference" << std::endl;
          for(auto j = current_amount/2+diff; j < current_amount - 1; ++j)
          {
						// first triangle
						geometryIndices[c_o++] = (p_o/3) - current_amount + j;
						geometryIndices[c_o++] = (p_o/3) - current_amount + j + 1;
						geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount - diff;
						// second triangle
						geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount - diff;
						geometryIndices[c_o++] = (p_o/3) - current_amount + j - last_amount + 1 - diff;
						geometryIndices[c_o++] = (p_o/3) - current_amount + j + 1;
          }
          // iterate through the midvein points
          // std::cout << "iterate through the midvein points" << std::endl;
          // note that this is a specific interpretation of what happens here and might not be correct
          // for a complete correct implementation of this, we would have to start with the phi array
          for(auto j = 0; j < diff - 1; ++j)
          {
            // we triangulate from one point against pairs of points
            // so we only need one triangle, but with the most recent point included
            geometryIndices[c_o++] = (p_o/3) - current_amount/2 + j;
            geometryIndices[c_o++] = (p_o/3) - current_amount/2 + j + 1;
            geometryIndices[c_o++] = (p_o/3);
          }
          // increase buffer
          p_o += 3;
        }
        else
        {
					// since we have more points in one of the sections
					// we have to construct triangles with the midvein in mind
					// we construct pairwise triangles for the two sections
					geometry[p_o + 0] = midpoint.x;
					geometry[p_o + 1] = midpoint.y;
					geometry[p_o + 2] = midpoint.z;
          std::cout << " which is smaller to the last one" << std::endl;
          auto diff = last_amount - current_amount;
          for(auto j = 0; j < last_amount/2 -diff - 1; ++j)
          {
            geometryIndices[c_o++] = (p_o/3) - current_amount - last_amount + j;
            geometryIndices[c_o++] = (p_o/3) - current_amount + j; 
            geometryIndices[c_o++] = (p_o/3) - current_amount - last_amount + j + 1;
            geometryIndices[c_o++] = (p_o/3) - current_amount - last_amount + j + 1;
            geometryIndices[c_o++] = (p_o/3) - current_amount + j;
            geometryIndices[c_o++] = (p_o/3) - current_amount + j + 1;
          }
          for(auto j = last_amount/2+diff; j < last_amount - 1; ++j)
          {
            geometryIndices[c_o++] = (p_o/3) - current_amount - last_amount + j;
            geometryIndices[c_o++] = (p_o/3) - current_amount + j - diff;; 
            geometryIndices[c_o++] = (p_o/3) - current_amount - last_amount + j + 1; 
            geometryIndices[c_o++] = (p_o/3) - current_amount + j - diff;
            geometryIndices[c_o++] = (p_o/3) - current_amount + j - diff + 1; 
            geometryIndices[c_o++] = (p_o/3) - current_amount - last_amount + j + 1;
          }
          for(auto j = 0; j < diff - 1; ++j)
          {
            geometryIndices[c_o++] = (p_o/3) - current_amount - last_amount/2 + j;
            geometryIndices[c_o++] = (p_o/3); 
            geometryIndices[c_o++] = (p_o/3) - current_amount - last_amount/2 + j + 1;
          }
          p_o += 3;
        }
      }
    }
    if(current_amount > 1)
    {
      last_non_petiole = i;
			last_amount = current_amount;
		}
		last_orientation = local_q;
		last_position = midpoint;
	}
	std::cout << "In the end I ended up adding " << c_o - start_c_o << " where I thought I'd add " << index_buffer << std::endl;
}

void MappedPlant::GenerateRadialLeafGeometryFromPhi(std::shared_ptr<Leaf> leaf, unsigned int p_o, unsigned int c_o)
{
}

} // namespace
