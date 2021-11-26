import sys; sys.path.append("../src/python_modules/"); sys.path.append("../")

import plantbox as pb
import rsml_reader
import xylem_flux

import numpy as np


class DataModel:
    """ 
    Model in the sense of model view controller (MVC), stores most that is presented in the view
    
    * manages the reading of rsml DataModel.open_rsml
    * can add an artifical shoot, see DataModel.add_artificial_shoot
    
    Rest of the methods should not be called directly          
    """

    def __init__(self):
        self.fname = None  # file name (str)
        self.polylines = None  # from rsml (list of list of list)
        self.properties = None  # from rsml (dict of list) value per root
        self.functions = None  # from rsml (dict of list of list) value per node
        self.metadata = None  # from rsml (rsml_writer.Metadata)
        self.radii, self.cts, self.types, self.tagnames = None, None, None, None  # selected from rsml
        self.analyser = None  # created by convert_to_xylem_flux_
        self.xylem_flux = None  # created by convert_to_xylem_flux_
        self.max_ct = 0.  # created by convert_to_xylem_flux_
        self.base_nodes = [0]  # base nodes indices (of roots or multiple plants)
        self.base_segs = [0]  # emerging segment indices from base nodes

        # add more where needed: key = cm *value (for length), and key = day * value (for time)
        self.scales_ = {"pixel":1, "px":1, "dots": 1,
                   "cm": 1, "mm": 0.1, "dm": 10., "m": 100,
                   "h": 1. / 24., "s": 1 / (24.*3600), "sec": 1 / (24.*3600), "d": 1., "day": 1, "days": 1}

    def exists(self):
        """ true if a rsml was set """
        return self.polylines is not None

    def set_rsml(self, polylines, properties, functions, metadata):
        """ setter for rsml fields """
        self.polylines = polylines
        self.properties = properties
        self.functions = functions
        self.metadata = metadata

    def set_selected(self, radii, cts, types, tagnames):
        """ setter for selected tags """
        self.radii = radii
        self.cts = cts
        self.types = types
        self.tagnames = tagnames

    def open_rsml(self, fname):
        """ 
        opens an rsml file into self.data, using rsml_reader (in CPlantBox/src/python_modules)  
             
        * converts units to cm and day 
        * if necessary converts 2d -> 3d, 
        * creates an analyser (pb.SegmentAnalyser) and a xylem_flux (pb.XylemFluxPython) object 
        """
        polylines, properties, functions, metadata = rsml_reader.read_rsml(fname)
        print("DataModel.open_rsml(): scale to cm", metadata.scale_to_cm)
        self.set_rsml(polylines, properties, functions, metadata)
        self.scale_polylines_()  # converts units
        self.check_polylines_2d_()  # 2d -> 3d
        radii, cts, types, tagnames = rsml_reader.get_parameter(polylines, functions, properties)  # paramter per node
        self.set_selected(radii, cts, types, tagnames)
        self.scale_selected_()  # converts units of special fields radii, cts, types
        self.convert_to_xylem_flux_()

    def scale_polylines_(self):
        """ 
        scales nodes, see rsml_writer.Metadata, and self.scale_to_cm 
        """
        scale = self.metadata.scale_to_cm  # default length scales
        for i in range(0, len(self.polylines)):
            for j in range(0, len(self.polylines[i])):
                for k in range(0, 3):
                    self.polylines[i][j][k] *= scale

    def check_polylines_2d_(self):
        """ 
        converts 2d image coordinates to 3d coordinates
        """
        nodes, segs = rsml_reader.get_segments(self.polylines, self.properties)  # fetch nodes and segments
        maxz = np.max(nodes[:, 2])
        minz = np.min(nodes[:, 2])
        if maxz >= 0 and minz >= 0:  # image coordinates in px often start in lower left corner
            print("DataModel.check_polylines_2d_(): assuming image coordinates, y-centered and z-flipped ")
            miny = np.min(nodes[:, 1])
            yy = np.max(nodes[:, 1]) - miny
            for pl in self.polylines:  # both are references
                for node in pl:
                    node[2] = -node[2]
                    node[1] = node[1] - miny - yy / 2

    def scale_selected_(self):
        """ 
        scales radius and creation times, see rsml_writer.Metadata, and self.scale_to_cm
        shifts types, in a way they start at zero
        """
        scale = self.metadata.scale_to_cm  # default length scales
        # radii
        extra_str = ""
        if self.tagnames[0]:
            if self.tagnames[0] in self.metadata.properties:
                r_scale = self.scales_[self.metadata.properties[self.tagnames[0]].unit]
            else:  # assume same scaling as polylines
                r_scale = scale
        else:  # assume same scaling as polylines
            r_scale = scale
        if self.metadata.software == "smartroot":
            r_scale = scale
            extra_str = " (smartroot)"
        print("DataModel.scale_selected_():radius length scale" + extra_str, r_scale)
        for i in range (0, len(self.radii)):
            self.radii[i] *= r_scale
        # creation times
        cts_scale = 1.  # assume it is in days
        if self.tagnames[1]:
            if self.tagnames[1] in self.metadata.properties:
                cts_scale = self.scales_[self.metadata.properties[self.tagnames[1]].unit]
                print("DataModel.scale_rsml() temporal scale", cts_scale)
        for i in range (0, len(self.cts)):
            self.cts[i] *= cts_scale
        # types
        min_types = np.min(self.types)
        for i in range (0, len(self.types)):
            self.types[i] -= min_types

    def convert_to_xylem_flux_(self):
        """ 
        converts the polylines to a SegmentAnalyser and a MappedSegments object, and stores max_ct   
        
        uses:
        properties["parent-poly"], properties["parent-nodes"]
        radii, cts, types                    
        """
        nodes, segs = rsml_reader.get_segments(self.polylines, self.properties)  # fetch nodes and segments
        segRadii = np.zeros((segs.shape[0], 1))  # convert to paramter per segment
        segCTs = np.zeros((segs.shape[0], 1))
        subTypes = np.zeros((segs.shape[0], 1))
        if np.isnan(self.cts[0]):  # nan indicates creation times not given in rsml
            self.cts = np.zeros((len(self.cts),))
        for i, s in enumerate(segs):
            segRadii[i] = self.radii[s[1]]  # seg to node index
            segCTs[i] = self.cts[s[1]]
            subTypes[i] = self.types[s[1]]
        if np.isnan(subTypes[0]):
            subTypes = np.ones((len(segs),), dtype = np.int64)
        self.max_ct = np.max(segCTs)
        segs_ = [pb.Vector2i(s[0], s[1]) for s in segs]  # convert to CPlantBox types
        nodes_ = [pb.Vector3d(n[0], n[1], n[2]) for n in nodes]
        self.analyser = pb.SegmentAnalyser(nodes_, segs_, segCTs, segRadii)
        self.analyser.addData("subType", subTypes)
        ms = pb.MappedSegments(self.analyser.nodes, np.array(self.cts), segs_, np.array(segRadii), np.array(subTypes))
        self.xylem_flux = xylem_flux.XylemFluxPython(ms)
        self.base_nodes = self.get_base_node_indices_()
        self.xylem_flux.neumann_ind = self.base_nodes  # needed for suf
        self.xylem_flux.dirichlet_ind = self.base_nodes  # needed for krs
        self.base_segs = self.xylem_flux.find_base_segments()

    def get_base_node_indices_(self):
        """
        get all node indices of base roots
        """
        c = 0
        bni = []
        for i, p in enumerate(self.polylines):
            if self.properties['parent-poly'][i] == -1:
                bni.append(c)
            c += len(p)  # next first node index
        return bni

    def add_artificial_shoot(self):
        """
        adds a 1 cm shoot element, connecting all base roots
        the type for looking up conductivities is set to 10 
        """
        nodes = self.analyser.nodes
        bni = self.base_nodes
        print("DataModel.add_artificial_shoot() base node indices are", bni)
        mid = np.zeros((3,))
        for i in bni:
            mid += np.array([nodes[i].x, nodes[i].y, nodes[i].z])
        mid /= len(bni)
        print("DataModel.add_artificial_shoot() mid point is", mid)

        rsml_reader.artificial_shoot(self.polylines, self.properties, self.functions)  # append artifial shoot (default values)
        radii, cts, types, tagnames = rsml_reader.get_parameter(self.polylines, self.functions, self.properties)  # paramter per node
        # change default values from artificial shoot
        collar = mid.copy()
        mid[2] += 0.1  # shift a bit up
        collar[2] += 1.1  # 1 cm shoot length
        self.polylines[0][0] = collar
        self.polylines[0][1] = mid
        self.set_selected(radii, cts, types, tagnames)
        self.scale_selected_()
        self.radii[0] = 0.1  # cm
        self.radii[1] = 0.1  # cm
        self.types[0] = 10
        self.types[1] = 10
        # print("after ADD SHOOT")
        # print(self.polylines)
        # print(self.properties["parent-poly"])
        # print(self.properties["parent-node"])
        # print("radii", radii)
        # print("cts", cts)
        # print("types", types)
        self.convert_to_xylem_flux_()
#         print("mid ", str(self.analyser.nodes[1]), " cm")
#         print("collar ", str(self.analyser.nodes[0]), " cm")
#         print("seg 0", str(self.analyser.segments[0]))
#         print("radius", str(self.analyser.data["radius"][0]))
#         print("type", str(self.analyser.data["subType"][0]))
#         print("seg 1", str(self.analyser.segments[1]))
#         print("radius", str(self.analyser.data["radius"][1]))
#         print("type", str(self.analyser.data["subType"][1]))

    def add_creation_times(self):
        """
        lineary interpolates creation times assuming a lateral delay time of one day
        """
        pl_ = self.polylines  # rename
        self.functions["creation_time"] = [None] * len(pl_)
        for i, pl in enumerate(pl_):
            if self.properties["parent-poly"][i] == -1:
                self.functions["creation_time"][i] = np.zeros((len(pl,)))
                if not (i == 0 and len(pl) == 2):  # not artifical shoot, else [0,0] is fine
                    ct = self.functions["creation_time"][i]  # rename
                    lt = self.get_length_(pl)
                    l = 0
                    ct[0] = 0
                    for j in range(0, len(pl) - 1):
                        l += np.linalg.norm(np.array(pl[j + 1]) - np.array(pl[j]))
                        ct[j + 1] = self.max_ct * (l / lt)
        for i, pl in enumerate(pl_):
            if self.functions["creation_time"][i] is None:
                self.add_creation_times_(i)
        radii, cts, types, tagnames = rsml_reader.get_parameter(self.polylines, self.functions, self.properties)  # paramter per node
        self.set_selected(radii, cts, types, tagnames)
        self.scale_selected_()
        self.convert_to_xylem_flux_()

    def add_creation_times_(self, i):
        """ 
        recursive funtion called by add_creation_times, 
        interpolates a single polyline assuming a lateral delay time of one day
        """
        ppi = self.properties["parent-poly"][i]
        if self.functions["creation_time"][ppi] is None:
            print(ppi, self.functions["creation_time"][ppi])
            self.add_creation_times_(ppi)  # label parent first
        pl = self.polylines[i]
        self.functions["creation_time"][i] = np.zeros((len(pl),))
        ct = self.functions["creation_time"][i]  # rename
        lt = self.get_length_(pl)
        pni = self.properties["parent-node"][i]
        # print("Root", i, "Parent", ppi, "node index", pni)
        ct[0] = self.functions["creation_time"][ppi][pni]  #  + 1  # 1 day
        l = 0.
        for j in range(0, len(pl) - 1):
            l += np.linalg.norm(np.array(pl[j + 1]) - np.array(pl[j]))
            ct[j + 1] = ct[0] + (self.max_ct - ct[0]) * (l / lt)

    def get_length_(self, pl):
        """ polyline length """
        lt = 0
        for j in range(0, len(pl) - 1):
            lt += np.linalg.norm(np.array(pl[j + 1]) - np.array(pl[j]))
        return lt
