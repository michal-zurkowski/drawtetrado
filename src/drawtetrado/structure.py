import json
import sys
import math
import subprocess

from drawtetrado.svg_painter import Point, ConnType, ConnFlow

class Nucleotide:
    def FindConnections(self, used):
        # Find to what it is connected.
        connected_to = ""
        curr_min = sys.maxsize
        for name, nucl in used.items():
            if nucl["index"] > self.index and \
               curr_min > nucl["index"] and \
               self.chain == nucl["chain"]:
                curr_min = nucl["index"]
                connected_to = name
        return connected_to

    def Block(self, width, height, angle):
        self.coords = []
        self.coords.append(Point(0.0, 0.0))
        shift = height / math.tan(math.radians(angle))
        self.coords.append(Point(shift, -height))
        self.coords.append(Point(shift + width, -height))
        self.coords.append(Point(width, 0.0))
        self.center = Point((shift + width) / 2.0, -height / 2.0)

    def ShiftBlock(self, shift):
        self.coords[0] += shift
        self.coords[1] += shift
        self.coords[2] += shift
        self.coords[3] += shift
        self.center += shift

    def CalculateCoordinates(self, conf):
        shift = Point(0, 0)
        tan_val = math.tan(math.radians(conf.angle))
        sin_val = math.sin(math.radians(conf.angle))
        cos_val = math.cos(math.radians(conf.angle))
        tan_val_rev = math.tan(math.radians(90.0 - conf.angle))

        shift.y -= self.tetrade_no * (sin_val * (conf.longer + conf.shorter + conf.spacing) + \
                                      conf.tetrade_spacing)
        if self.onz == "-":
            width_0 = conf.longer
            height_0 = sin_val * conf.shorter
            width_1 = conf.shorter
            height_1 = sin_val * conf.longer
        else: # Should be just this: self.onz == "?+": 
            width_0 = conf.shorter
            height_0 = sin_val * conf.longer
            width_1 = conf.longer
            height_1 = sin_val * conf.shorter

        if self.position == 0:
            self.Block(width_0, height_0, conf.angle)
        elif self.position == 1:
            self.Block(width_1, height_1, conf.angle)
            shift.y -= (width_1 + conf.spacing) * sin_val
            shift.x += (width_1 + conf.spacing) * cos_val
        elif self.position == 2:
            self.Block(width_0, height_0, conf.angle)
            shift.y -= (width_0 + conf.spacing) * sin_val
            shift.x += (width_1 + conf.spacing) + \
                       ((width_0 + conf.spacing) * sin_val) * tan_val_rev
        elif self.position == 3:
            self.Block(width_1, height_1, conf.angle)
            shift.x += width_0 + conf.spacing

        self.ShiftBlock(shift)

    def GetOnz(self):
        if self.onz_full[-1] == "-":
            return "-"
        return "+"

    def __init__(self, data, used_nucl, tetr_no, tetr_onz, pos):
        self.number = data["number"]
        self.short_name = data["shortName"]
        self.full_name = data["fullName"]
        self.chain = data["chain"]
        self.index = data["index"]
        self.onz_full = tetr_onz
        self.onz = self.GetOnz()
        self.bond = data["glycosidicBond"]
        self.tetrade_no = tetr_no
        self.position = pos
        self.connected_to = self.FindConnections(used_nucl)
        self.connected_from = ""
        self.coords = []
        self.center = Point(0, 0)
        self.connection_type = ConnType.UNKNOWN

        # Parameters regarding connection from this nucleotide to another.
        # out - flow starting the connection
        # in - flow ending the conneciton
        self.flow_out = ConnFlow.UNKNOWN
        self.flow_in = ConnFlow.UNKNOWN


        # Priority of drawing.
        self.priority_conn = -1
        self.priority_edge = -1
        self.priority_nucl = -1

    def UpdatePriorities(self, nucleotides):
        # 1
        if self.connection_type == ConnType.LEFT or \
           self.connection_type == ConnType.LEFT_CROSS or \
           self.connection_type == ConnType.FRONT_TO_BACK:
            self.priority_conn = 1
        elif self.connection_type == ConnType.FRONT_BACK_CROSS and \
             (self.position == 1 or self.position == 2):
            self.priority_conn = 1
        elif self.connection_type == ConnType.SAME_LEVEL and \
             (self.flow_out == ConnFlow.DOWN and self.flow_in == ConnFlow.DOWN):
            self.priority_conn = 1
        elif self.connection_type == ConnType.SIMPLE and self.position == 1:
            self.priority_conn = 1
        # 2 + 3
        if self.connection_type == ConnType.SAME_LEVEL and \
             (self.flow_out == ConnFlow.UP and self.flow_in == ConnFlow.DOWN):
            self.priority_nucl = 2
            self.priority_conn = 3
        elif self.connection_type == ConnType.SAME_LEVEL and \
             (self.flow_out == ConnFlow.DOWN and self.flow_in == ConnFlow.UP):
            if self.connected_to != "":
                nucleotides[self.connected_to].priority_nucl = 2
            self.priority_conn = 3
        # 4
        if self.priority_nucl == -1:
            self.priority_nucl = 4
        # 5 (Have to be later to be on top of the nucleotides)
        if self.priority_edge == -1:
            self.priority_edge = 5
        # 6 (all other connections are processed)
        if  self.priority_conn == -1:
            self.priority_conn = 6


class Quadruplex:
    def UsedNucleotides(self, tetrad, nucl):
        used = {}
        for name, data in tetrad.items():
            used[data["nt1"]] = nucl[data["nt1"]]
            used[data["nt2"]] = nucl[data["nt2"]]
            used[data["nt3"]] = nucl[data["nt3"]]
            used[data["nt4"]] = nucl[data["nt4"]]
        return used

    def PrepareNucleotides(self, structure, quadruplex_id, tetrad_id):
        nucl = structure.nucleotides
        if tetrad_id >= 0:
            tetrads = structure.single_tetrads[quadruplex_id][tetrad_id]
        else:
            tetrads = structure.tetrads[quadruplex_id]
        tetrads_order = structure.tetrads_order[quadruplex_id]

        used_nucl = self.UsedNucleotides(tetrads, nucl)

        tetr_no = 0
        for tetrad_name in tetrads_order:
            if tetrad_name in tetrads:
                tetrad = tetrads[tetrad_name]
            else:
                continue
            nt1 = tetrad["nt1"]
            nt2 = tetrad["nt2"]
            nt3 = tetrad["nt3"]
            nt4 = tetrad["nt4"]
            onz = tetrad["onz"]

            self.nucl_quad[nt1] = Nucleotide(nucl[nt1], used_nucl, tetr_no, onz, 0)
            self.nucl_quad[nt2] = Nucleotide(nucl[nt2], used_nucl, tetr_no, onz, 1)
            self.nucl_quad[nt3] = Nucleotide(nucl[nt3], used_nucl, tetr_no, onz, 2)
            self.nucl_quad[nt4] = Nucleotide(nucl[nt4], used_nucl, tetr_no, onz, 3)

            self.tetrads.append([nt1, nt2, nt3, nt4])

            tetr_no = tetr_no + 1


    def GetChainFirstLast(self):
        chains = {}
        for _, nucl in self.nucl_quad.items():
            if nucl.chain not in chains:
                chains[nucl.chain] = {"first": "",\
                                      "last": "", \
                                      "val": sys.maxsize }
        # Analyze 
        curr_name = ""
        curr_min = sys.maxsize
        for name, nucl in self.nucl_quad.items():
            if chains[nucl.chain]["val"] > nucl.index:
                chains[nucl.chain]["val"] = nucl.index
                chains[nucl.chain]["first"] = name

        # Full connected_from variables.
        for chain, data in chains.items():
            curr_name = data["first"]
            next_conn = self.nucl_quad[data["first"]].connected_to
            while next_conn != "":
                self.nucl_quad[next_conn].connected_from = curr_name
                curr_name = next_conn
                next_conn = self.nucl_quad[next_conn].connected_to
            data["last"] = curr_name

        return chains

    def __init__(self, structure, quadruplex_id, tetrad_id = -1):
        self.nucl_quad = {}
        self.tetrads = []
        self.PrepareNucleotides(structure, quadruplex_id, tetrad_id)
        self.chains = self.GetChainFirstLast()
        if tetrad_id >= 0:
            single_tracts = []
            single_tracts.append(structure.tracts[quadruplex_id][tetrad_id])
            self.tracts = single_tracts
        else:
            self.tracts = structure.tracts[quadruplex_id]

    def GetNucleotidesPositions(self):
        lst = [-1] * len(self.nucl_quad)
        for name, nucl in self.nucl_quad.items():
            if nucl.connected_to != "":
              conn_to = self.nucl_quad[nucl.connected_to]
              lst[nucl.tetrade_no * 4 + nucl.position] = conn_to.tetrade_no * 4 + conn_to.position
            else:
              lst[nucl.tetrade_no * 4 + nucl.position] = -1
        return lst

    def GetSameRotations(self):
        lst = [-1] * len(self.tetrads)
        for idx, tetrad in enumerate(self.tetrads):
            for name in tetrad:
                for group, tetrad_tracts in enumerate(self.tracts):
                    for names in tetrad_tracts:
                        for nucl_name in names:
                            if name == nucl_name:
                                lst[idx] = group
        return lst

    def GetAlignments(self):
        lst = [-1] * len(self.nucl_quad)
        for name, nucl in self.nucl_quad.items():
            for group_1, tetrad_tracts in enumerate(self.tracts):
                for group_2, names in enumerate(tetrad_tracts):
                    for nucl_name in names:
                        if name == nucl_name:
                            lst[nucl.tetrade_no * 4 + nucl.position] = group_1 * 4 + group_2
        return lst

    # Use C++ code to rotate tetrads for more readable output.
    def Optimize(self, optimizer = "./svg_optimizer"):
        import optimizer
        optimized = optimizer.solve(self.GetNucleotidesPositions(),
                                    self.GetSameRotations(),
                                    self.GetAlignments())

        # Update position for nucleotide.
        for _, nucl in self.nucl_quad.items():
            pos = nucl.position + nucl.tetrade_no * 4
            for x in range(4):
                if int(optimized[nucl.tetrade_no * 4 + x]) == nucl.position:
                    nucl.position = x
                    break

        # Update positions in tetrades. For Tetrade border
        level = 0
        for tetrad in self.tetrads:
            new_positions = []
            new_positions.append(tetrad[int(optimized[level * 4 + 0])])
            new_positions.append(tetrad[int(optimized[level * 4 + 1])])
            new_positions.append(tetrad[int(optimized[level * 4 + 2])])
            new_positions.append(tetrad[int(optimized[level * 4 + 3])])

            #print("{0} -> {1}".format(tetrad, new_positions))
            for i, val in enumerate(new_positions):
                tetrad[i] = new_positions[i]

            level = level + 1

    def DetermineConnectionTypes(self):
        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "":
                conn = self.nucl_quad[nucl.connected_to]
                # Tetrades are counted from bottom to top.
                # > 0 - connection from top to bottom
                # < 0 - connection from bottom to top
                # = 0 - COnnection on the same level
                level_difference = conn.tetrade_no - nucl.tetrade_no
                if conn.position == nucl.position and abs(level_difference) == 1:
                    nucl.connection_type = ConnType.SIMPLE
                elif level_difference == 0:
                    nucl.connection_type = ConnType.SAME_LEVEL
                elif (nucl.position == 2 and conn.position == 2) or \
                     (nucl.position == 3 and conn.position == 3):
                    nucl.connection_type = ConnType.RIGHT
                elif (nucl.position == 2 and conn.position == 3) or \
                     (nucl.position == 3 and conn.position == 2):
                    nucl.connection_type = ConnType.RIGHT_CROSS
                elif (nucl.position == 0 and conn.position == 0) or \
                     (nucl.position == 1 and conn.position == 1):
                    nucl.connection_type = ConnType.LEFT
                elif (nucl.position == 0 and conn.position == 1) or \
                     (nucl.position == 1 and conn.position == 0):
                    nucl.connection_type = ConnType.LEFT_CROSS
                elif (nucl.position == 0 and conn.position == 3) or \
                     (nucl.position == 3 and conn.position == 0) or \
                     (nucl.position == 1 and conn.position == 2) or \
                     (nucl.position == 2 and conn.position == 1):
                    nucl.connection_type = ConnType.FRONT_BACK_CROSS
                elif (nucl.position == 1 and conn.position == 3) or \
                     (nucl.position == 3 and conn.position == 1) or \
                     (nucl.position == 0 and conn.position == 2) or \
                     (nucl.position == 2 and conn.position == 0):
                    nucl.connection_type = ConnType.FRONT_TO_BACK


    def PrintFlow(self, chain):
        nucl = self.nucl_quad[chain["first"]]
        while nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]
            print(nucl.full_name, nucl.flow_out, nucl.flow_in)
            nucl = conn

    # TODO Optimize. Do I rly need so many iterations?
    # I HAVE TO FIX IT!! it "works"
    def CalculateFlow(self, chain):
        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "":
                conn = self.nucl_quad[nucl.connected_to]
                # Tetrades are counted from bottom to top.
                # > 0 - connection from top to bottom
                # < 0 - connection from bottom to top
                # = 0 - COnnection on the same level
                level_difference = nucl.tetrade_no - conn.tetrade_no
                if conn.position == nucl.position and abs(level_difference) == 1:
                    if level_difference == -1:
                        # Connection going bot to top
                        nucl.flow_out = ConnFlow.UP
                        nucl.flow_in = ConnFlow.UP
                        conn.flow_out = ConnFlow.UP
                    else:
                        # Connection top to bot
                        nucl.flow_out = ConnFlow.DOWN
                        nucl.flow_in = ConnFlow.DOWN
                        conn.flow_out = ConnFlow.DOWN


        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "" and nucl.flow_in == ConnFlow.UNKNOWN:
                conn = self.nucl_quad[nucl.connected_to]
                if conn.flow_out == ConnFlow.UP:
                    nucl.flow_in = ConnFlow.DOWN
                elif conn.flow_out == ConnFlow.DOWN:
                    nucl.flow_in = ConnFlow.UP

        nucl = self.nucl_quad[chain["first"]]
        while nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]

            if nucl.flow_in == ConnFlow.UNKNOWN:
                level_difference = nucl.tetrade_no - conn.tetrade_no
                if level_difference < 0:
                    # Connection going bot to top
                    nucl.flow_in = ConnFlow.DOWN
                    conn.flow_out = ConnFlow.UP
                elif level_difference > 0:
                    # Connection top to bot
                    nucl.flow_in = ConnFlow.UP
                    conn.flow_out = ConnFlow.DOWN
            nucl = conn

        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "" and nucl.flow_in == ConnFlow.UNKNOWN:
                conn = self.nucl_quad[nucl.connected_to]
                if conn.flow_out == ConnFlow.UP:
                    nucl.flow_in = ConnFlow.DOWN
                elif conn.flow_out == ConnFlow.DOWN:
                    nucl.flow_in = ConnFlow.UP

        nucl = self.nucl_quad[chain["first"]]
        while nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]
            if nucl.flow_in == ConnFlow.UNKNOWN:
                nucl.flow_in = nucl.flow_out
                if nucl.flow_out == ConnFlow.UP:
                    conn.flow_out = ConnFlow.DOWN
                elif nucl.flow_out == ConnFlow.DOWN:
                    conn.flow_out = ConnFlow.UP
            nucl = conn

        nucl = self.nucl_quad[chain["last"]]
        while nucl.connected_from != "":
            conn = self.nucl_quad[nucl.connected_from]
            if conn.flow_out == ConnFlow.UNKNOWN:
                level_difference = nucl.tetrade_no - conn.tetrade_no
                if level_difference < 0:
                    # Connection going bot to top
                    conn.flow_out = ConnFlow.DOWN
                elif level_difference > 0:
                    # Connection top to bot
                    conn.flow_out = ConnFlow.UP
                else:
                    conn.flow_out = conn.flow_in
            nucl = conn

        for _, nucl in self.nucl_quad.items():
            # Check if it is not the last nucleotide.
            if nucl.connected_to != "" and nucl.flow_in == ConnFlow.UNKNOWN:
                conn = self.nucl_quad[nucl.connected_to]
                if conn.flow_out == ConnFlow.UP:
                    nucl.flow_in = ConnFlow.DOWN
                elif conn.flow_out == ConnFlow.DOWN:
                    nucl.flow_in = ConnFlow.UP

        # Determine starting flow.
        nucl = self.nucl_quad[chain["first"]]
        if nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]
            level_difference = nucl.tetrade_no - conn.tetrade_no
            if level_difference < 0:
                # Connection going bot to top
                nucl.flow_out = ConnFlow.UP
            elif level_difference > 0:
                # Connection top to bot
                nucl.flow_out = ConnFlow.DOWN
            else:
                nucl.flow_out = nucl.flow_in

            for _, nucl in self.nucl_quad.items():
                # Check if it is not the last nucleotide.
                if nucl.connected_to != "" and nucl.flow_in == ConnFlow.UNKNOWN:
                    conn = self.nucl_quad[nucl.connected_to]
                    if conn.flow_out == ConnFlow.UP:
                        nucl.flow_in = ConnFlow.DOWN
                    elif conn.flow_out == ConnFlow.DOWN:
                        nucl.flow_in = ConnFlow.UP

        nucl = self.nucl_quad[chain["first"]]
        while nucl.connected_to != "":
            conn = self.nucl_quad[nucl.connected_to]
            if nucl.flow_in == ConnFlow.UNKNOWN:
                if nucl.flow_out == ConnFlow.UNKNOWN:
                    if nucl.connected_from == "":
                        nucl.flow_out = ConnFlow.UP
                    else:
                        conn_from = self.nucl_quad[nucl.connected_from]
                        if conn_from.flow_in == ConnFlow.DOWN:
                            nucl.flow_out = ConnFlow.UP
                        else:
                            nucl.flow_out = ConnFlow.DOWN
                nucl.flow_in = nucl.flow_out
                if nucl.flow_out == ConnFlow.UP:
                    conn.flow_out = ConnFlow.DOWN
                elif nucl.flow_out == ConnFlow.DOWN:
                    conn.flow_out = ConnFlow.UP
            nucl = conn

        #print("\n\n")
        #self.PrintFlow(chain)

class Structure:
    def __init__(self):
        self.nucleotides = {}
        # Index represents corepsponding quadruplex ID data.
        self.tetrads = []
        self.tetrads_order = []
        self.single_tetrads = []
        self.tracts = []



    """
    data - json dict of tetrade parameters.
    {
      "id":"A.DC2-A.DG10-B.DC9-B.DG3",
      "nt1":"A.DC2",
      "nt2":"A.DG10",
      "nt3":"B.DC9",
      "nt4":"B.DG3",
      "onz":"N+",
      "gbaClassification":"VIIIa",
      "planarityDeviation":0.1270031987786134,
      "ionsChannel":[
      ],
      "ionsOutside":[
      ]
    },
    """
    #def addTetrade(self, name, data, tetrade_no):
    #    self.tetrads[tetrade_no][name] = data

    """
    data - json dict of nucleotide parameters.
    {
      "index":22,
      "model":1,
      "chain":"B",
      "number":11,
      "icode":" ",
      "molecule":"DNA",
      "fullName":"B.DG11",
      "shortName":"G",
      "chi":-89.97,
      "glycosidicBond":"syn"
    }
    """
    def addNucleotide(self, name, data):
        self.nucleotides[name] = data

    def fromFile(self, path):
        with open(path) as file:
            data = json.load(file)
        return self.fromJsonDict(data)

    def fromString(self, json_string):
        data = json.loads(json_string)
        return self.fromJsonDict(data)

    def fromJsonDict(self, json_dict):
        for data in json_dict["nucleotides"]:
            self.addNucleotide(data["fullName"], data)


        # TODO What to do with other helices/quadruplexes?
        for index, helice in enumerate(json_dict["helices"]): #["quadruplexes"][0]["tetrads"].items():
            # Unordered tetrads. Order them from "tetrad_pairs"
            single_tetrads_local = []
            tetrad_unordered = {}
            tracts_all = []
            for _, quadruplex in enumerate(helice["quadruplexes"]):
                tetrad_unordered_local = {}
                for data in quadruplex["tetrads"]:
                    tetrad_unordered[data["id"]] = data
                    tetrad_unordered_local[data["id"]] = data
                single_tetrads_local.append(tetrad_unordered_local)
                if "tracts" in quadruplex:
                    tracts_all.append(quadruplex["tracts"])
                else:
                    tracts_all.append(list())

            # Order tetrads according to "tetrad_pairs" data.
            tetrad_pairs = helice["tetradPairs"]
            tetrad_ordered = []
            if tetrad_pairs != None and len(tetrad_pairs) > 0:
                # Add first pair.
                pair = tetrad_pairs.pop(0)
                tetrad_ordered.append(pair["tetrad1"])
                tetrad_ordered.append(pair["tetrad2"])
                to_temove = 0
                while len(tetrad_pairs) > 0:
                    found = False
                    for index_pair, pair in enumerate(tetrad_pairs):
                        for index, tetrad in enumerate(tetrad_ordered):
                            if tetrad == pair["tetrad1"]:
                                tetrad_ordered.insert(index + 1, pair["tetrad2"])
                                found = True
                                break
                            elif tetrad == pair["tetrad2"]:
                                tetrad_ordered.insert(index, pair["tetrad1"])
                                found = True
                                break
                        if found == True:
                            to_remove = index_pair
                            break
                    tetrad_pairs.pop(to_remove)
                    if len(tetrad_pairs) == 0:
                        break

            tetrad_ordered.reverse()

            # Do not add single tetrads as quadruplexes.
            if len(tetrad_ordered) > 1:
                self.tetrads_order.append(tetrad_ordered)
                self.tetrads.append(tetrad_unordered)

            self.tracts.append(tracts_all)
            self.single_tetrads.append(single_tetrads_local)

        # TODO check if input data was valid.
        # Are there any nucleotides and tetrads?
        return self
