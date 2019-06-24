from scipy.interpolate import RegularGridInterpolator
import gc
from scipy.interpolate import interp1d, interp2d, griddata
import numpy as np
import os
import pickle
from .filedialog import *
import h5py
import re

__author__ = "Daniel Winklehner"
__doc__ = """Class to load field data save as interpolation functions to get field data at point"""

label_dict = {"EX": "x",
              "EY": "y",
              "EZ": "z",
              "BX": "x",
              "BY": "y",
              "BZ": "z",
              "X": "x",
              "Y": "y",
              "Z": "z"}

dim_dict = {"mm": "mm",
            "MM": "mm",
            "cm": "cm",
            "CM": "cm",
            "m": "m",
            "M": "m",
            "LENGU": "cm",
            "Tesla": "Tesla",
            "TESLA": "Tesla",
            "Gauss": "Gauss",
            "GAUSS": "Gauss",
            "FLUXU": "Gauss"}

dim_numbers = {"x": 0,
               "y": 1,
               "z": 2}


class Field(object):

    def __init__(self,
                 label="Field",
                 dim=0,
                 dim_labels=None,
                 field=None,
                 scaling=1.0,
                 units="m",
                 debug=False):

        self._debug = debug

        self._label = label
        self._dim = dim  # Number of spatial dimensions. Can be 0, 1, 2, 3
        self._dim_labels = dim_labels  # list holding info on which spatial dimensions we are using
        self._filename = None
        self._unit = None
        self._scaling = scaling
        self._field = {"x": None,
                       "y": None,
                       "z": None}

        if field is not None:
            self._field = field

        self._dim_switch = {
            0: self._get_field_0d,
            1: self._get_field_1d,
            2: self._get_field_2d,
            3: self._get_field_3d
        }

        self._units_switch = {
            "m": 1.0,
            "cm": 0.01,
            "mm": 0.001
        }

        assert units in self._units_switch.keys(), "Unit not recognized!"
        self._unit_scale = self._units_switch[units]

    def __call__(self, pts):

        if not isinstance(pts, np.ndarray):
            # print("converting to numpy array")
            pts = np.array(pts)

        if not pts.ndim == 2:
            # print("converting to (1, 3) array")
            pts = np.array([pts])

        # assert len(pts) == 3, "'pts' has to be an array, list, or tuple of length three (x, y, z)!"
        assert self._field is not None, "No field data in memory!"

        return self._dim_switch[self._dim](pts)

    def __str__(self):

        if self._dim == 0:
            return "Field '{}' with {} dimensions and value {}.".format(self._label, self._dim,
                                                                        [item for _, item in self._field.items()])
        else:
            return "Field '{}' with {} dimensions.".format(self._label, self._dim)

    @staticmethod
    def cartesian_to_cylinder(_x, _y, _z):
        """
        Converts polar coordinates into cartesian coordinates.
        :param _x:
        :param _y:
        :param _z:
        :return _r, _t, _z:
        """

        _r = np.sqrt(_x * _x + _y * _y)
        _t = np.arctan2(_y, _x)

        return _r, _t, _z

    @staticmethod
    def cylinder_to_cartesian_field(_r, _t, _z, _br, _bt, _bz):
        """
        Converts polar coordinates into cartesian coordinates.
        :param _r:
        :param _t:
        :param _z:
        :param _br:
        :param _bt:
        :param _bz:
        :return _x, _y, _z:
        """

        _t_rad = np.deg2rad(_t)

        _x = _r * np.cos(_t_rad)
        _y = _r * np.sin(_t_rad)

        _bx = _br * np.cos(_t_rad) - _bt * np.sin(_t_rad)
        _by = _br * np.sin(_t_rad) + _bt * np.cos(_t_rad)

        return _x, _y, _z, _bx, _by, _bz

    @staticmethod
    def cylinder_to_cartesian(_r, _t, _z):
        """
        Converts polar coordinates into cartesian coordinates.
        :param _r:
        :param _t:
        :param _z:
        :return _x, _y, _z:
        """

        _x = _r * np.cos(np.deg2rad(_t))
        _y = _r * np.sin(np.deg2rad(_t))

        return _x, _y, _z

    def dimensions(self, dim=None):

        if dim is not None:

            assert isinstance(dim, int), "'dim' has to be integer!"
            self._dim = dim

        return self._dim

    def field(self):

        return self._field

    def _get_field_0d(self, pts):

        if len(pts) == 1:

            return self._scaling * np.array([self._field["x"],
                                             self._field["y"],
                                             self._field["z"]])

        else:

            return self._scaling * self._field["x"] * np.ones(len(pts)), \
                   self._scaling * self._field["y"] * np.ones(len(pts)), \
                   self._scaling * self._field["z"] * np.ones(len(pts))

    def _get_field_1d(self, pts):

        pts = pts[:, 2]

        if len(pts) == 1:

            return self._scaling * np.array([self._field["x"](pts)[0],
                                             self._field["y"](pts)[0],
                                             self._field["z"](pts)[0]])

        else:

            return self._scaling * self._field["x"](pts), \
                   self._scaling * self._field["y"](pts), \
                   self._scaling * self._field["z"](pts)

    def _get_field_2d(self, pts):

        pts = pts[:, :2]

        if len(pts) == 1:

            return self._scaling * np.array([self._field["x"](pts)[0],
                                             self._field["y"](pts)[0],
                                             self._field["z"](pts)[0]])

        else:

            return self._scaling * self._field["x"](pts), \
                   self._scaling * self._field["y"](pts), \
                   self._scaling * self._field["z"](pts)

    def _get_field_3d(self, pts):

        if len(pts) == 1:

            return self._scaling * np.array([self._field["x"](pts)[0],
                                             self._field["y"](pts)[0],
                                             self._field["z"](pts)[0]])

        else:

            return self._scaling * self._field["x"](pts), \
                   self._scaling * self._field["y"](pts), \
                   self._scaling * self._field["z"](pts)

    def label(self, label=None):

        if label is not None:

            assert isinstance(label, str), "'label' has to be a string!"
            self._label = label

        return self._label

    def load_field_from_file(self, filename=None, mirror=False, **kwargs):

        if filename is None:
            fd = FileDialog()
            filename = fd.get_filename("open")
            if filename is None:
                return 0

        assert os.path.exists(filename)
        print("Loading field from file '{}'".format(os.path.split(filename)[1]))
        self._filename = filename

        _, ext = os.path.splitext(filename)

        if ext == ".map":
            print("Detected AIMA Agora field file. All '*.map' files in the current folder will be used!")
            return self._load_from_agora_file(filename, mirror=mirror)

        elif ext == ".table":
            print("Detected OPERA table file from extension, loading...")
            return self._load_from_table_file(filename, **kwargs)

        elif ext == ".comsol":
            print("Detected COMSOL file from extension, loading...")
            return self._load_from_comsol_file(filename)

        elif ext == ".pickle":
            print("Detected internal 'pickle' file format, loading...")
            return self._load_from_pickle_file(filename)

        else:
            print("Did not recognize the extension. Must be one of '.map', '.table', '.comsol', '.pickle'.")
            return 1

    def _load_from_agora_file(self, filename, mirror=False):

        self._dim = 3
        interp_method = "linear"

        # Find all *.map files in the current folder
        _path = os.path.split(filename)[0]
        _all_files = []
        for entry in os.scandir(_path):
            if entry.is_file() and os.path.splitext(entry.name)[1] == ".map":
                _all_files.append(entry.name)

        _br = np.zeros([])
        _bf = np.zeros([])
        _bz = np.zeros([])

        _r = np.zeros([])
        _th = np.zeros([])
        _z = np.zeros([])

        _metadata = []
        _metadata_types = np.dtype([("fn", np.unicode_, 1024),   # filename
                                    ("fcomp", np.unicode_, 2),   # field component (bf, br, bz)
                                    ("zpos", float)])            # z position (mm)

        for _file in _all_files:
            _z_pos, res = _file.split("Z")[1].split("-")  # mm
            _b_comp = res.split(".")[0]

            if self._debug:
                print("File '{}' with component '{}' has Z position {} mm".format(_file, _b_comp, _z_pos))

            _metadata.append((os.path.join(_path, _file), _b_comp, _z_pos))

        _metadata = np.array(_metadata, dtype=_metadata_types)

        # Assert that there are an equal number of files for each component and each z position
        assert \
            np.array_equal(_metadata[np.where(_metadata["fcomp"] == "bz")]["zpos"],
                           _metadata[np.where(_metadata["fcomp"] == "br")]["zpos"]) and \
            np.array_equal(_metadata[np.where(_metadata["fcomp"] == "bz")]["zpos"],
                           _metadata[np.where(_metadata["fcomp"] == "bf")]["zpos"]), \
            "Not all components have the same z positions. Maybe a file is missing?"

        # Assert that r and theta resolution and limits are the same
        with open(_metadata["fn"][0]) as infile:
            infile.readline()
            _nth, _nr = [int(value) for value in infile.readline().strip().split()]
            _sr, _dr = [float(value) for value in infile.readline().strip().split()]
            _sth, _dth = [float(value) for value in infile.readline().strip().split()]
            old_limits = np.array([_nth, _nr, _sr, _dr, _sth, _dth])
        for _filename in _metadata["fn"][1:]:
            with open(_filename, 'r') as infile:
                infile.readline()
                _nth, _nr = [int(value) for value in infile.readline().strip().split()]
                _sr, _dr = [float(value) for value in infile.readline().strip().split()]
                _sth, _dth = [float(value) for value in infile.readline().strip().split()]
                new_limits = np.array([_nth, _nr, _sr, _dr, _sth, _dth])
                assert np.array_equal(old_limits, new_limits), "Limits for r, th in some of the files are not the same!"
                old_limits = new_limits

        _r_unique = np.round(np.linspace(_sr, _nr * _dr, _nr, endpoint=False), 10)  # cm
        _th_unique = np.round(np.linspace(_sth, _nth * _dth, _nth, endpoint=False), 10)  # degree
        _z_unique = np.sort(np.unique(_metadata["zpos"])) / 10.0  # mm --> cm
        _n_z_pos = len(_z_unique)

        _r, _th = np.meshgrid(_r_unique, _th_unique)

        _x_temp, _y_temp, _ = self.cylinder_to_cartesian(_r, _th, _z_unique)
        _x_unique = np.unique(np.round(_x_temp, 10))
        _y_unique = np.unique(np.round(_y_temp, 10))
        _x = np.round(np.linspace(min(_x_unique), max(_x_unique), 1000), 10)
        _y = np.round(np.linspace(min(_y_unique), max(_y_unique), 1000), 10)

        if mirror:
            _z_actual = np.concatenate((-_z_unique[:0:-1], _z_unique))
            idx_offset = len(_z_unique) - 1
        else:
            _z_actual = _z_unique
            idx_offset = 0

        grid_x, grid_y, grid_z = np.meshgrid(_x, _y, _z_actual, indexing='ij')
        data = {"BX": np.zeros(grid_x.shape),
                "BY": np.zeros(grid_x.shape),
                "BZ": np.zeros(grid_x.shape)}

        for j, _z in enumerate(_z_unique):

            i = j + idx_offset
            k = idx_offset - j

            _data = {}

            for _val in _metadata[np.where(_metadata["zpos"] == 10.0*_z)]:
                _filename, _fcomp, _zpos = _val

                print("Reading component {} at z = {} mm".format(_fcomp, _zpos))

                with open(_filename, 'r') as infile:
                    _raw_data = infile.readlines()[4:-1]

                _data[_fcomp] = np.array([np.fromstring(_line, sep=" ") for _line in _raw_data]).flatten()

            # TODO: All of this should be restructured to allow for cylindrical fields, optionally with symmetries
            _data["x"], _data["y"], _data["z"], _data["bx"], _data["by"], _data["bz"] = \
                self.cylinder_to_cartesian_field(_r.flatten(), _th.flatten(), _z, _data["br"], _data["bf"], _data["bz"])

            # Naturally, this is not a regular grid in x, y now...
            print()
            print("Starting griddata for slice {} ...".format(j)),

            data["BX"][:, :, i] = griddata((_data["x"], _data["y"]), _data["bx"],
                                           (grid_x[:, :, i], grid_y[:, :, i]),
                                           method=interp_method,
                                           fill_value=0.0)

            data["BY"][:, :, i] = griddata((_data["x"], _data["y"]), _data["by"],
                                           (grid_x[:, :, i], grid_y[:, :, i]),
                                           method=interp_method,
                                           fill_value=0.0)

            data["BZ"][:, :, i] = griddata((_data["x"], _data["y"]), _data["bz"],
                                           (grid_x[:, :, i], grid_y[:, :, i]),
                                           method=interp_method,
                                           fill_value=0.0)

            if mirror and j != 0:
                data["BX"][:, :, k] = -data["BX"][:, :, i]
                data["BY"][:, :, k] = -data["BY"][:, :, i]
                data["BZ"][:, :, k] = data["BZ"][:, :, i]

            # from matplotlib import pyplot as plt
            # plt.subplot(221)
            # plt.title("Bx")
            # plt.xlabel("x (cm)")
            # plt.ylabel("y (cm)")
            # cont = plt.contourf(grid_x, grid_y, bx)
            # plt.colorbar(cont)
            # plt.subplot(222)
            # plt.title("By")
            # plt.xlabel("x (cm)")
            # plt.ylabel("y (cm)")
            # cont = plt.contourf(grid_x, grid_y, by)
            # plt.colorbar(cont)
            # plt.subplot(223)
            # plt.title("Bz")
            # plt.xlabel("x (cm)")
            # plt.ylabel("y (cm)")
            # cont = plt.contourf(grid_x, grid_y, bz)
            # plt.colorbar(cont)
            # plt.show()

            print("Done!")

        for key in ["BX", "BY", "BZ"]:
            self._field[label_dict[key]] = RegularGridInterpolator(
                points=[_x * self._unit_scale,
                        _y * self._unit_scale,
                        _z_actual * self._unit_scale],
                values=data[key],
                bounds_error=False,
                fill_value=0.0)

        return 0

    def _load_from_comsol_file(self, filename):

        self._dim = 3

        label_map = {"mir1x": "X",
                     "mir1y": "Y",
                     "mir1z": "Z",
                     "mf.Bx": "BX",
                     "mf.By": "BY",
                     "mf.Bz": "BZ"}

        with open(filename, 'r') as infile:

            nlines = 9  # Number of header lines
            data = {}

            for i in range(nlines):
                line = infile.readline().strip()
                sline = line.split()
                if i == 3:
                    # self._dim = int(sline[2])
                    spatial_dims = 3
                    efield_dims = 0
                    bfield_dims = 3
                elif i == 4:
                    array_len = int(sline[2])
                elif i == 7:
                    l_unit = sline[3]
                elif i == 8:
                    if "(T)" in sline:
                        b_unit = "T"
                    j = 0
                    for label in sline:
                        if label not in ["%", "(T)"]:
                            nlabel = label_map[label]
                            data[nlabel] = {"column": j}
                            if nlabel in ["X", "Y", "Z"]:
                                data[nlabel]["unit"] = l_unit
                            elif nlabel in ["BX", "BY", "BZ"]:
                                data[nlabel]["unit"] = b_unit
                            # Another statement for EX, EY, EZ
                            if self._debug:
                                print("Header: '{}' recognized as unit {}".format(nlabel, data[nlabel]["unit"]))
                            j += 1

            _data = np.zeros([array_len, spatial_dims + efield_dims + bfield_dims])

            # Read data
            for i, line in enumerate(infile):
                _data[i] = [float(item) for item in line.split()]

        if self._debug:
            print("Data Info:")
            for key in sorted(data.keys()):
                print(key, data[key])

        _data = _data[np.lexsort(_data.T[[1]])]
        _data = _data[np.lexsort(_data.T[[0]])].T

        n = {"X": len(np.unique(_data[0])),
             "Y": len(np.unique(_data[1])),
             "Z": len(np.unique(_data[2]))}

        for key, item in data.items():
            if item["unit"] == "mm":
                item["data"] = np.unique(_data[item["column"]])

            else:
                item["data"] = np.reshape(_data[item["column"]], [n["X"], n["Y"], n["Z"]])

        del _data
        gc.collect()

        print("Creating Interpolator for {}D field...".format(self._dim))

        self._dim_labels = []

        # Process the data into the class variables
        # TODO: For now assume that field labels are either "BX", ... or "EX", ...
        label_selector = [["BX", "BY", "BZ"], ["EX", "EY", "EZ"]]
        field_type = len(np.unique([0 for key, item in data.items() if item["unit"] == "V/cm"]))

        for key in label_selector[field_type]:

            if key in data.keys():

                item = data[key]
                self._unit = item["unit"]
                self._dim_labels.append(label_dict[key])

                if self._dim == 3:
                    #  Create a Interpolator object for each of the field dimensions
                    # TODO: For now we assume that columns are labeled X, Y, Z and existed in the file

                    self._field[label_dict[key]] = RegularGridInterpolator(
                        points=[data["X"]["data"] * self._unit_scale,
                                data["Y"]["data"] * self._unit_scale,
                                data["Z"]["data"] * self._unit_scale],
                        values=item["data"],
                        bounds_error=False,
                        fill_value=0.0)

        return 0

    def _load_from_table_file(self, filename, extents=None, extents_dims=None):

        if extents is not None:
            assert extents_dims is not None, "Need to identify the extent directions"

        with open(filename, 'r') as infile:

            # Read first line with lengths (first three values after number label)
            first_line = infile.readline().strip().split()
            _n = np.array([int(val) - 1 for val in first_line])[:3]

            spatial_dims = 0
            efield_dims = 0
            bfield_dims = 0

            data = {}

            while True:

                line = infile.readline().strip()

                if line == "0":
                    break

                col_no, label, unit = line.split()
                # print(label)
                data[label] = {"column": int(col_no) - 1}

                data[label]["unit"] = dim_dict[re.split('\[|\]', unit)[1]]

                if data[label]["unit"] in ["mm", "cm", "m"]:
                    spatial_dims += 1
                    self._unit_scale = self._units_switch[data[label]["unit"]]

                elif data[label]["unit"] in ["Gauss", "Tesla"]:
                    bfield_dims += 1

                elif data[label]["unit"] in ["V/cm", "V/mm", "V/m"]:
                    efield_dims += 1

                if self._debug:
                    print("Header: '{}' recognized as unit {}".format(label, data[label]["unit"]))

            if extents is not None:
                for i, label in enumerate(extents_dims):
                    data[label] = {}
                    data[label]["column"] = int(col_no) + i
                    data[label]["unit"] = "cm"  # TODO: Pass this from the loading function - PW
                    spatial_dims += 1
                    self._unit_scale = self._units_switch[data[label]["unit"]]

            # print(spatial_dims, bfield_dims, efield_dims)

            n = {}
            # if extents_dims is None:
            if "X" in data.keys():
                n["X"] = _n[0] + 1
                if "Y" in data.keys():
                    n["Y"] = _n[1] + 1
                    n["Z"] = 1
                    if "Z" in data.keys():
                        n["Z"] = _n[2] + 1
                elif "Z" in data.keys():
                    n["Y"] = 1
                    n["Z"] = _n[1] + 1
            elif "Y" in data.keys():
                n["X"] = 1
                n["Y"] = _n[0] + 1
                n["Z"] = _n[1] + 1
            elif "Z" in data.keys():
                n["X"] = 1
                n["Y"] = 1
                n["Z"] = _n[0] + 1
            # else:
            #     for i, dim in enumerate(extents_dims):
            #         n[dim] = _n[i] + 1
            #
            #     for dim in ["X", "Y", "Z"]:
            #         if dim not in extents_dims:
            #             n[dim] = 1

            print("Sorted spatial dimension lengths: {}".format(n))

            # Generate numpy arrays holding the data:
            array_len = n["X"] * n["Y"] * n["Z"]
            print("array_len", array_len)
            # Do some assertions
            self._dim = len(np.where(_n > 1)[0])

            assert self._dim == spatial_dims, "Mismatch between detected number of spatial dimensions and " \
                                              "values in dataset.\n Potentially, only field values are saved " \
                                              "in file. \n This mode will be implemented in the future."

            assert not (efield_dims > 0 and bfield_dims > 0), "Mixed E and B fields are not supported!"

            _data = np.zeros([array_len, spatial_dims + efield_dims + bfield_dims])

            if extents is not None:
                # TODO: Assuming 3D for now... -PW
                xlims, ylims, zlims = extents[0], extents[1], extents[2]
                _x = np.linspace(xlims[0], xlims[1], n["X"])
                _y = np.linspace(ylims[0], ylims[1], n["Y"])
                _z = np.linspace(zlims[0], zlims[1], n["Z"])

                xv, yv, zv = np.meshgrid(_x, _y, _z, indexing='ij')
                _data[:, 3] = xv.ravel()
                _data[:, 4] = yv.ravel()
                _data[:, 5] = zv.ravel()

            # Read data
            if extents is None:
                for i, line in enumerate(infile):
                    _data[i] = [float(item) for item in line.split()]
            else:
                # TODO: Assuming 3D for now... -PW
                for i, line in enumerate(infile):
                    _data[i, :3] = [float(item) for item in line.split()]

        # Get limits and resolution from spatial columns if they exist
        for key, item in data.items():
            if item["unit"] in ["m", "cm", "mm"]:
                if extents is not None:
                    if key == "X":
                        rmin, rmax = extents[0]
                    elif key == "Y":
                        rmin, rmax = extents[1]
                    elif key == "Z":
                        rmin, rmax = extents[2]
                else:
                    rmin = _data[0][item["column"]]
                    rmax = _data[-1][item["column"]]

                item["limits"] = {"lower": rmin,
                                  "upper": rmax}

                item["dr"] = (rmax - rmin) / (n[key] - 1)

        # TODO: If no spatial columns exist, the user will have to provide the limits!

        if self._debug:
            print("Data Info:")
            for key in sorted(data.keys()):
                print(key, data[key])

        _data = _data.T

        for key, item in data.items():
            if item["unit"] in ["m", "cm", "mm"]:
                item["data"] = np.unique(_data[item["column"]])
            else:
                item["data"] = np.reshape(_data[item["column"]], [n["X"], n["Y"], n["Z"]])

        del _data
        gc.collect()

        print("Creating Interpolator for {}D field...".format(self._dim))

        self._dim_labels = []

        # Process the data into the class variables
        # TODO: For now assume that field labels are either "BX", ... or "EX", ...
        label_selector = [["BX", "BY", "BZ"], ["EX", "EY", "EZ"]]
        field_type = len(np.unique([0 for key, item in data.items() if item["unit"] in ["V/cm", "V/m", "V/mm"]]))

        for key in label_selector[field_type]:

            if key in data.keys():

                item = data[key]
                # print(item["data"].shape)
                self._unit = item["unit"]
                self._dim_labels.append(label_dict[key])

                if self._dim == 1:
                    #  Create a Interpolator object for each of the field dimensions
                    # TODO: For now we assume that column is labeled Z and existed in the file

                    self._field[label_dict[key]] = RegularGridInterpolator(
                        points=[data["Z"]["data"] * self._unit_scale],
                        values=item["data"].squeeze(),
                        bounds_error=False,
                        fill_value=0.0)

                elif self._dim == 2:
                    #  Create a Interpolator object for each of the field dimensions
                    # TODO: For now we assume that columns are labeled X, Y and existed in the file

                    self._field[label_dict[key]] = RegularGridInterpolator(
                        points=[data["X"]["data"] * self._unit_scale,
                                data["Y"]["data"] * self._unit_scale],
                        values=item["data"].squeeze(),
                        bounds_error=False,
                        fill_value=0.0)

                elif self._dim == 3:
                    #  Create a Interpolator object for each of the field dimensions
                    # TODO: For now we assume that columns are labeled X, Y, Z and existed in the file

                    self._field[label_dict[key]] = RegularGridInterpolator(
                        points=[data["X"]["data"] * self._unit_scale,
                                data["Y"]["data"] * self._unit_scale,
                                data["Z"][
                                    "data"] * self._unit_scale],
                        values=item["data"],
                        bounds_error=False,
                        fill_value=0.0)

            else:

                self._dim_labels.append(label_dict[key])

                if self._dim == 1:

                    _z = np.array([-1e20, 1e20])
                    _f = np.zeros(2)

                    self._field[label_dict[key]] = RegularGridInterpolator(points=[_z],
                                                                           values=_f,
                                                                           bounds_error=False,
                                                                           fill_value=0.0)

                elif self._dim == 2:

                    _x = np.array([-1e20, 1e20])
                    _y = np.array([-1e20, 1e20])
                    _f = np.zeros([2, 2])

                    self._field[label_dict[key]] = RegularGridInterpolator(points=[_x, _y],
                                                                           values=_f,
                                                                           bounds_error=False,
                                                                           fill_value=0.0)

                elif self._dim == 3:

                    _x = np.array([-1e20, 1e20])
                    _y = np.array([-1e20, 1e20])
                    _z = np.array([-1e20, 1e20])
                    _f = np.zeros([2, 2, 2])

                    self._field[label_dict[key]] = RegularGridInterpolator(points=[_x, _y, _z],
                                                                           values=_f,
                                                                           bounds_error=False,
                                                                           fill_value=0.0)

        return 0

    def _load_from_pickle_file(self, filename):

        with open(filename, "rb") as infile:
            _data = pickle.load(infile)

        for _key in self.__dict__.keys():

            if _key not in _data.keys():
                print("Not all attributes were found in the pickle file. Maybe from a different version?")
                return 1

            setattr(self, _key, _data[_key])

        return 0

    def save_to_file(self, filename=None, spacing=None, r_min=None, r_max=None):

        if filename is None:
            fd = FileDialog()
            filename = fd.get_filename("save")
            if filename is None:
                return 0

        print("Saving field to file '{}'".format(os.path.split(filename)[1]))

        _, ext = os.path.splitext(filename)

        if ext == ".h5part":
            self._save_to_h5part(filename, spacing, r_min, r_max)

        elif ext == ".pickle":
            self._save_to_pickle(filename)

        else:
            print("Did not recognize the extension, must be one of '.pickle', '.h5part'.")

        return 0

    def _save_to_pickle(self, filename):

        with open(filename, "wb") as outfile:

            _data = {}
            for _key in self.__dict__.keys():
                _data[_key] = getattr(self, _key)

            pickle.dump(_data, outfile)

        return 0

    def _save_to_h5part(self, filename, spacing, r_min, r_max):

        assert spacing is not None and r_min is not None and r_max is not None,\
            "Have to specify spacing and limits for h5part!"
        assert self._dim == 3, "Saving to h5part only implemented for 3D fields at the moment!"

        # Generate list of points to get fields at
        nr = np.array((r_max - r_min) / spacing + 1, int)
        print("r_min = ", r_min, " r_max = ", r_max)
        print("spacing = ", spacing)
        print("size = ", nr)

        x_mesh, y_mesh, z_mesh = np.meshgrid(np.linspace(r_min[0], r_max[0], nr[0]),
                                             np.linspace(r_min[1], r_max[1], nr[1]),
                                             np.linspace(r_min[2], r_max[2], nr[2]),
                                             indexing='ij', sparse=False)

        _bx = self._scaling * self._field["x"]((x_mesh, y_mesh, z_mesh))
        _by = self._scaling * self._field["y"]((x_mesh, y_mesh, z_mesh))
        _bz = self._scaling * self._field["z"]((x_mesh, y_mesh, z_mesh))

        _data = {"null": np.zeros([nr[2], nr[1], nr[0]]),
                 "hx": np.zeros([nr[2], nr[1], nr[0]]),
                 "hy": np.zeros([nr[2], nr[1], nr[0]]),
                 "hz": np.zeros([nr[2], nr[1], nr[0]])}

        for i in range(nr[2]):
            for j in range(nr[1]):
                for k in range(nr[0]):
                    _data["hx"][i, j, k] = _bx[k, j, i]
                    _data["hy"][i, j, k] = _by[k, j, i]
                    _data["hz"][i, j, k] = _bz[k, j, i]

        # Create new h5 file
        try:
            os.remove(filename)
        except Exception as _e:
            print(_e)

        h5_file = h5py.File(filename, "w")

        # Create the zeroth step and the Block inside of it
        h5_file.attrs.__setitem__("Resonance Frequency(Hz)", np.array([32800000.0]))
        step0 = h5_file.create_group("Step#0")
        block = step0.create_group("Block")

        # Create the E Field group
        e_field = block.create_group("Efield")

        # Store the x, y, and z data for the E Field
        e_field.create_dataset("0", data=_data["null"])
        e_field.create_dataset("1", data=_data["null"])
        e_field.create_dataset("2", data=_data["null"])

        # Set the spacing and origin attributes for the E Field group
        e_field.attrs.__setitem__("__Spacing__", spacing)
        e_field.attrs.__setitem__("__Origin__", r_min)

        # Create the H Field group
        h_field = block.create_group("Hfield")

        # Store the x, y, and z data points for the H Fiend
        h_field.create_dataset("0", data=_data["hx"])
        h_field.create_dataset("1", data=_data["hy"])
        h_field.create_dataset("2", data=_data["hz"])

        # Set the spacing and origin attributes for the H Field group
        h_field.attrs.__setitem__("__Spacing__", spacing)
        h_field.attrs.__setitem__("__Origin__", r_min)

        print(_data["hx"].shape)

        # Close the file
        h5_file.close()

        return 0

    def scaling(self, scaling=None):

        if scaling is not None:

            assert isinstance(scaling, float), "'scaling factor' has to be float!"
            self._scaling = scaling

        return self._scaling


# Set up some tests
if __name__ == "__main__":

    mydebug = True
    bfield1 = Field(label="Test importing AIMA field 3D",
                    debug=mydebug)
    bfield1.load_field_from_file()

    x = np.zeros(2000)
    y = np.zeros(2000)
    z = np.linspace(-15, 5, 2000)

    points = np.vstack([x, y, z]).T

    _, _, bz = bfield1(points)

    import matplotlib.pyplot as plt
    plt.plot(z, bz)
    plt.show()

    exit()
    # if platform.node() == "Mailuefterl":
    #     folder = r"D:\Daniel\Dropbox (MIT)\Projects" \
    #              r"\RFQ Direct Injection\Cyclotron"
    # elif platform.node() == "TARDIS":
    #     folder = r"D:\Dropbox (MIT)\Projects" \
    #              r"\RFQ Direct Injection\Cycloton"
    # else:
    #     folder = r"C:\Users\Daniel Winklehner\Dropbox (MIT)\Projects" \
    #              r"\RFQ Direct Injection\Cyclotron"

    folder = "D:\Daniel\Dropbox(MIT)\Projects\IsoDAR\Ion Source\MIST - 1_Magnet"

    # Test manual 0D field creation (constant everywhere)
    # bfield = Field(label="Constant Cyclotron B-Field",
    #                dim=0,
    #                field={'x': 0, 'y': 0, 'z': 1.04},
    #                debug=mydebug)
    #
    # print(bfield)
    # print(bfield(np.array([0.0, 0.0, 0.0])))

    # Test OPERA 1D field import
    bfield1 = Field(label="Test Cyclotron B-Field 1D",
                    debug=mydebug,
                    scaling=-1.0)
    bfield1.load_field_from_file(
        filename=os.path.join(folder, "BZ_vs_Z_along_z_axis.table"))

    print(bfield1)
    print(bfield1(np.array([0.0, 0.0, 0.0])))

    x = np.zeros(2000)
    y = np.zeros(2000)
    z = np.linspace(-15, 5, 2000)

    points = np.vstack([x, y, z]).T

    _, _, bz = bfield1(points)

    import matplotlib.pyplot as plt
    plt.plot(z, bz)
    plt.show()

    # Test OPERA 2D field import
    # bfield = Field(label="Test Cyclotron B-Field 2D",
    #                debug=mydebug)
    # bfield.load_field_from_file(
    #     filename=os.path.join(folder, "Bz_of_x_and_y.table"))
    #
    # print(bfield)
    # print(bfield(np.array([0.0, 0.0, 0.0])))
    #
    # x = np.linspace(-10, 10, 100)
    # y = np.linspace(-10, 10, 100)
    #
    # mesh_x, mesh_y = meshgrid(x, y, indexing='ij')
    # points = np.vstack([mesh_x.flatten(), mesh_y.flatten(), np.zeros(10000)]).T
    #
    # _, _, bz = bfield(points)
    #
    # plt.contour(x, y, bz.reshape([100, 100]), 40)
    # plt.colorbar()
    # plt.show()

    # Test OPERA 3D field import
    # bfield = Field(label="Test Cyclotron B-Field 3D",
    #                debug=mydebug,
    #                scaling=-1.0)
    # bfield.load_field_from_file(
    #     filename=os.path.join(folder, "Bx_By_Bz_of_x_y_z.table"))
    #
    # print(bfield)
    # print(bfield(np.array([0.0, 0.0, 0.0])))

    # x = np.linspace(-10, 10, 100)
    # y = np.linspace(-10, 10, 100)
    #
    # from scipy import meshgrid
    # mesh_x, mesh_y = meshgrid(x, y, indexing='ij')
    # points = np.vstack([mesh_x.flatten(), mesh_y.flatten(), np.zeros(10000)]).T
    #
    # _, _, bz = bfield(points)
    #
    # plt.contour(x, y, bz.reshape([100, 100]), 40)
    # plt.colorbar()
    # plt.show()

    # x = np.zeros(200)
    # y = np.zeros(200)
    # z = np.linspace(-15, 5, 200)
    #
    # points = np.vstack([x, y, z]).T
    #
    # _, _, bz = bfield(points)
    #
    # plt.plot(z, bz)
    # plt.show()
