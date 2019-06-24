import h5py
import numpy as np


class TableToH5(object):
    def __init__(self, spacing, r_min, r_max, filename='my_h5'):

        self.spacing = spacing
        self.r_min = r_min
        self.r_max = r_max
        self.filename = filename

        # Create a new h5 file
        self.h5_file = h5py.File(filename + '.h5part', )

        # Calculates the size of data arrays
        # noinspection PyTypeChecker
        self._size = np.array((r_max - r_min) / spacing + 1, int)

        # Initialize the h5 data
        # Data Format:
        # Dictionary with keys "ex", "ey", "ez", "hx", "hy", and "hz", which correspond to the vector components
        # of the electric field and the H field.
        self.data = {"ex": np.zeros(self._size),
                     "ey": np.zeros(self._size),
                     "ez": np.zeros(self._size),
                     "hx": np.zeros(self._size),
                     "hy": np.zeros(self._size),
                     "hz": np.zeros(self._size)}

    def set_data(self):

        with open(self.filename + '.table', 'r') as table:
            h = 0
            _s = self._size[0] * self._size[1] * self._size[2]
            _lines = table.readlines()[5:]
            _ex, _ey, _ez = np.zeros(_s), np.zeros(_s), np.zeros(_s)
            # _hx, _hy, _hz = np.zeros(_s), np.zeros(_s), np.zeros(_s)
            for line in _lines:
                _tmp = line.lstrip().rstrip().split()
                _ex[h], _ey[h], _ez[h] = float(_tmp[0]), float(_tmp[1]), float(_tmp[2])
                # _hx, _hy, _hz = float(_tmp[0]), float(_tmp[1]), float(_tmp[2])
                h += 1

        for i in range(self._size[2]):
            for j in range(self._size[1]):
                for k in range(self._size[0]):
                    self.data["ex"][k, j, i] = _ex[k + j * self._size[2] + i * self._size[2] * self._size[1]] * 1e-4
                    self.data["ey"][k, j, i] = _ey[k + j * self._size[2] + i * self._size[2] * self._size[1]] * 1e-4
                    self.data["ez"][k, j, i] = _ez[k + j * self._size[2] + i * self._size[2] * self._size[1]] * 1e-4
                    # self.data["hx"][i, j, k] = _hx[i + j * self._size[2] + k * self._size[2] * self._size[1]] * 1e-3
                    # self.data["hy"][i, j, k] = _hy[i + j * self._size[2] + k * self._size[2] * self._size[1]] * 1e-3
                    # self.data["hz"][i, j, k] = _hz[i + j * self._size[2] + k * self._size[2] * self._size[1]] * 1e-3

    def generate(self):

        # Create the zeroth step and the Block inside of it
        self.h5_file.attrs.__setitem__("Resonance Frequency(Hz)", np.array([49200000.0]))
        step0 = self.h5_file.create_group("Step#0")
        block = step0.create_group("Block")

        # Create the E Field group
        e_field = block.create_group("Efield")

        # Store the x, y, and z data for the E Field
        e_field.create_dataset("0", data=self.data["ex"])
        e_field.create_dataset("1", data=self.data["ey"])
        e_field.create_dataset("2", data=self.data["ez"])

        # Set the spacing and origin attributes for the E Field group
        e_field.attrs.__setitem__("__Spacing__", self.spacing)
        e_field.attrs.__setitem__("__Origin__", self.r_min)

        # Create the H Field group
        h_field = block.create_group("Hfield")

        # Store the x, y, and z data points for the H Fiend
        h_field.create_dataset("0", data=self.data["hx"])
        h_field.create_dataset("1", data=self.data["hy"])
        h_field.create_dataset("2", data=self.data["hz"])

        # Set the spacing and origin attributes for the H Field group
        h_field.attrs.__setitem__("__Spacing__", self.spacing)
        h_field.attrs.__setitem__("__Origin__", self.r_min)

        # Close the file
        self.h5_file.close()

    def _set_uniform_bfield(self, tesla=None, kgauss=None):
        # Set the magnetic field in the "hz" direction
        if tesla is not None:
            self.data["hz"][:, :, :] = tesla * 10.0

        elif kgauss is not None:
            self.data["hz"][:, :, :] = kgauss


class COMSOLToH5(object):
    def __init__(self, spacing, r_min, r_max, filename='my_h5'):

        self.spacing = spacing
        self.r_min = r_min
        self.r_max = r_max
        self.filename = filename

        # Create a new h5 file
        self.h5_file = h5py.File(filename + '.h5part', )

        # Calculates the size of data arrays
        # noinspection PyTypeChecker
        self._size = np.array((r_max - r_min) / spacing + 1, int)
        print("Size = ", self._size)
        # Initialize the h5 data
        # Data Format:
        # Dictionary with keys "ex", "ey", "ez", "hx", "hy", and "hz", which correspond to the vector components
        # of the electric field and the H field.
        self.data = {"ex": np.zeros([self._size[2], self._size[1], self._size[0]]),
                     "ey": np.zeros([self._size[2], self._size[1], self._size[0]]),
                     "ez": np.zeros([self._size[2], self._size[1], self._size[0]]),
                     "hx": np.zeros([self._size[2], self._size[1], self._size[0]]),
                     "hy": np.zeros([self._size[2], self._size[1], self._size[0]]),
                     "hz": np.zeros([self._size[2], self._size[1], self._size[0]])}

    def set_data(self):

        with open(self.filename + '.txt', 'r') as table:

            h = 0
            _s = self._size[0] * self._size[1] * self._size[2]
            _lines = table.readlines()[9:]
            # _x, _y, _z = np.zeros(_s), np.zeros(_s), np.zeros(_s)
            _ex, _ey, _ez = np.zeros(_s), np.zeros(_s), np.zeros(_s)
            # _hx, _hy, _hz = np.zeros(_s), np.zeros(_s), np.zeros(_s)

            for line in _lines:
                # [X] [Y] [Z] [EX] [EY] [EZ]
                _tmp = line.lstrip().rstrip().split()
                # _xy_values = [float(_tmp[0]), float(_tmp[1]), float(_tmp[2])]
                _values = [float(_tmp[3]), float(_tmp[4]), float(_tmp[5])]
                for i in range(3):
                    if np.isnan(_values[i]):
                        _values[i] = 0.0
                    # if np.isnan(_xy_values[i]):
                    #     _xy_values[i] = 0.0
                # _x[h], _y[h], _z[h] = _xy_values
                _ex[h], _ey[h], _ez[h] = _values
                # _hx, _hy, _hz = float(_tmp[0]), float(_tmp[1]), float(_tmp[2])
                h += 1

        for i in range(self._size[2]):
            for j in range(self._size[1]):
                for k in range(self._size[0]):
                    self.data["ex"][i, j, k] = _ex[k + j * self._size[0] + i * self._size[0] * self._size[1]] * 1e-6
                    self.data["ey"][i, j, k] = _ey[k + j * self._size[0] + i * self._size[0] * self._size[1]] * 1e-6
                    self.data["ez"][i, j, k] = _ez[k + j * self._size[0] + i * self._size[0] * self._size[1]] * 1e-6
                    # self.data["hx"][i, j, k] = _hx[i + j * self._size[2] + k * self._size[2] * self._size[1]] * 1e-3
                    # self.data["hy"][i, j, k] = _hy[i + j * self._size[2] + k * self._size[2] * self._size[1]] * 1e-3
                    # self.data["hz"][i, j, k] = _hz[i + j * self._size[2] + k * self._size[2] * self._size[1]] * 1e-3

        print(self.data["ex"].shape)

    def generate(self):
        # Create the zeroth step and the Block inside of it
        self.h5_file.attrs.__setitem__("Resonance Frequency(Hz)", np.array([49200000.0]))
        step0 = self.h5_file.create_group("Step#0")
        block = step0.create_group("Block")

        # Create the E Field group
        e_field = block.create_group("Efield")

        # Store the x, y, and z data for the E Field
        e_field.create_dataset("0", data=self.data["ex"])
        e_field.create_dataset("1", data=self.data["ey"])
        e_field.create_dataset("2", data=self.data["ez"])

        # Set the spacing and origin attributes for the E Field group
        e_field.attrs.__setitem__("__Spacing__", self.spacing)
        e_field.attrs.__setitem__("__Origin__", self.r_min)

        # Create the H Field group
        h_field = block.create_group("Hfield")

        # Store the x, y, and z data points for the H Fiend
        h_field.create_dataset("0", data=self.data["hx"])
        h_field.create_dataset("1", data=self.data["hy"])
        h_field.create_dataset("2", data=self.data["hz"])

        # Set the spacing and origin attributes for the H Field group
        h_field.attrs.__setitem__("__Spacing__", self.spacing)
        h_field.attrs.__setitem__("__Origin__", self.r_min)

        # Close the file
        self.h5_file.close()

    def _set_uniform_bfield(self, tesla=None, kgauss=None):

        # Set the magnetic field in the "hz" direction
        if tesla is not None:
            self.data["hz"][:, :, :] = tesla * 10.0

        elif kgauss is not None:
            self.data["hz"][:, :, :] = kgauss


class createBFieldMap(object):

    def __init__(self, spacing, r_min, r_max, filename='dummy_field'):
        self.spacing = spacing
        self.r_min = r_min
        self.r_max = r_max
        self.filename = filename

        # Create a new h5 file
        self.h5_file = h5py.File(filename + '.h5part', )

        # Calculates the size of data arrays
        # noinspection PyTypeChecker
        self._size = np.array((r_max - r_min) / spacing + 1, int)

        # Initialize the h5 data
        # Data Format:
        # Dictionary with keys "ex", "ey", "ez", "hx", "hy", and "hz", which correspond to the vector components
        # of the electric field and the H field.
        self.data = {"ex": np.zeros(self._size),
                     "ey": np.zeros(self._size),
                     "ez": np.zeros(self._size),
                     "hx": np.zeros(self._size),
                     "hy": np.zeros(self._size),
                     "hz": np.zeros(self._size)}

    def generate(self):

        # Create the zeroth step and the Block inside of it
        self.h5_file.attrs.__setitem__("Resonance Frequency(Hz)", np.array([49200000.0]))
        step0 = self.h5_file.create_group("Step#0")
        block = step0.create_group("Block")

        # Create the E Field group
        e_field = block.create_group("Efield")

        # Store the x, y, and z data for the E Field
        e_field.create_dataset("0", data=self.data["ex"])
        e_field.create_dataset("1", data=self.data["ey"])
        e_field.create_dataset("2", data=self.data["ez"])

        # Set the spacing and origin attributes for the E Field group
        e_field.attrs.__setitem__("__Spacing__", self.spacing)
        e_field.attrs.__setitem__("__Origin__", self.r_min)

        # Create the H Field group
        h_field = block.create_group("Hfield")

        # Store the x, y, and z data points for the H Fiend
        h_field.create_dataset("0", data=self.data["hx"])
        h_field.create_dataset("1", data=self.data["hy"])
        h_field.create_dataset("2", data=self.data["hz"])

        # Set the spacing and origin attributes for the H Field group
        h_field.attrs.__setitem__("__Spacing__", self.spacing)
        h_field.attrs.__setitem__("__Origin__", self.r_min)

        # Close the file
        self.h5_file.close()

    def _set_uniform_bfield(self, tesla=None, kgauss=None):
        # Set the magnetic field in the "hz" direction
        if tesla is not None:
            self.data["hz"][:, :, :] = tesla * 10.0

        elif kgauss is not None:
            self.data["hz"][:, :, :] = kgauss


if __name__ == '__main__':

    # Spacing and origin attributes
    # spacing = np.array([20.0, 20.0, 20.0])
    # r_min = np.array([-100.0, -100.0, -100.0])
    # r_max = np.array([100.0, 100.0, 100.0])

    spacing = np.array([1.0, 1.0, 1.0])
    r_min = np.array([-250.0, -250.0, -30.0])
    r_max = np.array([250.0, 250.0, 30.0])
    filename = r"C:\Users\Daniel Winklehner\Dropbox (MIT)\Projects\RFQ Direct" \
               r" Injection\Comsol\AIMA CR Design\AIMA_80kV_RF_no_SI_2mm"

    # Assumes that the .table filename is the same as the filename you want to save the h5 to.
    # filename = '/home/philip/src/dans_pymodules/dans_pymodules/test_fieldmaps/plate_capacitor_11x11x11_test'

    # my_h5 = TableToH5(spacing=spacing, r_min=r_min, r_max=r_max, filename=filename)
    # my_h5.set_data()
    # # my_h5 = createBFieldMap(spacing, r_min, r_max, filename=filename)
    # my_h5._set_uniform_bfield(tesla=1.041684)
    my_h5 = COMSOLToH5(spacing=spacing, r_min=r_min, r_max=r_max, filename=filename)
    my_h5.set_data()
    my_h5.generate()
