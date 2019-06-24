import sys

_gui_libs = []

try:
    import gi
    gi.require_version('Gtk', '3.0')  # nopep8
    _gui_libs.append("gi")

except Exception as e:
    print("import dans_pymodules: Exception caught when trying to import gi: {}".format(e))
    print("label_combo, and mpl_canvas_wrapper will not be available!")

try:
    import PyQt5
    _gui_libs.append("qt")

except Exception as e:
    print("import dans_pymodules: Exception caught when trying to import PyQt5: {}".format(e))

if sys.version_info.major == 3:

    from .particles import *
    from .mycolors import *
    from .pylatex import *
    from .particle_distribution import *
    from .vector2d import *
    from .vector import *
    from .read_igun import *
    from .power_of_two import power_of_two
    from .coordinate_transformation_3d import *
    from .field import *
    from .particle_pusher import *

    try:
        import tkinter

        _gui_libs.append("tk")

    except Exception as e:
        print("import dans_pymodules: Exception caught when trying to import tkinter: {}".format(e))

    if "qt" in _gui_libs:
        from .filedialog_qt import *
    elif "gi" in _gui_libs:
        from .filedialog_gtk import *
    else:
        from .filedialog_tk import *

    if "gi" in _gui_libs:
        from .label_combo import *
        from .mpl_canvas_wrapper import *

elif sys.version_info.major == 2:

    from particles import *
    from mycolors import *
    from pylatex import *
    from particle_distribution import *
    from vector2d import *
    from vector import *
    from read_igun import *
    from power_of_two import *
    from coordinate_transformation_3d import *
    from field import *
    from particle_pusher import *

    try:
        import Tkinter

        _gui_libs.append("tk")

    except Exception as e:
        print("import dans_pymodules: Exception caught when trying to import Tkinter: {}".format(e))

    if "qt" in _gui_libs:
        from filedialog_qt import *
    elif "gi" in _gui_libs:
        from filedialog_gtk import *
    else:
        from filedialog_tk import *

    if "gi" in _gui_libs:
        from label_combo import *
        from mpl_canvas_wrapper import *

else:
    raise Exception("This version of python is not supported!")
