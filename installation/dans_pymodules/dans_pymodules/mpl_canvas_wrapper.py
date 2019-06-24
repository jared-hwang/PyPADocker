import gi
gi.require_version('Gtk', '3.0')  # nopep8
from gi.repository import Gtk, GObject
import matplotlib
from matplotlib.patches import Rectangle
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.artist import ArtistInspector
from matplotlib.lines import Line2D
from matplotlib.collections import PathCollection
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg
import os
import numpy as np

__author__ = "Daniel Winklehner"
__doc__ = "Convenience class to generate a GUI window for matplotlib plots"


def generate_default_settings():
    """
    Generates a set of default values to use as entries in the Settings
    Dialog
    """
    entries = {"title": "Title",
               "x_label": "x",
               "y1_label": "y",
               "y2_label": "secondary y",
               "x_min": float(-1),
               "y1_min": float(-1),
               "y2_min": float(-1),
               "x_max": float(1),
               "y1_max": float(1),
               "y2_max": float(1)}

    adjustments = {"major_ticks": 15,
                   "label_fs": 30,
                   "legend_ms": 1,
                   "legend_fs": 18,
                   "legend_lw": 1,
                   "title_fs": 30,
                   "line_width": 1,
                   "marker_size": 25
                   }

    legend = []

    defaults = {"entry": entries,
                "adj": adjustments,
                "legend_entries": legend,
                "legend_global_flag": False,
                "legend_entries_flags": np.array([], 'bool'),
                "autoscalex": True,
                "autoscaley1": True,
                "autoscaley2": True,
                "xtime": False}

    return defaults


def timestr(s):
    """
    Takes a value in seconds (s) and converts it into a string of the
    format hh:mm:ss in 24 hour time.
    If the number of seconds exceeds a day, it starts over from the
    beginning. This assumes that days is the largest unit of time in the
    context.
    """

    # Number of days (if any)
    days = int(s / 86400)

    if days > 0:
        s -= days * 86400

    h = int(s / 3600)
    m = int((s - 3600 * h) / 60)
    s = int(s - 3600 * h - 60 * m)

    if h < 10:
        hstr = "0%i" % h
    else:
        hstr = "%i" % h

    if m < 10:
        mstr = "0%i" % m
    else:
        mstr = "%i" % m

    if s < 10:
        sstr = "0%i" % s
    else:
        sstr = "%i" % s

    ts = hstr + ":" + mstr + ":" + sstr

    return ts


class SettingsDialog(object):
    """
    Dialog box for all the different plot settings in the WARP postprocessor
    """

    def __init__(self):
        """
        Constructor of the SettingsDialog class
        """
        self.gladefile = os.path.join(os.path.dirname(__file__), "PlotSettingsDialog.glade")

        # Uncomment this line for cxFreeze binary building
        # self.gladefile  = "PlotSettingsDialog.glade"

        # --- Generate Dialog GUI from gladefile --- #
        self.builder = Gtk.Builder()
        self.builder.add_from_file(self.gladefile)
        self.dialog = self.builder.get_object("settings_dialog")

        self.settings = None

    def deactivate_autoscale_x(self, *_):
        """
        """

        self.builder.get_object("autoscalex_cb").set_active(False)

        return 0

    def deactivate_autoscale_y(self, *_):
        """
        """

        self.builder.get_object("autoscaley1_cb").set_active(False)

        return 0

    def deactivate_secondary_autoscale_y(self, *_):
        """
        """

        self.builder.get_object("autoscaley2_cb").set_active(False)

        return 0

    def run(self, parent=None, old_settings=None):
        """
        Function that shows the dialog and gives back the dict with all the
        new settings
        """
        if old_settings is None:

            old_settings = generate_default_settings()

        self.settings = old_settings.copy()

        if parent is not None:

            self.dialog.set_transient_for(parent)

        self.builder.connect_signals({"deactivate_autoscale_x": self.deactivate_autoscale_x,
                                      "deactivate_autoscale_y1": self.deactivate_autoscale_y,
                                      "deactivate_autoscale_y2": self.deactivate_secondary_autoscale_y})

        # --- Fill in the entries and spinboxes with the current values --- #
        for key in old_settings["entry"].keys():
            
            self.builder.get_object(key + "_entry").set_text(str(old_settings["entry"][key]))

        for key in old_settings["adj"].keys():
            
            self.builder.get_object(key + "_adj").set_value(old_settings["adj"][key])

        # --- Make a VBox for the legend entries --- #
        legend_table = self.builder.get_object("legend_table")

        legend_entries = []
        legend_entry_show = []

        old_legend = self.settings["legend_entries"]
        old_legend_flags = self.settings["legend_entries_flags"]

        n_legend_entries = len(old_legend)

        legend_table.resize(n_legend_entries + 1, 4)

        for i in range(len(old_legend)):
            l_entry = Gtk.Entry()
            l_flag_cb = Gtk.CheckButton()
            l_flag_cb.set_active(old_legend_flags[i])
            l_flag_cb.show()
            # l_flag_cb.set_sensitive(False)
            l_entry.set_size_request(300, -1)
            l_entry.set_text(old_legend[i])
            # l_entry.set_activates_default(True)

            legend_entries.append(l_entry)
            legend_entry_show.append(l_flag_cb)

            legend_table.attach(l_entry, 0, 1, i + 1, i + 2)
            legend_table.attach(l_flag_cb, 1, 2, i + 1, i + 2)

        legend_table.show_all()

        if old_settings["legend_global_flag"]:
            self.builder.get_object("legend_cb").set_active(True)
        if old_settings["autoscalex"]:
            self.builder.get_object("autoscalex_cb").set_active(True)
        if old_settings["autoscaley1"]:
            self.builder.get_object("autoscaley1_cb").set_active(True)
        if old_settings["autoscaley2"]:
            self.builder.get_object("autoscaley2_cb").set_active(True)

        if old_settings["xtime"]:

            self.builder.get_object("x_min_entry").set_text(
                timestr(float(self.builder.get_object("x_min_entry").get_text())))

            self.builder.get_object("x_max_entry").set_text(
                timestr(float(self.builder.get_object("x_max_entry").get_text())))

        response = self.dialog.run()

        if response == 1:

            for key in old_settings["entry"].keys():

                if key in ["title", "x_label", "y1_label", "y2_label"]:

                    self.settings["entry"][key] = self.builder.get_object(key + "_entry").get_text()

                else:

                    text = self.builder.get_object(key + "_entry").get_text()

                    if old_settings["xtime"]:

                        if key in ["x_min", "x_max"]:

                            try:

                                h, m, s = text.split(":")
                                self.settings["entry"][key] = int(h) * 3600 + int(m) * 60 + float(s)

                            except Exception as e:

                                print("Exception caught: {}".format(e))

                                return -1, None

                        else:

                            if text == "":
                                text = "0.0"

                            self.settings["entry"][key] = float(text)

                    else:

                        if text == "":
                            text = "0.0"

                        self.settings["entry"][key] = float(text)

            for key in old_settings["adj"].keys():

                self.settings["adj"][key] = self.builder.get_object(key + "_adj").get_value()

            if self.builder.get_object("legend_cb").get_active():
                self.settings["legend_global_flag"] = True
            else:
                self.settings["legend_global_flag"] = False
            if self.builder.get_object("autoscalex_cb").get_active():
                self.settings["autoscalex"] = True
            else:
                self.settings["autoscalex"] = False
            if self.builder.get_object("autoscaley1_cb").get_active():
                self.settings["autoscaley1"] = True
            else:
                self.settings["autoscaley1"] = False
            if self.builder.get_object("autoscaley2_cb").get_active():
                self.settings["autoscaley2"] = True
            else:
                self.settings["autoscaley2"] = False

            for i in range(len(legend_entries)):
                self.settings["legend_entries"][i] = legend_entries[i].get_text()
                self.settings["legend_entries_flags"][i] = legend_entry_show[i].get_active()

            self.dialog.destroy()

            return response, self.settings

        else:

            self.dialog.destroy()

            return response, None


class MPLCanvasWrapper(Gtk.VBox):
    """
    Class that wraps around the Matplotlib Canvas/Figure/Subplot to be used with
    GTK3.
    Supports the following:
            -) Zoom (rectangle, scroll wheel)
            -) Display coordinates
            -) pipes changing the plot settings through to the mpl figure
    """
    __gsignals__ = {"update_request": (GObject.SIGNAL_RUN_FIRST, GObject.TYPE_NONE, (int,)), }

    def append_legend_entries_flag(self, flag=True):

        self.show_legend_entries = np.append(self.show_legend_entries, [flag])

    def reset_legend_entries_flags(self):

        self.show_legend_entries = np.array([], 'bool')

    def change_settings(self, *_):
        """
        Callback for the settings button.
        Calls a popup window from Dialogs.py and asks the user to specify
        plot settings
        """

        self.zoom_out(self)

        sd = SettingsDialog()

        response, settings = sd.run(parent=self._main_window, old_settings=self.get_settings())

        if response == 1:

            self.set_settings(settings)
            self.emit("update_request", self.nbp)

            return settings

        else:

            return None

    def draw_event(self, *_):
        """
        Event handler for any draw operation on the canvas
        Used here to set dummy x-axis to the same scaling as the primary x-axis
        (workaround for coordinate display with mouse-over events)
        """

        self.set_dummy_xlim()

        return 0

    def unhide_secondary_axis(self):
        """
        Secondary axis cannot be deleted without remakeing the whole canvas,
        therefore, it will only be hidden if it is not used...
        """
        if self.show_secondary_axis:
            return 0

        self.show_secondary_axis = True
        self.secondary_axis.yaxis.set_visible(True)
        # self.secondary_axis.set_frame_on(True)

        self.draw_idle()

        return 0

    def hide_secondary_axis(self):
        """
        Secondary axis cannot be deleted without remakeing the whole canvas,
        therefore, it will only be hidden if it is not used...
        """
        # self.secondary_axis.set_ylabel("")
        # self.secondary_axis.set_yticks([])
        if not self.show_secondary_axis:
            return 0

        self.show_secondary_axis = False
        self.secondary_axis.yaxis.set_visible(False)
        # self.secondary_axis.set_frame_on(False)

        self.draw_idle()

        return 0

    def zoom_out(self, *_):
        """
        """

        if self.zoom_flag:

            self.zoom_flag = False

            self.axis.set_xlim(self.maxlimits[:2])
            self.axis.set_ylim(self.maxlimits[2:])

            self.set_autoscale(self.autoscale_old[0], "x")
            self.set_autoscale(self.autoscale_old[1], "y")

            if self.xtime:
                self.set_xtime()

            self.canvas.draw_idle()

        return 0

    def draw_zoom_box(self, *_):
        """
        Draw the zoom to rect rubberband rectangle
        """

        if not (self.box_start is None or self.box_end is None):

            if self.zoom_box is not None:
                self.zoom_box.remove()

            xy = (min([self.box_start[0], self.box_end[0]]),
                  min([self.box_start[1], self.box_end[1]]))

            width = abs(self.box_start[0] - self.box_end[0])
            height = abs(self.box_start[1] - self.box_end[1])

            self.zoom_box = Rectangle(xy=xy,
                                      width=width, height=height,
                                      linestyle='dashed', linewidth=1,
                                      color='grey', fill=False)

            self.axis.add_artist(self.zoom_box)

            self.canvas.draw_idle()

        return 0

    def mouse_move(self, event):
        """
        """

        if event.inaxes:

            x, y = event.xdata, event.ydata

            self.z_label.set_text("| %.4e" % x)
            self.r_label.set_text("%.4e | " % y)

            if self.holddown:

                if self.box_start_px[0] != event.x and \
                                self.box_start_px[1] != event.y:
                    self.box_end = x, y
                    self.draw_zoom_box(self)

            self.inaxes = True
            self.mouse_x = x
            self.mouse_y = y

        else:
            self.inaxes = False

        return 0

    # TODO: Unused?
    @staticmethod
    def axes_enter_callback(event):

        x, y = event.xdata, event.ydata
        del x, y

        return 0

    @staticmethod
    def axes_leave_callback(event):

        x, y = event.xdata, event.ydata
        del x, y

        return 0

    def button_pressed(self, event):
        """
        """

        leg = self.get_legend()

        if event.button == 3:
            # Action on pressing right mouse button inside canvas
            if leg is not None:

                start, end = leg.get_window_extent().get_points()
                x_px, y_px = event.x, event.y

                if start[0] <= x_px <= end[0] and start[1] <= y_px <= end[1]:
                    return 1

                else:
                    self.zoom_out(self)

            else:
                self.zoom_out(self)

        elif event.button == 1:
            # Action on pressing left mouse button inside canvas
            if leg is not None:

                start, end = leg.get_window_extent().get_points()
                x_px, y_px = event.x, event.y

                if start[0] <= x_px <= end[0] and start[1] <= y_px <= end[1]:
                    return 1

            self.box_start = None

            if event.inaxes:
                self.holddown = True
                self.box_start = event.xdata, event.ydata
                self.box_start_px = event.x, event.y

        elif event.button == 2:
            # Action on pressing mouse wheel button inside canvas
            pass

        return 0

    def zoom(self, *_):
        """
        """

        if not self.zoom_flag:
            self.maxlimits = self.get_limits()[:-2]
            self.zoom_flag = True
            self.autoscale_old = [self.get_autoscalex(), self.get_autoscaley1()]

        self.zoomlimits = [min([self.box_start[0], self.box_end[0]]),
                           max([self.box_start[0], self.box_end[0]]),
                           min([self.box_start[1], self.box_end[1]]),
                           max([self.box_start[1], self.box_end[1]])]

        self.axis.set_xlim(self.zoomlimits[:2])
        self.axis.set_ylim(self.zoomlimits[2:])

        if self.xtime:
            self.set_xtime()

        self.canvas.draw_idle()

        return 0

    def button_released(self, event):
        """
        """

        if event.button == 1:

            if self.holddown:

                if self.zoom_box is not None:
                    self.zoom_box.remove()

                if self.box_start_px[0] != event.x and self.box_start_px[1] != event.y:
                    self.zoom(self)

                self.holddown = False
                self.box_start = None
                self.box_start_px = None
                self.box_end = None
                self.zoom_box = None

        return 0

    def mouse_scroll(self, event):
        """
        Cave: Something wrong with the aspect ratio!
        """
        zoom_perc = 0.05

        if not event.inaxes:
            return 0

        x = event.xdata
        y = event.ydata

        lim = self.get_limits()

        # first zoom - set maxlimits
        if not self.zoom_flag:
            self.maxlimits = self.get_limits()[:-2]
            self.zoom_flag = True
            self.autoscale_old = [self.get_autoscalex(), self.get_autoscaley1()]

        # Zoom out
        if event.step < 0:

            zoom_perc = 1 + zoom_perc

        # Zoom in
        elif event.step > 0:

            zoom_perc = 1 - zoom_perc

        xmax = x + zoom_perc * (lim[1] - x)
        xmin = x + zoom_perc * (lim[0] - x)
        ymax = y + zoom_perc * (lim[3] - y)
        ymin = y + zoom_perc * (lim[2] - y)

        self.axis.set_xlim([xmin, xmax])
        self.axis.set_ylim([ymin, ymax])

        if self.xtime:
            self.set_xtime()

        self.canvas.draw_idle()

        return 0

    def plot(self, *args, **kwargs):
        """
        """
        if "show_in_legend" in kwargs:

            self.append_legend_entries_flag(kwargs["show_in_legend"])
            kwargs.pop("show_in_legend")

        else:

            self.append_legend_entries_flag(True)

        # If there are no previous plots, get the marker size from saved main_plot_settings
        plots, dummy = self.axis.get_legend_handles_labels()

        if len(plots) == 0:
            kwargs["linewidth"] = self.main_plot_settings["adj"]["line_width"]
            kwargs["markersize"] = self.main_plot_settings["adj"]["marker_size"]

        # Handle secondary axis
        if "secondary_axis" in kwargs:

            if kwargs["secondary_axis"]:
                kwargs.pop("secondary_axis")

                self.unhide_secondary_axis()

                return self.secondary_axis.plot(*args, **kwargs)

            kwargs.pop("secondary_axis")

        return self.axis.plot(*args, **kwargs)

    def scatter(self, *args, **kwargs):
        """
        """
        if "show_in_legend" in kwargs:

            self.append_legend_entries_flag(kwargs["show_in_legend"])
            kwargs.pop("show_in_legend")

        else:

            self.append_legend_entries_flag(True)

        # If there are no previous plots, get the marker size from saved main_plot_settings
        plots, dummy = self.axis.get_legend_handles_labels()

        if len(plots) == 0:
            kwargs["s"] = self.main_plot_settings["adj"]["marker_size"]

        if "secondary_axis" in kwargs:

            if kwargs["secondary_axis"]:
                kwargs.pop("secondary_axis")

                self.unhide_secondary_axis()
                return self.secondary_axis.scatter(*args, **kwargs)

            kwargs.pop("secondary_axis")

        return self.axis.scatter(*args, **kwargs)

    def toggle_legend(self):
        """
        """
        leg = self.get_legend()

        if leg is not None:

            if leg.get_visible():

                leg.set_visible(False)

            else:

                leg.set_visible(True)

        else:

            self.set_legend(flag=True)

        self.draw_idle()

        return 0

    def reset_secondary_axes(self):
        """
        """
        print(self.figure.axes)

        self.figure.delaxes(self.secondary_axis)
        self.figure.delaxes(self.dummy_axis)

        # --- The secondary y-axis (hidden by default)
        self.secondary_axis = self.axis.twinx()
        self.secondary_axis.set_frame_on(False)
        self.show_secondary_axis = True
        self.hide_secondary_axis()

        # --- A dummy secondary x axis for mousover events
        self.dummy_axis = self.axis.twiny()
        self.dummy_axis.set_frame_on(False)
        self.dummy_axis.set_xticks([])
        self.set_dummy_xlim()

        print(self.figure.axes)

        return 0

    def set_settings(self, settings=None):
        """
        Sets the settings according to the values in Settings-array.
        """
        if settings is not None:

            # Get the containers for the different values
            ent = settings["entry"]
            adj = settings["adj"]
            self.main_plot_settings = settings

            # Set the plots' settings (line width, colors etc.)
            self.set_linewidth(adj["line_width"])
            self.set_markersize(adj["marker_size"])

            # Set the Canvas (axes, fontsizes, legend, etc.) settings
            self.set_title(ent["title"], adj["title_fs"])
            self.set_xlabel(ent["x_label"], adj["label_fs"])
            self.set_ylabel(ent["y1_label"], adj["label_fs"])
            self.set_secondary_ylabel(ent["y2_label"], adj["label_fs"])

            self.set_ticksize(adj["major_ticks"])

            # Set the plot labels (this is where the legend entries are stored)
            handles1, labels1 = self.axis.get_legend_handles_labels()
            handles2, labels2 = self.secondary_axis.get_legend_handles_labels()

            handles = handles1 + handles2
            # labels = labels1 + labels2

            if "legend_entries_flags" in settings:
                self.show_legend_entries = settings["legend_entries_flags"]

            if "legend_entries" in settings:

                for handle, label in zip(handles, settings["legend_entries"]):
                    handle.set_label(label)

                self.set_legend(settings["legend_global_flag"],
                                settings["legend_entries"],
                                adj["legend_fs"],
                                adj["legend_ms"],
                                adj["legend_lw"])

            if settings["autoscalex"]:

                self.set_autoscale(True, "x")

            else:

                self.set_xlim(ent["x_min"], ent["x_max"])

            if settings["autoscaley1"]:

                self.set_autoscale(True, "y")

            else:

                self.set_ylim(ent["y1_min"], ent["y1_max"])

            if settings["autoscaley2"]:

                self.set_secondary_autoscale(True, "y")

            else:

                self.set_secondary_ylim(ent["y2_min"], ent["y2_max"])

            if self.xtime:
                self.set_xtime()

            if not self.show_secondary_axis:
                self.hide_secondary_axis()

            self.canvas.draw_idle()

        return 0

    def set_legend(self, flag=False, entries=None, fs=18, ms=1, lw=1):
        """
        Handles the creation and destruction of the legend object
        flag......whether or not the legend is to be displayed
        entries...the actual legend entries
        ms........the legend marker size
        fs........the legend entries' font size
        lw........the legend linewidth
        """
        del entries

        handles1, labels1 = self.axis.get_legend_handles_labels()
        handles2, labels2 = self.secondary_axis.get_legend_handles_labels()

        handles = np.array([handles1 + handles2])[0]
        labels = np.array([labels1 + labels2])[0]

        handles = handles[self.show_legend_entries]
        labels = labels[self.show_legend_entries]

        prop = FontProperties(size=fs)

        self.figure.sca(self.axis)

        leg = self.dummy_axis.legend(handles, labels,
                                     loc=3,
                                     markerscale=ms,
                                     prop=prop,
                                     scatterpoints=1,
                                     numpoints=1)

        leg.draggable()

        for entry in leg.legendHandles:
            entry.set_linewidth(lw)

        if not (flag and len(labels) > 0):
            leg.set_visible(False)

        return 0

    def set_ticksize(self, ticksize=20):
        """
        Set the major tick-size of all axes
        """
        for tick in self.axis.xaxis.get_major_ticks():
            tick.label1.set_fontsize(ticksize)

        for tick in self.axis.yaxis.get_major_ticks():
            tick.label1.set_fontsize(ticksize)

        for tick in self.secondary_axis.yaxis.get_major_ticks():
            tick.label2.set_fontsize(ticksize)

        # If there is a colorbar, set the size of those ticks as well
        if self.cb is not None:
            self.cb.ax.tick_params(labelsize=ticksize)

        return 0

    def set_title(self, title="Title", fontsize=None):

        if fontsize is None:
            fontsize = self.get_titlesize()

        return self.axis.set_title(title, fontsize=fontsize)

    def set_xlabel(self, label="x", fontsize=None):

        if fontsize is None:
            fontsize = self.get_xaxis_labelsize()

        self.axis.set_xlabel(label, fontsize=fontsize)

        return 0

    def set_ylabel(self, label="y", fontsize=None):

        if fontsize is None:
            fontsize = self.get_yaxis_labelsize()

        self.axis.set_ylabel(label, fontsize=fontsize)

        return 0

    def set_yscale(self, scale):
        """
        Set scale of primary y_axis
        scale can be 'linear', 'log', 'symlog'
        """

        self.axis.set_yscale(scale)

        return 0

    def set_secondary_ylabel(self, label="secondary y", fontsize=None):

        if fontsize is None:
            fontsize = self.get_yaxis_labelsize()

        self.secondary_axis.set_ylabel(label, fontsize=fontsize)

        return 0

    def set_xlim(self, xmin, xmax):

        self.axis.set_xlim(xmin, xmax)
        self.set_dummy_xlim()

        return 0

    def set_ylim(self, ymin, ymax):

        self.axis.set_ylim(ymin, ymax)

        return 0

    def set_dummy_xlim(self):

        limits = self.get_limits()
        self.dummy_axis.set_xlim(limits[0], limits[1])

        return 0

    def set_secondary_ylim(self, ymin, ymax):

        self.secondary_axis.set_ylim(ymin, ymax)

        return 0

    def set_aspect(self, aspect='equal'):

        self.axis.set_aspect(aspect)
        # self.secondary_axis.set_aspect(aspect)
        # self.dummy_axis.set_aspect(aspect)

        return 0

    def set_linewidth(self, lw=1):
        """
        """
        plots1, dummy = self.axis.get_legend_handles_labels()
        plots2, dummy = self.secondary_axis.get_legend_handles_labels()

        for plot in plots1 + plots2:
            plot.set_lw(lw)

        return 0

    def set_scatter_markersize(self, artist, new_ms, axis=1):
        """
        """
        prop = ArtistInspector(artist).properties()

        if len(prop["edgecolors"]) > 0:

            edgecolor = prop["edgecolors"][0]

        else:

            edgecolor = []

        facecolor = prop["facecolors"][0]
        # linewidth = prop["linewidths"][0]
        label = prop["label"]

        x = np.array(prop["offsets"])[:, 0]
        y = np.array(prop["offsets"])[:, 1]

        artist.remove()

        if axis == 1:

            self.axis.scatter(x, y, s=new_ms, c=facecolor, edgecolor=edgecolor, label=label)

        elif axis == 2:

            self.secondary_axis.scatter(x, y, s=new_ms, c=facecolor, edgecolor=edgecolor, label=label)

        return 0

    def set_markersize(self, ms=18):
        """
        """
        plots1, dummy = self.axis.get_legend_handles_labels()
        plots2, dummy = self.secondary_axis.get_legend_handles_labels()

        for plot in plots1:

            if type(plot) == matplotlib.lines.Line2D:

                plot.set_ms(ms)

            else:

                self.set_scatter_markersize(plot, ms, axis=1)

        for plot in plots2:

            if type(plot) == matplotlib.lines.Line2D:

                plot.set_ms(ms)

            else:

                self.set_scatter_markersize(plot, ms, axis=2)

        return 0

    def set_autoscale(self, autoscale=True, axis="both", tight=None):

        self.axis.autoscale(autoscale, axis, tight)

        return 0

    def set_secondary_autoscale(self, autoscale=True, axis="y", tight=None):

        self.secondary_axis.autoscale(autoscale, axis, tight)

        return 0

    def set_xtime(self):
        """
        Assumes the x axis values are given in s and changes the axis labels to
        have 5 tickmarks in hh:mm:ss format rather than in s.
        """

        xmin, xmax, ymin, ymax, y_sec_min, y_sec_max = self.get_limits()

        xvalues = np.linspace(xmin, xmax, 5)

        xlabels = []

        for val in xvalues:
            xlabels.append(self.timestr(val))

        self.axis.set_xticks(xvalues)
        self.axis.set_xticklabels(xlabels)

        return 0

    def get_autoscalex(self):
        """
        """
        return self.axis.get_autoscalex_on()

    def get_autoscaley1(self):
        """
        """
        return self.axis.get_autoscaley_on()

    def get_autoscaley2(self):
        """
        """
        return self.secondary_axis.get_autoscaley_on()

    def get_scale(self):
        """
        """
        return self.axis.yaxis.get_scale()

    def get_settings(self):
        """
        Gets the settings from the figure object and stores them in the
        Settings-array
        """
        limits = self.get_limits()
        labels = self.get_axis_labels()
        legend = self.get_legend()

        if legend is None:

            legend_global_flag = False

        else:

            legend_global_flag = legend.get_visible()

        handles1, legend_labels1 = self.axis.get_legend_handles_labels()
        handles2, legend_labels2 = self.secondary_axis.get_legend_handles_labels()

        # handles = handles1 + handles2
        legend_labels = legend_labels1 + legend_labels2

        entries = {"title": self.get_title(),
                   "x_label": labels[0],
                   "x_min": limits[0],
                   "x_max": limits[1],
                   "y1_label": labels[1],
                   "y1_min": limits[2],
                   "y1_max": limits[3],
                   "y2_label": labels[2],
                   "y2_min": limits[4],
                   "y2_max": limits[5]}

        adjustments = {"major_ticks": self.get_xaxis_ticksize(),
                       "label_fs": self.get_xaxis_labelsize(),
                       "legend_ms": self.get_legend_ms(),
                       "legend_fs": self.get_legend_fs(),
                       "legend_lw": self.get_legend_lw(),
                       "title_fs": self.get_titlesize(),
                       "line_width": self.get_linewidth(),
                       "marker_size": self.get_markersize()}

        settings = {"entry": entries,
                    "adj": adjustments,
                    "legend_entries": legend_labels,
                    "legend_global_flag": legend_global_flag,
                    "legend_entries_flags": self.show_legend_entries,
                    "autoscalex": self.get_autoscalex(),
                    "autoscaley1": self.get_autoscaley1(),
                    "autoscaley2": self.get_autoscaley2(),
                    "xtime": self.xtime}

        return settings

    def get_limits(self):
        """
        returns the limits in the format [xmin, xmax, ymin, ymax]
        """

        lim = list(self.axis.get_xlim())
        lim.extend(list(self.axis.get_ylim()))
        lim.extend(list(self.secondary_axis.get_ylim()))

        return lim

    def get_title(self):
        """
        Returns the current plot title
        """

        return self.axis.get_title()

    def get_titlesize(self):
        """
        """

        return self.axis.title.get_size()

    def get_axis_labels(self):
        """
        Returns the current axis labels (x,y)
        """

        return self.axis.xaxis.get_label_text(), \
            self.axis.yaxis.get_label_text(), \
            self.secondary_axis.yaxis.get_label_text()

    def get_xaxis_labelsize(self):
        """
        Returns the current axis labels font size (x)
        """

        return self.axis.xaxis.get_label().get_size()

    def get_yaxis_labelsize(self):
        """
        Returns the current axis labels font size (y)
        """

        return self.axis.yaxis.get_label().get_size()

    def get_xaxis_ticksize(self):
        """
        Get the major tick-size of xaxis
        """

        return self.axis.xaxis.get_major_ticks()[0].label1.get_fontsize()

    def get_yaxis_ticksize(self):
        """
        Get the major tick-size of yaxis
        """

        return self.axis.yaxis.get_major_ticks()[0].label1.get_fontsize()

    def get_legend(self):
        """
        """

        return self.dummy_axis.get_legend()

    def get_legend_ms(self):
        """
        Only line2D plots and scatter plots (as collections) are considered for
        the legend
        """
        legend = self.get_legend()

        ms = 1.0

        if legend is not None and len(legend.get_texts()) > 0:

            for handle in legend.legendHandles:

                if type(handle) == PathCollection:
                    ms = np.sqrt(handle.get_sizes()[0] / self.get_markersize())

        return ms

    def get_legend_lw(self):
        """
        Only line2D plots and scatter plots (as collections) are considered for
        the legend
        """
        legend = self.get_legend()

        lw = 1.0

        if legend is not None and len(legend.get_texts()) > 0:

            for handle in legend.legendHandles:

                if type(handle) == matplotlib.lines.Line2D:
                    lw = handle.get_linewidth()

        return lw

    def get_legend_fs(self):
        """
        """
        legend = self.get_legend()

        if legend is not None and len(legend.get_texts()) > 0:

            legend_fs = legend.get_texts()[0].get_size()
            return legend_fs

        else:

            return 18

    def get_markersize(self):
        """
        """

        plots, dummy = self.axis.get_legend_handles_labels()

        if len(plots) == 0:
            return self.main_plot_settings["adj"]["marker_size"]

        plot = plots[0]

        if type(plot) == matplotlib.lines.Line2D:

            ms = plot.get_ms()

        else:

            ms = plot.get_sizes()[0]

        return ms

    def get_linewidth(self):
        """
        """
        plots, dummy = self.axis.get_legend_handles_labels()

        if len(plots) == 0:
            return self.main_plot_settings["adj"]["line_width"]

        plot = plots[0]

        if type(plot) == matplotlib.lines.Line2D:

            lw = plot.get_linewidth()

        else:

            lw = plot.get_linewidth()[0]

        return lw

    def draw_idle(self):

        self.canvas.draw_idle()

        if self.xtime:
            self.set_xtime()

        return 0

    def clear(self):

        self.axis.clear()
        self.secondary_axis.clear()
        self.hide_secondary_axis()

        return 0

    def __init__(self, main_window=None, nbp=0):

        Gtk.VBox.__init__(self)

        self._main_window = main_window

        self.figure = Figure(facecolor='#474747', edgecolor='#474747')
        self.canvas = FigureCanvasGTK3Agg(self.figure)  # a Gtk.DrawingArea

        # --- Flag for primary axis time format
        self.xtime = False

        # --- The primary axis
        self.axis = self.figure.add_subplot(111)

        # --- The secondary y-axis (hidden by default)
        self.secondary_axis = self.axis.twinx()
        self.secondary_axis.set_frame_on(False)
        self.show_secondary_axis = True
        self.hide_secondary_axis()

        # --- A dummy secondary x axis for mousover events
        self.dummy_axis = self.axis.twiny()
        self.dummy_axis.set_frame_on(False)
        self.dummy_axis.set_xticks([])
        self.set_dummy_xlim()

        self.cb = None
        self.nbp = nbp

        # Init some globals for getting the mouse position in the client
        # application
        self.inaxes = False
        self.mouse_x = 0
        self.mouse_y = 0

        hbox = Gtk.HBox(spacing=4)

        self.r_label = Gtk.Label("0.0 | ")
        self.z_label = Gtk.Label("| 0.0")
        self.spacer = Gtk.Label("|")
        self.reset_b = Gtk.Button("Reset Zoom")
        self.settings_b = Gtk.Button("Settings")

        self.pack_start(self.canvas, True, True, 0.0)
        self.pack_start(hbox, False, False, 0.0)
        hbox.pack_end(self.r_label, False, False, 0.0)
        hbox.pack_end(self.spacer, False, False, 0.0)
        hbox.pack_end(self.z_label, False, False, 0.0)
        hbox.pack_end(self.reset_b, False, False, 0.0)
        hbox.pack_start(self.settings_b, False, False, 0.0)

        # --- Variables for zooming --- #
        self.holddown = False
        self.zoom_flag = False
        self.box_start = None
        self.box_start_px = None
        self.box_end = None
        self.zoom_box = None
        self.zoomlimits = [0., 1., 0., 1.]
        self.maxlimits = []
        self.autoscale_old = [True, True]

        # --- Container for flags to show certain legend entries -- #
        self.show_legend_entries = np.array([], 'bool')

        # --- Connections --- #
        self.canvas.mpl_connect("motion_notify_event", self.mouse_move)
        self.canvas.mpl_connect("scroll_event", self.mouse_scroll)
        self.canvas.mpl_connect("button_press_event", self.button_pressed)
        self.canvas.mpl_connect("button_release_event", self.button_released)
        self.reset_b.connect("clicked", self.zoom_out)
        self.settings_b.connect("clicked", self.change_settings)
        self.canvas.mpl_connect("axes_enter_event", self.axes_enter_callback)
        self.canvas.mpl_connect("axes_leave_event", self.axes_leave_callback)
        self.canvas.mpl_connect("draw_event", self.draw_event)

        self.set_settings(generate_default_settings())

        # Create a plotsettings dictionary as a copy of the current plot settings
        # Redundand, but a quick workaround for keeping the marker size and linewidth
        # even if settings were loaded without any plots
        self.main_plot_settings = self.get_settings()


class TestWindow(object):

    @staticmethod
    def destroy(*_):
        Gtk.main_quit()

    @staticmethod
    def main():
        Gtk.main()

    def test_function(self, axis):

        x = np.arange(100)
        y = x ** 2.0

        axis.plot(x, y, c="blue", label="test0", linestyle="solid")
        self.canvas.append_legend_entries_flag(True)

        axis.plot(x, y, c="green", label="test1", linestyle="solid")
        self.canvas.append_legend_entries_flag(True)

        axis.scatter(x, y * 1.2, s=5 ** 2, c="red", edgecolor="red", label="test2", marker="o")
        self.canvas.append_legend_entries_flag(True)

        return 0

    def __init__(self):
        window = Gtk.Window()
        window.set_title("MPLCanvasWrapper TestWindow")
        window.set_size_request(600, 600)
        window.connect("destroy", self.destroy)
        self.canvas = MPLCanvasWrapper(main_window=window)
        window.add(self.canvas)
        self.test_function(self.canvas.axis)
        self.canvas.set_autoscale(True)
        window.show_all()


if __name__ == '__main__':
    test = TestWindow()
    test.main()
