import gi
gi.require_version('Gtk', '3.0')  # nopep8
from gi.repository import Gtk


__author__ = "Daniel Winklehner"
__doc__ = "Simple file dialog based on the GTK+3 FileChooser"


class FileDialog(object):

    def __init__(self):

        self.filename = None
        self.parent = None
        self.checkbutton = None
        self.cb_value = None
        self.destroy_at_end = False

    def get_filename(self, action='open', old_path=None, icon=None, parent=None):

        # If there is no parent given (i.e. the Chooser is called on its own) we create an empty main window
        if parent is None:

            self.parent = Gtk.Window()
            self.destroy_at_end = True

        else:

            self.parent = parent

        if action == 'open':

            chooser = Gtk.FileChooserDialog(action=Gtk.FileChooserAction.OPEN,
                                            buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                                     Gtk.STOCK_OPEN, Gtk.ResponseType.OK))

            chooser.set_icon(icon)
            chooser.set_title('Open file...')

            if old_path is not None:
                chooser.set_current_folder(old_path)

        elif action == 'save':

            chooser = Gtk.FileChooserDialog(action=Gtk.FileChooserAction.SAVE,
                                            buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                                     Gtk.STOCK_SAVE, Gtk.ResponseType.OK))
            chooser.set_icon(icon)
            chooser.set_title('Save as...')

            if old_path is not None:
                chooser.set_current_folder(old_path)

        elif action == 'savefig':

            chooser = Gtk.FileChooserDialog(action=Gtk.FileChooserAction.SAVE,
                                            buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                                     Gtk.STOCK_SAVE, Gtk.ResponseType.OK))
            chooser.set_icon(icon)
            chooser.set_title('Save as...')

            hbox = Gtk.HBox()
            chooser.vbox.pack_end(hbox, False, False, 0)
            self.checkbutton = Gtk.CheckButton("Render in latex")
            hbox.pack_end(self.checkbutton, False, False, 0)

            if old_path is not None:
                chooser.set_current_folder(old_path)

            hbox.show()
            self.checkbutton.show()

        elif action == 'folder':

            chooser = Gtk.FileChooserDialog(action=Gtk.FileChooserAction.SELECT_FOLDER,
                                            buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                                     Gtk.STOCK_OK, Gtk.ResponseType.OK))
            chooser.set_icon(icon)
            chooser.set_title('Pick a folder...')
            hbox = Gtk.HBox()
            chooser.vbox.pack_end(hbox, False, False, 0)
            self.checkbutton = Gtk.CheckButton("Descend into subfolders")
            hbox.pack_end(self.checkbutton, False, False, 0)

            if old_path is not None:
                chooser.set_current_folder(old_path)

            hbox.show()
            self.checkbutton.show()

        else:

            raise Exception("action must be 'open', 'save' or 'folder' (got '%s')" % action)

        chooser.set_transient_for(self.parent)

        if self.filename:
            chooser.select_filename(self.filename)

        response = chooser.run()
        filename = chooser.get_filename()

        if action in ['folder', 'savefig']:
            self.cb_value = self.checkbutton.get_active()

        chooser.destroy()

        if self.destroy_at_end:
            self.parent.destroy()

        # By default, the GTK loop would wait until the process is
        # idle to process events. Now, it is very probable that file
        # I/O will be performed right after this method call and that
        # would delay hiding the dialog until I/O are done. So,
        # process pending events to hide the dialog right now.
        while Gtk.events_pending():
            Gtk.main_iteration_do(False)

        if action == "folder":

            if response == Gtk.ResponseType.OK:

                self.filename = filename

                return filename, self.cb_value

            else:

                return None, None

        elif action == "savefig":

            if response == Gtk.ResponseType.OK:

                self.filename = filename

                return filename, self.cb_value

            else:

                return None, None

        else:

            if response == Gtk.ResponseType.OK:

                self.filename = filename
                return filename

            else:

                return None


if __name__ == '__main__':
    # Test the dialog
    fd = FileDialog()
    print(fd.get_filename(action="save"))
