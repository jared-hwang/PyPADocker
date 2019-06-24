import sys
if sys.version_info.major == 3:
    from tkinter import filedialog
    from tkinter import *
    py_ver = 3
elif sys.version_info.major == 2:
    from Tkinter import *
    import tkFileDialog
    py_ver = 2
else:
    raise Exception("This version of python is not recognized!")
import os

__author__ = "Daniel Winklehner"
__doc__ = """A wrapper around tkinter filedialogs.
Note: This is mostly for legacy handling of file dialogs in my personal python scripts,
for actual GUI application development, this is probably not very useful... -DW"""


class FileDialog(object):

    def __init__(self):

        self._filename = None
        self._parent = None
        self._icon = None
        self._root = None
        self._cb_state = False

    def get_filename(self, action='open', old_path=None, icon=None, parent=None):

        self._parent = parent
        self._icon = icon
        self._root = Tk()
        self._root.withdraw()

        if py_ver == 2:
            myfiledialog = tkFileDialog
        else:
            myfiledialog = filedialog

        if action == 'open':

            self._filename = myfiledialog.askopenfilename(initialdir=old_path,
                                                          title='Open file...',
                                                          parent=self._root)

        elif action == 'save':

            self._filename = myfiledialog.asksaveasfilename(initialdir=old_path,
                                                            title='Save as...')
        elif action == 'savefig':
            self._filename = myfiledialog.asksaveasfilename(initialdir=old_path,
                                                            title='Save figure as...')
            # TODO: savefig should have a checkbutton to switch on latex rendering
        elif action == 'folder':
            self._filename = myfiledialog.askdirectory(initialdir=old_path,
                                                       title='Select folder (manually write new one in path string)')

            if self._filename is not None and not os.path.exists(self._filename) and self._filename != '':
                os.makedirs(self._filename)
            # TODO: folder should have a checkbutton to switch on going into subdirs
        else:

            raise Exception("'action' must be 'open', 'save', 'savefig', or 'folder' (got '%s')" % action)

        self._root.destroy()

        if not self._filename:
            self._filename = None
            self._cb_state = None

        if action in ["open", "save"]:

            return self._filename

        else:

            return self._filename, self._cb_state


if __name__ == '__main__':
    # Test the dialog
    fd = FileDialog()
    print(fd.get_filename(action="open", old_path="D:\Daniel\Dropbox (MIT)\Projects\RFQ Direct Injection\RFQ_Tests"))
