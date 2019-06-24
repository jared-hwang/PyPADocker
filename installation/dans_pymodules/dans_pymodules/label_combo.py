import gi
gi.require_version('Gtk', '3.0')  # nopep8
from gi.repository import Gtk

__author__ = "Daniel Winklehner"
__doc__ = "Convenience Class for a Gtk.Label plus Gtk.ComboEntry widget"


class LabelCombo(Gtk.HBox):
    """
    Convenience Class for a Label plus ComboEntry widget
    """

    def __init__(self, mode=1, label="", combovalues=None, tooltip=None):
        """
        Type reference: mode = 1 is Label + ComboBox
                        mode = 2 is Label + ComboBoxEntry
                        mode = 3 is Label + Entry
        """
        self.mode = mode
        Gtk.HBox.__init__(self)
        combolabel = Gtk.Label(label)

        if mode == 1:
            
            self.comboboxentry = Gtk.ComboBoxText()
            
        elif mode == 2:
            
            self.comboboxentry = Gtk.ComboBoxText.new_with_entry()
            
        else:
            
            self.comboboxentry = Gtk.Entry()
            self.comboboxentry.set_width_chars(12)
            
        self.pack_start(combolabel, False, False, 2)
        self.pack_end(self.comboboxentry, False, False, 2)
        
        if combovalues is not None:

            if type(combovalues).__name__ not in ["list", "ndarray"]:
                
                combovalues = [combovalues]

            if self.mode in [1, 2]:

                for i in range(len(combovalues)):

                    if type(combovalues[i]) != 'str':

                        self.comboboxentry.append_text(str(combovalues[i]))

                    else:

                        self.comboboxentry.append_text(combovalues[i])

                self.comboboxentry.set_active(0)

            else:

                if type(combovalues[0]) == 'str':

                    self.comboboxentry.set_text(combovalues[0])
                else:

                    self.comboboxentry.set_text(str(combovalues[0]))

        if tooltip is not None:
            self.comboboxentry.set_tooltip_text(tooltip)

    def get_active(self):

        if self.mode in [1, 2]:

            return self.comboboxentry.get_active()

        else:

            print("Tried to get_active from Type 3 LabelCombo")

    def get_active_text(self):

        if self.mode == 1:

            model = self.comboboxentry.get_model()
            active = self.comboboxentry.get_active()

            if active < 0:
                return None

            return model[active][0]

        elif self.mode == 2:

            return self.comboboxentry.child.get_text()

        else:

            return self.comboboxentry.get_text()

    def append_text(self, text):

        if self.mode in [1, 2]:

            self.comboboxentry.append_text(text)

        else:

            self.comboboxentry.set_text(str(text))

    def prepend_text(self, text):

        if self.mode in [1, 2]:

            self.comboboxentry.prepend_text(text)

        else:

            self.comboboxentry.set_text(text)

    def insert_text(self, position, text):

        if self.mode in [1, 2]:

            self.comboboxentry.insert_text(position, text)

        else:

            self.comboboxentry.set_text(text)

    def remove_text(self, position):

        if self.mode in [1, 2]:

            self.comboboxentry.remove_text(position)

        else:

            self.comboboxentry.set_text("")

    def set_active(self, index):

        if self.mode in [1, 2]:
            
            self.comboboxentry.set_active(index)
