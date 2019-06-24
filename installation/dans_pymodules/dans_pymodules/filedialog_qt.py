print("QT5 File Dialog not yet implemented, falling back on tkinter...")
import sys
if sys.version_info >= (3, 0):
    from .filedialog_tk import *
else:
    from filedialog_tk import *
