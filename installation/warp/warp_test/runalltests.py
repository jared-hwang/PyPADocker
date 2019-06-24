"""Run all of the available tests.
This will run all files that have a name with the format *_test.py.
It will print out each file name and whether all of the tests passed or not.
"""
import glob
import sys
from subprocess import Popen, PIPE

textRed = '\033[0;31m'    # Red
textGreen = '\033[0;32m'  # Green
textColor_Off = '\033[0m' # Text Reset

def colored(text, textcolor):
    return textcolor + text + textColor_Off

# --- Run each of the files independently.
for f in glob.glob('*_test.py'):
    p = Popen([sys.executable, f], stdout=PIPE, stderr=PIPE, close_fds=True, universal_newlines=True)
    serr = p.stderr.readlines()
    if serr[-1] == 'OK\n':
        print('%s %s'%(f, colored('OK', textGreen)))
    else:
        print(colored(f, textRed))
        for e in serr[:-1]:
            print(e)
        # --- Use ASCII codes to print this is red
        print('%s %s'%(colored(f, textRed), colored(serr[-1], textRed)))

