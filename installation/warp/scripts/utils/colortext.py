"""Escape codes to print text in color.
For example, to print something in red:

print textRed, 'text', textColor_Off

or

print coloredtext('text', textRed)

"""

def textColor(i,t=None):
    """Set text color to one of 256 possibilities.
    A second optional argument will change the style. Possible values are:
    0: normal
    1: bold
    2: pale color
    4: underlined
    5: flashing
    7: inverse color (swaps foreground and background)
    8: background colored text

    For example, this will change the text color to green and flashing

    print textColor(35, 5)
    """
    if t is None:
        return '\033[38;5;%dm'%i
    else:
        return '\033[%d;38;5;%dm'%(t,i)

def coloredtext(text, color):
    """Returns a string, that when printed, will show the text in the specified color."""
    return color + text + textColor_Off

textColor_Off = '\033[0m'       # Text Reset

# Regular Colors
textBlack = '\033[0;30m'        # Black
textRed = '\033[0;31m'          # Red
textGreen = '\033[0;32m'        # Green
textYellow = '\033[0;33m'       # Yellow
textBlue = '\033[0;34m'         # Blue
textPurple = '\033[0;35m'       # Purple
textCyan = '\033[0;36m'         # Cyan
textWhite = '\033[0;37m'        # White

# Bold
textBBlack = '\033[1;30m'       # Black
textBRed = '\033[1;31m'         # Red
textBGreen = '\033[1;32m'       # Green
textBYellow = '\033[1;33m'      # Yellow
textBBlue = '\033[1;34m'        # Blue
textBPurple = '\033[1;35m'      # Purple
textBCyan = '\033[1;36m'        # Cyan
textBWhite = '\033[1;37m'       # White

# Underline
textUBlack = '\033[4;30m'       # Black
textURed = '\033[4;31m'         # Red
textUGreen = '\033[4;32m'       # Green
textUYellow = '\033[4;33m'      # Yellow
textUBlue = '\033[4;34m'        # Blue
textUPurple = '\033[4;35m'      # Purple
textUCyan = '\033[4;36m'        # Cyan
textUWhite = '\033[4;37m'       # White

# Background
textOn_Black = '\033[40m'       # Black
textOn_Red = '\033[41m'         # Red
textOn_Green = '\033[42m'       # Green
textOn_Yellow = '\033[43m'      # Yellow
textOn_Blue = '\033[44m'        # Blue
textOn_Purple = '\033[45m'      # Purple
textOn_Cyan = '\033[46m'        # Cyan
textOn_White = '\033[47m'       # White

# High Intensity
textIBlack = '\033[0;90m'       # Black
textIRed = '\033[0;91m'         # Red
textIGreen = '\033[0;92m'       # Green
textIYellow = '\033[0;93m'      # Yellow
textIBlue = '\033[0;94m'        # Blue
textIPurple = '\033[0;95m'      # Purple
textICyan = '\033[0;96m'        # Cyan
textIWhite = '\033[0;97m'       # White

# Bold High Intensity
textBIBlack = '\033[1;90m'      # Black
textBIRed = '\033[1;91m'        # Red
textBIGreen = '\033[1;92m'      # Green
textBIYellow = '\033[1;93m'     # Yellow
textBIBlue = '\033[1;94m'       # Blue
textBIPurple = '\033[1;95m'     # Purple
textBICyan = '\033[1;96m'       # Cyan
textBIWhite = '\033[1;97m'      # White

# High Intensity backgrounds
textOn_IBlack = '\033[0;100m'   # Black
textOn_IRed = '\033[0;101m'     # Red
textOn_IGreen = '\033[0;102m'   # Green
textOn_IYellow = '\033[0;103m'  # Yellow
textOn_IBlue = '\033[0;104m'    # Blue
textOn_IPurple = '\033[0;105m'  # Purple
textOn_ICyan = '\033[0;106m'    # Cyan
textOn_IWhite = '\033[0;107m'   # White
