import os
import subprocess
import time

__author__ = "Daniel Winklehner"
__doc__ = "Script to automatically compile LaTeX documents with very simple layout. " \
          "Requires a working installation of LaTeX like TeXLive or MikTeX" \
          "Make sure 'header.tex' is present in this file's directory."


# The following code to find an executable in PATH is from Jay's answer to
# http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):

    import os

    def is_exe(_fpath):

        return os.path.isfile(_fpath) and os.access(_fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    if fpath:

        if is_exe(program):

            return program

    else:

        for path in os.environ["PATH"].split(os.pathsep):

            path = path.strip('"')
            exe_file = os.path.join(path, program)

            if is_exe(exe_file):

                return exe_file

    return None


class Object:
    """
    The lowest level object.
    Figure, Text, Table, etc. inherit from here.
    """

    def __init__(self):

        self.allowed_parents = ['document', 'section', 'subsection', 'subsubsection']
        self.type = None
        self.parent = None
        self.next_object = None
        self.prev_object = None

    def append_to_section(self, section):

        if self.parent is None:

            if section.type in self.allowed_parents:

                self.parent = section
                section.append_child(self)

            return 0

        else:

            print("Already part of a section, doing nothing")
            return 1


class TextObject(Object):
    """
    A snippet of text to go in the LaTeX document
    """

    def __init__(self, text=None):

        Object.__init__(self)

        self.type = 'text'
        self.text = text

    def set_text(self, text):
        self.text = text

    def get_text(self):
        return self.text

    def get_latex_snippet(self):
        return self.text


class FigureObject(Object):
    """
    A snippet of text to go in the LaTeX document
    """

    def __init__(self, filename=None, caption=None):

        Object.__init__(self)

        self.type = 'figure'
        self.caption = caption
        self.filename = filename.replace('\\', '/').replace('_', '\string_')
        self.width = 1.0
        self.landscape = False

    def get_landscape(self):
        return self.landscape

    def get_caption(self):
        return self.caption

    def get_filename(self):
        return self.filename

    def get_latex_snippet(self):

        if self.filename is None:

            return ""

        else:
            if self.landscape:
                snippet = """\\begin{sidewaysfigure}
\centering
    \includegraphics[width=%.2f\\textwidth]{%s}""" % (self.width, self.filename)
            else:
                snippet = """\\begin{figure}
\centering
    \includegraphics[width=%.2f\\textwidth]{%s}""" % (self.width, self.filename)

            if self.caption is not None:
                snippet += """
    \caption{%s}""" % self.caption

            if self.landscape:
                snippet += """
    \label{fig:%s}
\end{sidewaysfigure}\n""" % os.path.split(self.filename)[1]

            else:
                snippet += """
    \label{fig:%s}
\end{figure}\n""" % os.path.split(self.filename)[1]

            return snippet

    def set_landscape(self, landscape=False):
        self.landscape = landscape

    def set_caption(self, caption):
        self.caption = caption

    def set_filename(self, filename):
        self.filename = filename.replace('\\', '/').replace('_', '\string_')

    def set_width(self, width):
        """
        Set the relative width ofthe image to the textwidth
        LaTeX will multiply with textwidth as in: 0.5*textwidth
        """
        self.width = width


class SectionObject(Object):
    """
    The Section class contains everything necessary to make a section in the LaTeX document.
    Subsections and subsubsections inherit from the Section().
    """

    def __init__(self, heading='Heading'):

        Object.__init__(self)

        self.allowed_parents = ['document']
        self.type = 'section'
        self.heading = heading
        self.children = []

    def append_child(self, child):

        if len(self.children) > 0:

            child.prev_object = self.children[-1]
            self.children[-1].next_object = child

        self.children.append(child)


class SubSectionObject(Object):
    """
    The Section class contains everything necessary to make a section in the LaTeX document.
    Subsections and subsubsections inherit from the Section().
    """

    def __init__(self, heading='Heading'):

        Object.__init__(self)

        self.allowed_parents = ['document', 'section']
        self.type = 'subsection'
        self.heading = heading
        self.children = []

    def append_child(self, child):

        if len(self.children) > 0:

            child.prev_object = self.children[-1]
            self.children[-1].next_object = child

        self.children.append(child)


class SubSubSectionObject(Object):
    """
    The Section class contains everything necessary to make a section in the LaTeX document.
    Subsections and subsubsections inherit from the Section().
    """

    def __init__(self, heading='Heading'):

        Object.__init__(self)

        self.allowed_parents = ['document', 'section', 'subsection']
        self.type = 'subsubsection'
        self.heading = heading
        self.children = []

    def append_child(self, child):

        if len(self.children) > 0:

            child.prev_object = self.children[-1]
            self.children[-1].next_object = child

        self.children.append(child)


class PyLatexDocument(SectionObject):
    """
    Class to handle generation of simple LaTeX documents for documentation of
    WARP postprocessing output
    """

    def __init__(self):
        """
        Create a new PyLatexDocument
        :return:
        """
        SectionObject.__init__(self)

        # --- Define variables
        self.type = 'document'
        self.author = "John Doe"
        self.title = "Generic Title"

        self.pdflatex_exe = which("pdflatex.exe")
        if self.pdflatex_exe is None:
            print("Could not find pdflatex.exe in PATH, please specify executable by using set_pdflatex_exe(<path>)")

        self.pdfreader_exe = which("")
        self.output_path = None

        self.header = self.read_header()

        self.output_stream = None

    def compile_latex(self, texfilename, mode='pdflatex'):
        """
        Compiles the latex file
        :param texfilename:
        :param mode:
        :return pdffilename:
        """

        if self.output_path is None:

            output_path = os.path.split(texfilename)[0]

        else:

            output_path = self.output_path

        if mode == 'pdflatex':

            args = [self.pdflatex_exe, "-output-directory=%s" % output_path, texfilename]

            pdffilename = os.path.join(output_path, os.path.split(texfilename)[1]).replace('\\', '/')
            pdffilename = os.path.splitext(pdffilename)[0]+'.pdf'

            subprocess.call(args)

            return pdffilename

        else:

            return None

    def get_author(self):
        return self.author

    def get_header(self):
        return self.header

    def get_output_path(self):
        return self.output_path

    def get_pdflatex_exe(self):
        return self.pdflatex_exe

    def get_title(self):
        return self.title

    @staticmethod
    def read_header():
        """
        Read the default latex file header from file
        :return text:
        """
        with open(os.path.join(os.path.dirname(__file__), 'header.tex'), 'rb') as infile:

            text = infile.read()

        return text

    def set_author(self, author):
        self.author = author

    def set_header(self, header):
        """
        Set the latex file header containing document type, margins, newcommands, etc.
        Be very careful when doing this manually, other parts of the class might depend
        on certain user defined commands and loaded libraries!
        :param header:
        :return:
        """
        self.header = header

    def set_output_path(self, output_path):
        """
        Set the path for the PDFLaTeX additional output files (.aux, etc.)
        :param output_path:
        :return:
        """
        self.output_path = output_path

    def set_pdflatex_exe(self, pdflatex_exe):
        """
        Set the path to the pdflatex executable
        :param pdflatex_exe:
        :return:
        """
        self.pdflatex_exe = pdflatex_exe

    def set_pdfreader_exe(self, pdfreader_exe):
        """
        Set the path to the executable of a pdf reader (like Adobe Reader)
        :param pdfreader_exe:
        :return:
        """
        self.pdfreader_exe = pdfreader_exe

    def set_title(self, title):
        """
        Set the document title
        :param title:
        :return:
        """
        self.title = title

    def show_pdf(self, pdffilename):

        pdffn_adj = '"{}"'.format(pdffilename)
        print(pdffn_adj)

        if self.pdfreader_exe is None:

            os.system(pdffn_adj)

        else:

            args = [self.pdfreader_exe, pdffilename]
            subprocess.call(args)

        return 0

    def test_me(self):
        """
        Function that tests the pdflatex and pdf viewer functionalities
        :return:
        """
        # self.set_pdflatex_exe('C:/texlive/2014/bin/win32/pdflatex.exe')
        self.set_pdfreader_exe('C:\Program Files (x86)\Adobe\Acrobat DC\Acrobat\Acrobat.exe')

        figpath = os.path.join(os.path.dirname(__file__), 'vitruvian.jpg')
        texfilename = os.path.join(os.path.dirname(__file__), 'Test.tex').replace('\\', '/')

        section1 = SectionObject(heading="Section 1")
        section1.append_to_section(self)

        section2 = SectionObject(heading="Section 2")
        section2.append_to_section(self)

        section3 = SectionObject(heading="Section 3")
        section3.append_to_section(self)

        subsection1 = SubSectionObject(heading="Subsection 1")
        subsection1.append_to_section(section2)

        subsection2 = SubSectionObject(heading="Subsection 2")
        subsection2.append_to_section(section2)

        text1 = TextObject(text="First text in my LaTeX document\n")
        text1.append_to_section(section1)

        text2 = TextObject(text="Second text in my LaTeX document\n")
        text2.append_to_section(subsection1)

        text3 = TextObject(text="Third text in my LaTeX document\n")
        text3.append_to_section(subsection1)
        text3.append_to_section(section1)

        figure1 = FigureObject(filename=figpath, caption="Current Variation w/o offset")
        figure1.set_width(0.5)
        figure1.append_to_section(subsection2)

        self.write_tex_file(texfilename)

        pdffilename = self.compile_latex(texfilename, mode='pdflatex')

        if pdffilename is not None:

            time.sleep(2)
            self.show_pdf(pdffilename)

        return 0

    def update_output_stream(self):
        """
        """

        self.output_stream = self.header
        self.output_stream += """
\\title{%s}

\\author{%s}

\\begin{document}

\\maketitle

""" % (self.title, self.author)

        # --- Add the main body from the different sources here ------------------------------------------------------ #
        for level0 in self.children:

            if level0.type in ['text', 'figure']:

                self.output_stream += level0.get_latex_snippet()
                self.output_stream += "\n"

            elif level0.type in ['section', 'subsection', 'subsubsections']:

                self.output_stream += '\\%s{%s}\n\n' % (level0.type, level0.heading)

                for level1 in level0.children:

                    if level1.type in ['text', 'figure']:

                        self.output_stream += level1.get_latex_snippet()
                        self.output_stream += "\n"

                    elif level1.type in ['subsection', 'subsubsections']:

                        self.output_stream += '\\%s{%s}\n\n' % (level1.type, level1.heading)

                        for level2 in level1.children:

                            if level2.type in ['text', 'figure']:

                                self.output_stream += level2.get_latex_snippet()
                                self.output_stream += "\n"

                            elif level2.type in ['subsubsections']:

                                self.output_stream += '\\%s{%s}\n\n' % (level2.type, level2.heading)

                                for level3 in level2.children:

                                    if level3.type in ['text', 'figure']:

                                        self.output_stream += level3.get_latex_snippet()
                                        self.output_stream += "\n"

        self.output_stream += """
\\end{document}
"""

        return 0

    def write_tex_file(self, texfilename):
        """
        Write the output stream to a tex file
        :param texfilename:
        :return:
        """
        self.update_output_stream()

        with open(texfilename, 'wb') as outfile:

            outfile.write(self.output_stream)

        return 0

if __name__ == '__main__':
    # Tests
    pld = PyLatexDocument()
    pld.test_me()
