#Boa:FramePanel:panel

import wx
#from wx import *
from wx.lib.anchors import LayoutAnchors
from ..envelope import matchenv
import newstdout
import warp
import sys
import __main__

[wxID_PANEL, wxID_PANELENDMATCHOUTPUT, wxID_PANELFINALA, 
 wxID_PANELFINALALABEL, wxID_PANELFINALAP, wxID_PANELFINALAPLABEL, 
 wxID_PANELFINALB, wxID_PANELFINALBLABEL, wxID_PANELFINALBP, 
 wxID_PANELFINALBPLABEL, wxID_PANELMATCHEND, wxID_PANELMATCHINGTYPES, 
 wxID_PANELPANEL1, wxID_PANELPANEL2, wxID_PANELPERIODICLABEL, 
 wxID_PANELPLOTENDMATCH, wxID_PANELSETQUAD0, wxID_PANELSETQUAD1, 
 wxID_PANELSETQUAD2, wxID_PANELSETQUAD3, wxID_PANELUSEEMLT, wxID_PANELUSEHELE, 
 wxID_PANELUSEMMLT, wxID_PANELUSEQUAD, wxID_PANELVARYQUADS, 
] = map(lambda _init_ctrls: wx.NewId(), range(25))

class panel(wx.Panel):
    def _init_coll_MatchingTypes_Pages(self, parent):
        # generated method, don't edit

        parent.AddPage(imageId=-1, page=self.panel1, select=False,
              text='Periodic')
        parent.AddPage(imageId=-1, page=self.panel2, select=True,
              text='End Conditions')

    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Panel.__init__(self, id=wxID_PANEL, name='MatchingGUI', parent=prnt,
              pos=wx.Point(0, 0), size=wx.Size(556, 349), style=wx.TAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wx.Size(548, 315))
        wx.EVT_PAINT(self, self.OnMatchingguiPaint)

        self.MatchingTypes = wx.Notebook(id=wxID_PANELMATCHINGTYPES,
              name='MatchingTypes', parent=self, pos=wx.Point(0, 0),
              size=wx.Size(548, 315), style=0)
        self.MatchingTypes.SetToolTipString('Envelope matching')

        self.panel1 = wx.Panel(id=wxID_PANELPANEL1, name='panel1',
              parent=self.MatchingTypes, pos=wx.Point(0, 0), size=wx.Size(540,
              289), style=wx.TAB_TRAVERSAL)

        self.PeriodicLabel = wx.StaticText(id=wxID_PANELPERIODICLABEL,
              label='Matching to a periodic lattice', name='PeriodicLabel',
              parent=self.panel1, pos=wx.Point(4, 4), size=wx.Size(238, 22),
              style=0)
        self.PeriodicLabel.SetToolTipString('')
        self.PeriodicLabel.SetFont(wx.Font(14, wx.SWISS, wx.NORMAL, wx.NORMAL,
              False, ''))

        self.panel2 = wx.Panel(id=wxID_PANELPANEL2, name='panel2',
              parent=self.MatchingTypes, pos=wx.Point(0, 0), size=wx.Size(540,
              289), style=wx.TAB_TRAVERSAL)

        self.FinalaLabel = wx.StaticText(id=wxID_PANELFINALALABEL,
              label='Final a', name='FinalaLabel', parent=self.panel2,
              pos=wx.Point(16, 12), size=wx.Size(39, 16), style=0)
        self.FinalaLabel.SetToolTipString('')

        self.FinalbLabel = wx.StaticText(id=wxID_PANELFINALBLABEL,
              label='Final b', name='FinalbLabel', parent=self.panel2,
              pos=wx.Point(16, 42), size=wx.Size(39, 16), style=0)
        self.FinalbLabel.SetToolTipString('')

        self.Finala = wx.TextCtrl(id=wxID_PANELFINALA, name='Finala',
              parent=self.panel2, pos=wx.Point(70, 4), size=wx.Size(80, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='afinal')
        self.Finala.SetToolTipString('Final value of a to match to')
        wx.EVT_TEXT_ENTER(self.Finala, wxID_PANELFINALA, self.OnFinalaTextEnter)

        self.Finalb = wx.TextCtrl(id=wxID_PANELFINALB, name='Finalb',
              parent=self.panel2, pos=wx.Point(70, 34), size=wx.Size(80, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='bfinal')
        self.Finalb.SetToolTipString('Final value of b to match to')
        wx.EVT_TEXT_ENTER(self.Finalb, wxID_PANELFINALB, self.OnFinalbTextEnter)

        self.FinalapLabel = wx.StaticText(id=wxID_PANELFINALAPLABEL,
              label="Final a'", name='FinalapLabel', parent=self.panel2,
              pos=wx.Point(16, 72), size=wx.Size(42, 16), style=0)
        self.FinalapLabel.SetToolTipString('')

        self.FinalbpLabel = wx.StaticText(id=wxID_PANELFINALBPLABEL,
              label="Final b'", name='FinalbpLabel', parent=self.panel2,
              pos=wx.Point(16, 102), size=wx.Size(42, 16), style=0)
        self.FinalbpLabel.SetToolTipString('')

        self.Finalap = wx.TextCtrl(id=wxID_PANELFINALAP, name='Finalap',
              parent=self.panel2, pos=wx.Point(70, 64), size=wx.Size(80, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='apfinal')
        self.Finalap.SetToolTipString("Final value of a' to match to")
        wx.EVT_TEXT_ENTER(self.Finalap, wxID_PANELFINALAP, self.OnFinalapTextEnter)

        self.Finalbp = wx.TextCtrl(id=wxID_PANELFINALBP, name='Finalbp',
              parent=self.panel2, pos=wx.Point(70, 94), size=wx.Size(80, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='bpfinal')
        self.Finalbp.SetToolTipString("Final value of b' to match to")
        wx.EVT_TEXT_ENTER(self.Finalbp, wxID_PANELFINALBP, self.OnFinalbpTextEnter)

        self.VaryQuads = wx.StaticText(id=wxID_PANELVARYQUADS,
              label='Quads to vary', name='VaryQuads', parent=self.panel2,
              pos=wx.Point(12, 120), size=wx.Size(80, 16), style=0)

        self.SetQuad0 = wx.TextCtrl(id=wxID_PANELSETQUAD0, name='SetQuad0',
              parent=self.panel2, pos=wx.Point(8, 140), size=wx.Size(30, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        wx.EVT_TEXT_ENTER(self.SetQuad0, wxID_PANELSETQUAD0,
              self.OnSetquad0TextEnter)

        self.SetQuad1 = wx.TextCtrl(id=wxID_PANELSETQUAD1, name='SetQuad1',
              parent=self.panel2, pos=wx.Point(42, 140), size=wx.Size(30, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        wx.EVT_TEXT_ENTER(self.SetQuad1, wxID_PANELSETQUAD1,
              self.OnSetquad1TextEnter)

        self.SetQuad2 = wx.TextCtrl(id=wxID_PANELSETQUAD2, name='SetQuad2',
              parent=self.panel2, pos=wx.Point(76, 140), size=wx.Size(30, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        self.SetQuad2.SetToolTipString('')
        wx.EVT_TEXT_ENTER(self.SetQuad2, wxID_PANELSETQUAD2,
              self.OnSetquad2TextEnter)

        self.SetQuad3 = wx.TextCtrl(id=wxID_PANELSETQUAD3, name='SetQuad3',
              parent=self.panel2, pos=wx.Point(110, 140), size=wx.Size(30, 22),
              style=wx.TAB_TRAVERSAL | wx.TE_PROCESS_TAB | wx.TE_PROCESS_ENTER,
              value='')
        self.SetQuad3.SetToolTipString('')
        wx.EVT_TEXT_ENTER(self.SetQuad3, wxID_PANELSETQUAD3,
              self.OnSetquad3TextEnter)

        self.UseQuad = wx.RadioButton(id=wxID_PANELUSEQUAD, label='Quad',
              name='UseQuad', parent=self.panel2, pos=wx.Point(16, 170),
              size=wx.Size(94, 24), style=wx.TAB_TRAVERSAL)
        self.UseQuad.SetValue(True)
        self.UseQuad.SetToolTipString('Use hard edged quadrupoles')
        wx.EVT_RADIOBUTTON(self.UseQuad, wxID_PANELUSEQUAD,
              self.OnUsequadRadiobutton)

        self.UseHele = wx.RadioButton(id=wxID_PANELUSEHELE, label='Hele',
              name='UseHele', parent=self.panel2, pos=wx.Point(16, 188),
              size=wx.Size(94, 24), style=wx.TAB_TRAVERSAL)
        self.UseHele.SetValue(False)
        self.UseHele.SetToolTipString('Use hard edged elements')
        wx.EVT_RADIOBUTTON(self.UseHele, wxID_PANELUSEHELE,
              self.OnUseheleRadiobutton)

        self.UseEmlt = wx.RadioButton(id=wxID_PANELUSEEMLT, label='Emlt',
              name='UseEmlt', parent=self.panel2, pos=wx.Point(86, 170),
              size=wx.Size(94, 24), style=wx.TAB_TRAVERSAL)
        self.UseEmlt.SetValue(False)
        self.UseEmlt.SetToolTipString('Use axially varying electric element')
        wx.EVT_RADIOBUTTON(self.UseEmlt, wxID_PANELUSEEMLT,
              self.OnUseemltRadiobutton)

        self.UseMmlt = wx.RadioButton(id=wxID_PANELUSEMMLT, label='Mmlt',
              name='UseMmlt', parent=self.panel2, pos=wx.Point(86, 188),
              size=wx.Size(94, 24), style=wx.TAB_TRAVERSAL)
        self.UseMmlt.SetValue(False)
        self.UseMmlt.SetToolTipString('Use axially varying magnetic elements')
        wx.EVT_RADIOBUTTON(self.UseMmlt, wxID_PANELUSEMMLT,
              self.OnUsemmltRadiobutton)

        self.MatchEnd = wx.Button(id=wxID_PANELMATCHEND, label='Match',
              name='MatchEnd', parent=self.panel2, pos=wx.Point(12, 236),
              size=wx.Size(80, 22), style=wx.TAB_TRAVERSAL)
        self.MatchEnd.SetToolTipString('')
        wx.EVT_BUTTON(self.MatchEnd, wxID_PANELMATCHEND, self.OnMatchendButton)

        self.PlotEndMatch = wx.CheckBox(id=wxID_PANELPLOTENDMATCH,
              label='Plot envelope after match', name='PlotEndMatch',
              parent=self.panel2, pos=wx.Point(12, 258), size=wx.Size(168, 24),
              style=wx.TAB_TRAVERSAL)
        self.PlotEndMatch.SetValue(True)
        wx.EVT_CHECKBOX(self.PlotEndMatch, wxID_PANELPLOTENDMATCH,
              self.OnPlotendmatchCheckbox)

        self.EndMatchOutput = wx.TextCtrl(id=wxID_PANELENDMATCHOUTPUT,
              name='EndMatchOutput', parent=self.panel2, pos=wx.Point(180, 24),
              size=wx.Size(360, 250), style=wx.TE_READONLY | wx.TE_MULTILINE,
              value='')
        self.EndMatchOutput.SetToolTipString('')

        self._init_coll_MatchingTypes_Pages(self.MatchingTypes)

    def __init__(self, parent, id=0, pos=0, size=0, style=0, name=0):
        self._init_ctrls(parent)
        self.Move(wx.Point(0,0))
        self.plotafterendmatch = self.PlotEndMatch.GetValue()
        self.endmatchquads = [None,None,None,None]
        self.usequad = self.UseQuad.GetValue()
        self.usehele = self.UseHele.GetValue()
        self.useemlt = self.UseEmlt.GetValue()
        self.usemmlt = self.UseMmlt.GetValue()

    def OnFinalaTextEnter(self, event,defval=''):
        try:
            self.afinal = eval(self.Finala.GetValue(),__main__.__dict__)
            self.Finala.SetValue(str(self.afinal))
        except:
            self.Finala.SetValue(defval)

    def OnFinalbTextEnter(self, event,defval=''):
        try:
            self.bfinal = eval(self.Finalb.GetValue(),__main__.__dict__)
            self.Finalb.SetValue(str(self.bfinal))
        except:
            self.Finalb.SetValue(defval)

    def OnFinalapTextEnter(self, event,defval=''):
        try:
            self.apfinal = eval(self.Finalap.GetValue(),__main__.__dict__)
            self.Finalap.SetValue(str(self.apfinal))
        except:
            self.Finalap.SetValue(defval)

    def OnFinalbpTextEnter(self, event,defval=''):
        try:
            self.bpfinal = eval(self.Finalbp.GetValue(),__main__.__dict__)
            self.Finalbp.SetValue(str(self.bpfinal))
        except:
            self.Finalbp.SetValue(defval)

    def OnMatchendButton(self, event):
        savestdout = sys.stdout
        sys.stdout = newstdout.newstdout(self.EndMatchOutput)
        matchenv.matchenv(self.endmatchquads,self.afinal,self.bfinal,
                          self.apfinal,self.bpfinal,
                          usequad=self.usequad,usehele=self.usehele,
                          useemlt=self.useemlt,usemmlt=self.usemmlt)
        sys.stdout = savestdout
        if self.plotafterendmatch:
            warp.fma()
            warp.penv()

    def OnPlotendmatchCheckbox(self, event):
        self.plotafterendmatch = self.PlotEndMatch.GetValue()

    def OnSetquad0TextEnter(self, event):
        self.endmatchquads[0] = eval(self.SetQuad0.GetValue())
        if self.endmatchquads[1] is None:
            self.endmatchquads[1] = self.endmatchquads[0] + 1
            self.endmatchquads[2] = self.endmatchquads[1] + 1
            self.endmatchquads[3] = self.endmatchquads[2] + 1
            self.SetQuad1.SetValue(str(self.endmatchquads[1]))
            self.SetQuad2.SetValue(str(self.endmatchquads[2]))
            self.SetQuad3.SetValue(str(self.endmatchquads[3]))

    def OnSetquad1TextEnter(self, event):
        self.endmatchquads[1] = eval(self.SetQuad1.GetValue())

    def OnSetquad2TextEnter(self, event):
        self.endmatchquads[2] = eval(self.SetQuad2.GetValue())

    def OnSetquad3TextEnter(self, event):
        self.endmatchquads[3] = eval(self.SetQuad3.GetValue())

    def OnMatchingguiPaint(self, event):
        warp.package('env')
        self.OnFinalaTextEnter(None,'afinal')
        self.OnFinalbTextEnter(None,'bfinal')
        self.OnFinalapTextEnter(None,'apfinal')
        self.OnFinalbpTextEnter(None,'bpfinal')
        event.Skip()

    def OnUseheleRadiobutton(self, event):
        self.usehele = self.UseHele.GetValue()
        event.Skip()

    def OnUsequadRadiobutton(self, event):
        self.usequad = self.UseQuad.GetValue()
        event.Skip()

    def OnUseemltRadiobutton(self, event):
        self.useemlt = self.UseEmlt.GetValue()
        event.Skip()

    def OnUsemmltRadiobutton(self, event):
        self.usemmlt = self.UseMmlt.GetValue()
        event.Skip()
