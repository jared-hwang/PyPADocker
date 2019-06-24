#Boa:FramePanel:panel

from wx import *
from wx.lib.anchors import LayoutAnchors
from wx.grid import *
from warp import *
from Opyndx import *
import __main__

[wxID_PANEL, wxID_PANELFILENAME, wxID_PANELFORMAT, wxID_PANELMODE, 
 wxID_PANELNEWWINDOW, wxID_PANELRENDERING, wxID_PANELSAVE, 
 wxID_PANELSTATICBOX1, wxID_PANELSTATICBOX2, wxID_PANELSTATICBOX3, 
 wxID_PANELSTATICTEXT1, wxID_PANELSTATICTEXT10, wxID_PANELSTATICTEXT11, 
 wxID_PANELSTATICTEXT2, wxID_PANELSTATICTEXT3, wxID_PANELSTATICTEXT4, 
 wxID_PANELSTATICTEXT5, wxID_PANELSTATICTEXT6, wxID_PANELSTATICTEXT7, 
 wxID_PANELSTATICTEXT8, wxID_PANELSTATICTEXT9, wxID_PANELWINDOW, 
 wxID_PANELXSCALE, wxID_PANELXSCALE_MAX, wxID_PANELXSCALE_MIN, 
 wxID_PANELXSCALE_VALUE, wxID_PANELYSCALE, wxID_PANELYSCALE_MAX, 
 wxID_PANELYSCALE_MIN, wxID_PANELYSCALE_VALUE, wxID_PANELZSCALE, 
 wxID_PANELZSCALE_MAX, wxID_PANELZSCALE_MIN, wxID_PANELZSCALE_VALUE, 
] = map(lambda _init_ctrls: wx.NewId(), range(34))

class panel(wx.Panel):
    def _init_utils(self):
        # generated method, don't edit
        pass

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Panel.__init__(self, id=wxID_PANEL, name='OpenDX', parent=prnt,
              pos=wx.Point(433, 233), size=wx.Size(365, 246),
              style=wx.TAB_TRAVERSAL)
        self._init_utils()
        self.SetClientSize(wx.Size(357, 222))
        EVT_ENTER_WINDOW(self, self.OnPanelEnterWindow)

        self.staticText1 = wx.StaticText(id=wxID_PANELSTATICTEXT1, label='Mode',
              name='staticText1', parent=self, pos=wx.Point(16, 44),
              size=wx.Size(44, 16), style=0)

        self.Mode = wx.Choice(choices=['rotate', 'translate', 'zoom'],
              id=wxID_PANELMODE, name='Mode', parent=self, pos=wx.Point(76, 40),
              size=wx.Size(90, 21), style=0, validator=wx.DefaultValidator)
        self.Mode.SetLabel('')
        self.Mode.Show(True)
        EVT_CHOICE(self.Mode, wxID_PANELMODE, self.OnModeChoice)

        self.staticText2 = wx.StaticText(id=wxID_PANELSTATICTEXT2,
              label='Filename', name='staticText2', parent=self,
              pos=wx.Point(192, 20), size=wx.Size(50, 13), style=0)

        self.Format = wx.Choice(choices=['eps', 'eps grey', 'miff', 'ps',
              'ps grey', 'rgb', 'r+g+b', 'tiff', 'yuv', 'dx', 'vrml'],
              id=wxID_PANELFORMAT, name='Format', parent=self, pos=wx.Point(244,
              40), size=wx.Size(90, 21), style=0, validator=wx.DefaultValidator)

        self.Filename = wx.TextCtrl(id=wxID_PANELFILENAME, name='Filename',
              parent=self, pos=wx.Point(244, 16), size=wx.Size(90, 21), style=0,
              value='Image')

        self.staticText3 = wx.StaticText(id=wxID_PANELSTATICTEXT3,
              label='Format', name='staticText3', parent=self, pos=wx.Point(192,
              44), size=wx.Size(40, 13), style=0)

        self.Save = wx.Button(id=wxID_PANELSAVE, label='Save', name='Save',
              parent=self, pos=wx.Point(244, 64), size=wx.Size(90, 21), style=0)
        EVT_BUTTON(self.Save, wxID_PANELSAVE, self.OnSaveButton)

        self.staticBox1 = wx.StaticBox(id=wxID_PANELSTATICBOX1, label='Display',
              name='staticBox1', parent=self, pos=wx.Point(8, 0),
              size=wx.Size(168, 96), style=0)

        self.staticText4 = wx.StaticText(id=wxID_PANELSTATICTEXT4,
              label='Rendering', name='staticText4', parent=self,
              pos=wx.Point(16, 68), size=wx.Size(55, 13), style=0)

        self.Rendering = wx.Choice(choices=['software', 'hardware'],
              id=wxID_PANELRENDERING, name='Rendering', parent=self,
              pos=wx.Point(76, 64), size=wx.Size(90, 21), style=0,
              validator=wx.DefaultValidator)
        EVT_CHOICE(self.Rendering, wxID_PANELRENDERING, self.OnRenderingChoice)

        self.staticBox2 = wx.StaticBox(id=wxID_PANELSTATICBOX2, label='Output',
              name='staticBox2', parent=self, pos=wx.Point(184, 0),
              size=wx.Size(164, 96), style=0)

        self.staticBox3 = wx.StaticBox(id=wxID_PANELSTATICBOX3, label='Scale',
              name='staticBox3', parent=self, pos=wx.Point(8, 96),
              size=wx.Size(340, 96), style=0)

        self.staticText5 = wx.StaticText(id=wxID_PANELSTATICTEXT5, label='X',
              name='staticText5', parent=self, pos=wx.Point(16, 126),
              size=wx.Size(16, 13), style=0)

        self.staticText6 = wx.StaticText(id=wxID_PANELSTATICTEXT6, label='Y',
              name='staticText6', parent=self, pos=wx.Point(16, 148),
              size=wx.Size(16, 13), style=0)

        self.staticText7 = wx.StaticText(id=wxID_PANELSTATICTEXT7, label='Z',
              name='staticText7', parent=self, pos=wx.Point(16, 170),
              size=wx.Size(16, 13), style=0)

        self.staticText8 = wx.StaticText(id=wxID_PANELSTATICTEXT8, label='value',
              name='staticText8', parent=self, pos=wx.Point(40, 108),
              size=wx.Size(28, 13), style=0)

        self.staticText9 = wx.StaticText(id=wxID_PANELSTATICTEXT9, label='min',
              name='staticText9', parent=self, pos=wx.Point(116, 108),
              size=wx.Size(18, 13), style=0)

        self.staticText10 = wx.StaticText(id=wxID_PANELSTATICTEXT10, label='max',
              name='staticText10', parent=self, pos=wx.Point(308, 108),
              size=wx.Size(21, 13), style=0)

        self.Xscale = wx.Slider(id=wxID_PANELXSCALE, maxValue=1000, minValue=0,
              name='Xscale', parent=self, pos=wx.Point(150, 122),
              size=wx.Size(139, 20), style=wx.SL_HORIZONTAL,
              validator=wx.DefaultValidator, value=1)
        self.Xscale.SetLabel('')
        EVT_SCROLL(self.Xscale, self.OnXscaleScroll)

        self.Xscale_min = wx.TextCtrl(id=wxID_PANELXSCALE_MIN, name='Xscale_min',
              parent=self, pos=wx.Point(104, 122), size=wx.Size(40, 21), style=0,
              value='0')

        self.Xscale_Value = wx.TextCtrl(id=wxID_PANELXSCALE_VALUE,
              name='XScale_Value', parent=self, pos=wx.Point(32, 122),
              size=wx.Size(48, 21), style=0, value='1.')
        EVT_TEXT_ENTER(self.Xscale_Value, wxID_PANELXSCALE_VALUE,
              self.OnXscale_valueTextEnter)

        self.Xscale_max = wx.TextCtrl(id=wxID_PANELXSCALE_MAX, name='XScale_max',
              parent=self, pos=wx.Point(296, 122), size=wx.Size(40, 21), style=0,
              value='10.')

        self.Yscale = wx.Slider(id=wxID_PANELYSCALE, maxValue=1000, minValue=0,
              name='Yscale', parent=self, pos=wx.Point(150, 144),
              size=wx.Size(139, 20), style=wx.SL_HORIZONTAL,
              validator=wx.DefaultValidator, value=1)
        self.Yscale.SetLabel('')
        EVT_SCROLL(self.Yscale, self.OnYscaleScroll)

        self.Yscale_min = wx.TextCtrl(id=wxID_PANELYSCALE_MIN, name='Yscale_min',
              parent=self, pos=wx.Point(104, 144), size=wx.Size(40, 21), style=0,
              value='0')

        self.Yscale_Value = wx.TextCtrl(id=wxID_PANELYSCALE_VALUE,
              name='YScale_Value', parent=self, pos=wx.Point(32, 144),
              size=wx.Size(48, 21), style=0, value='1.')
        EVT_TEXT_ENTER(self.Yscale_Value, wxID_PANELYSCALE_VALUE,
              self.OnYscale_valueTextEnter)

        self.Yscale_max = wx.TextCtrl(id=wxID_PANELYSCALE_MAX, name='YScale_max',
              parent=self, pos=wx.Point(296, 144), size=wx.Size(40, 21), style=0,
              value='10.')

        self.Zscale = wx.Slider(id=wxID_PANELZSCALE, maxValue=1000, minValue=0,
              name='Zscale', parent=self, pos=wx.Point(150, 166),
              size=wx.Size(139, 20), style=wx.SL_HORIZONTAL,
              validator=wx.DefaultValidator, value=1)
        self.Zscale.SetLabel('')
        EVT_SCROLL(self.Zscale, self.OnZscaleScroll)

        self.Zscale_min = wx.TextCtrl(id=wxID_PANELZSCALE_MIN, name='Zscale_min',
              parent=self, pos=wx.Point(104, 166), size=wx.Size(40, 21), style=0,
              value='0')

        self.Zscale_Value = wx.TextCtrl(id=wxID_PANELZSCALE_VALUE,
              name='ZScale_Value', parent=self, pos=wx.Point(32, 166),
              size=wx.Size(48, 21), style=0, value='1.')
        EVT_TEXT_ENTER(self.Zscale_Value, wxID_PANELZSCALE_VALUE,
              self.OnZscale_valueTextEnter)

        self.Zscale_max = wx.TextCtrl(id=wxID_PANELZSCALE_MAX, name='ZScale_max',
              parent=self, pos=wx.Point(296, 166), size=wx.Size(40, 21), style=0,
              value='10.')

        self.staticText11 = wx.StaticText(id=wxID_PANELSTATICTEXT11,
              label='Window', name='staticText11', parent=self, pos=wx.Point(16,
              20), size=wx.Size(42, 13), style=0)

        self.Window = wx.Choice(choices=['All'], id=wxID_PANELWINDOW,
              name='Window', parent=self, pos=wx.Point(76, 16), size=wx.Size(90,
              21), style=0, validator=wx.DefaultValidator)
        EVT_CHOICE(self.Window, wxID_PANELWINDOW, self.OnWindowChoice)

        self.NewWindow = wx.CheckBox(id=wxID_PANELNEWWINDOW,
              label='Next plot in new window', name='NewWindow', parent=self,
              pos=wx.Point(16, 200), size=wx.Size(160, 16), style=0)
        self.NewWindow.SetValue(False)
        EVT_CHECKBOX(self.NewWindow, wxID_PANELNEWWINDOW,
              self.OnNewwindowCheckbox)

    def __init__(self, parent):
        self._init_ctrls(parent)
        self.Move(wx.Point(0,0))

    def OnPanelEnterWindow(self, event):
        try:
            self.Mode.SetSelection(__main__.dxwindow.interactor)
            self.Rendering.SetSelection(__main__.dxwindow.l_hardware_acceleration)
            self.NewWindow.SetValue(__main__.l_dxnewwindow)
            dxwindows = __main__.DXWindows.copy()
            toremove=[]
            for i in range(1,self.Window.GetCount()):
                win_name = self.Window.GetString(i)
                if dxwindows.has_key(win_name):
                    dxwindows.pop(win_name)
                else:
                    toremove+=[i]
            toremove.reverse()
            for i in toremove:
                self.Window.Delete(i)
            for k in dxwindows.iterkeys():
                self.Window.Append(k)       
            if __main__.l_dxupdate_all_windows:
                self.Window.SetStringSelection('All')
            else:
                self.Window.SetStringSelection(__main__.dxwindow.name)
            event.Skip()
        except:
            event.Skip()
    
    def OnModeChoice(self, event):
        __main__.dxwindow.interactor = self.Mode.GetSelection()
        event.Skip()

    def OnSaveButton(self, event):
        self.filename = str(self.Filename.GetValue())
        self.format = str(self.Format.GetStringSelection())
        if self.format=='jpeg':self.format='ImageMagick supported format gamma=2.2 compression=JPEG quality=90'
        DXWriteImage(self.filename,__main__.dxwindow.dxobject,__main__.dxwindow.dcamera,None,self.format)
        event.Skip()

    def OnRenderingChoice(self, event):
        if __main__.l_dxupdate_all_windows:
            windows=__main__.DXWindows.values()
        else:
            windows=[__main__.dxwindow]
        rendering=str(self.Rendering.GetStringSelection())
        for w in windows:
            DXRendering(w,rendering)
        event.Skip()

    def OnXscaleScroll(self, event):
        xmin = float(self.Xscale_min.GetValue())
        xmax = float(self.Xscale_max.GetValue())
        __main__.dxwindow.dxscale[0]=xmin+(xmax-xmin)*self.Xscale.GetValue()/1000.
        self.Xscale_Value.SetValue('%g'%__main__.dxwindow.dxscale[0])
        __main__.dxwindow.dxobject=DXScale(__main__.dxwindow.dxobject_init,__main__.dxwindow.dxscale)
        __main__.dxwindow.l_dxrescaled=1
        __main__.l_dxforceupdate=1
        event.Skip()

    def OnYscaleScroll(self, event):
        ymin = float(self.Yscale_min.GetValue())
        ymax = float(self.Yscale_max.GetValue())
        __main__.dxwindow.dxscale[1]=ymin+(ymax-ymin)*self.Yscale.GetValue()/1000.
        self.Yscale_Value.SetValue('%g'%__main__.dxwindow.dxscale[1])
        __main__.dxwindow.dxobject=DXScale(__main__.dxwindow.dxobject_init,__main__.dxwindow.dxscale)
        __main__.dxwindow.l_dxrescaled=1
        __main__.l_dxforceupdate=1
        event.Skip()

    def OnZscaleScroll(self, event):
        zmin = float(self.Zscale_min.GetValue())
        zmax = float(self.Zscale_max.GetValue())
        __main__.dxwindow.dxscale[2]=zmin+(zmax-zmin)*self.Zscale.GetValue()/1000.
        self.Zscale_Value.SetValue('%g'%__main__.dxwindow.dxscale[2])
        __main__.dxwindow.dxobject=DXScale(__main__.dxwindow.dxobject_init,__main__.dxwindow.dxscale)
        __main__.dxwindow.l_dxrescaled=1
        __main__.l_dxforceupdate=1
        event.Skip()

    def OnXscale_valueTextEnter(self, event):
        __main__.dxwindow.dxscale[0]=float(self.Xscale_Value.GetValue())
        __main__.dxwindow.dxobject=DXScale(__main__.dxwindow.dxobject_init,__main__.dxwindow.dxscale)
        __main__.dxwindow.l_dxrescaled=1
        __main__.l_dxforceupdate=1
        event.Skip()

    def OnYscale_valueTextEnter(self, event):
        __main__.dxwindow.dxscale[1]=float(self.Xscale_Value.GetValue())
        __main__.dxwindow.dxobject=DXScale(__main__.dxwindow.dxobject_init,__main__.dxwindow.dxscale)
        __main__.dxwindow.l_dxrescaled=1
        __main__.l_dxforceupdate=1
        event.Skip()

    def OnZscale_valueTextEnter(self, event):
        __main__.dxwindow.dxscale[2]=float(self.Xscale_Value.GetValue())
        __main__.dxwindow.dxobject=DXScale(__main__.dxwindow.dxobject_init,__main__.dxwindow.dxscale)
        __main__.dxwindow.l_dxrescaled=1
        __main__.l_dxforceupdate=1
        event.Skip()

    def OnWindowChoice(self, event):
        selection=self.Window.GetStringSelection()
        if selection=='All':
            __main__.l_dxupdate_all_windows=1
        else:
            __main__.l_dxupdate_all_windows=0
            __main__.dxwindow=__main__.DXWindows[selection]
        event.Skip()

    def OnNewwindowCheckbox(self, event):
        __main__.l_dxnewwindow=self.NewWindow.GetValue()
        event.Skip()

