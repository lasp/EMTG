#EMTG: Evolutionary Mission Trajectory Generator
#An open-source global optimization tool for preliminary mission design
#Provided by NASA Goddard Space Flight Center
#
#Copyright (c) 2014 - 2024 United States Government as represented by the
#Administrator of the National Aeronautics and Space Administration.
#All Other Rights Reserved.
#
# Copyright (c) 2024 The Regents of the University of Colorado.
# All Other Rights Reserved.

#Licensed under the NASA Open Source License (the "License"); 
#You may not use this file except in compliance with the License. 
#You may obtain a copy of the License at:
#https://opensource.org/license/nasa1-3-php
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
#express or implied.   See the License for the specific language
#governing permissions and limitations under the License.

import OrbitElementsPanel

import wx
import platform

class DepartureElementsPanel(OrbitElementsPanel.OrbitElementsPanel):
    def __init__(self, parent, missionoptions):
        OrbitElementsPanel.OrbitElementsPanel.__init__(self, parent, missionoptions)

        #custom departure elements
        self.boxdeparture_elements = wx.StaticBox(self, -1, "Journey departure elements")
        

        self.DepartureElementsBox = wx.StaticBoxSizer(self.boxdeparture_elements, wx.VERTICAL)
        self.DepartureElementsBox.AddMany([self.elements_representation_box, self.elements_frame_box, self.elements_epoch_box, self.ElementsSizer])
        
        font = self.GetFont()
        font.SetWeight(wx.FONTWEIGHT_BOLD)
        self.boxdeparture_elements.SetFont(font)

        
        self.SetSizer(self.DepartureElementsBox)

        self.cmbelements_representation.Bind(wx.EVT_COMBOBOX, self.Changedeparture_elements_representation)
        self.cmbelements_frame.Bind(wx.EVT_COMBOBOX, self.Changedeparture_elements_frame)
        self.txtelements_reference_epoch.Bind(wx.EVT_KILL_FOCUS, self.Changedeparture_elements_reference_epoch)
        self.chkAllowJourneyFreePointToPropagate.Bind(wx.EVT_CHECKBOX, self.ChangeAllowJourneyFreePointDepartureToPropagate)
        self.chkSMA.Bind(wx.EVT_CHECKBOX,self.ChangevarySMA_departure)
        self.chkECC.Bind(wx.EVT_CHECKBOX,self.ChangevaryECC_departure)
        self.chkINC.Bind(wx.EVT_CHECKBOX,self.ChangevaryINC_departure)
        self.chkRAAN.Bind(wx.EVT_CHECKBOX,self.ChangevaryRAAN_departure)
        self.chkAOP.Bind(wx.EVT_CHECKBOX,self.ChangevaryAOP_departure)
        self.chkMA.Bind(wx.EVT_CHECKBOX,self.ChangevaryMA_departure)
        self.txtSMA.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA_departure)
        self.txtECC.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC_departure)
        self.txtINC.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC_departure)
        self.txtRAAN.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN_departure)
        self.txtAOP.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP_departure)
        self.txtMA.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA_departure)
        self.txtSMA0.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA0)
        self.txtECC0.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC0)
        self.txtINC0.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC0)
        self.txtRAAN0.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN0)
        self.txtAOP0.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP0)
        self.txtMA0.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA0)
        self.txtSMA1.Bind(wx.EVT_KILL_FOCUS,self.ChangeSMA1)
        self.txtECC1.Bind(wx.EVT_KILL_FOCUS,self.ChangeECC1)
        self.txtINC1.Bind(wx.EVT_KILL_FOCUS,self.ChangeINC1)
        self.txtRAAN1.Bind(wx.EVT_KILL_FOCUS,self.ChangeRAAN1)
        self.txtAOP1.Bind(wx.EVT_KILL_FOCUS,self.ChangeAOP1)
        self.txtMA1.Bind(wx.EVT_KILL_FOCUS,self.ChangeMA1)

    def update(self):
        self.cmbelements_representation.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_state_representation)
        self.cmbelements_frame.SetSelection(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_frame)
        self.txtelements_reference_epoch.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_reference_epoch))

        self.chkAllowJourneyFreePointToPropagate.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AllowJourneyFreePointDepartureToPropagate)
        self.chkSMA.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[0])
        self.chkECC.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[1])
        self.chkINC.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[2])
        self.chkRAAN.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[3])
        self.chkAOP.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[4])
        self.chkMA.SetValue(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[5])
        self.txtSMA.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[0]))
        self.txtECC.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[1]))
        self.txtINC.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[2]))
        self.txtRAAN.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[3]))
        self.txtAOP.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[4]))
        self.txtMA.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[5]))
        self.txtSMA0.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[0]))
        self.txtECC0.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[2]))
        self.txtINC0.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[4]))
        self.txtRAAN0.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[6]))
        self.txtAOP0.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[8]))
        self.txtMA0.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[10]))
        self.txtSMA1.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[1]))
        self.txtECC1.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[3]))
        self.txtINC1.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[5]))
        self.txtRAAN1.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[7]))
        self.txtAOP1.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[9]))
        self.txtMA1.SetValue(str(self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[11]))
        
        # Create labels for the state representations.
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_state_representation == 0:
            self.lblSMA.SetLabel("x (km)")
            self.lblECC.SetLabel("y (km)")
            self.lblINC.SetLabel("z (km)")
            self.lblRAAN.SetLabel("vx (km/s)")
            self.lblAOP.SetLabel("vy (km/s)")
            self.lblMA.SetLabel("vz (km/s)")
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_state_representation == 1:
            self.lblSMA.SetLabel("r (km)")
            self.lblECC.SetLabel("RA (degrees)")
            self.lblINC.SetLabel("DEC (degrees)")
            self.lblRAAN.SetLabel("v (km/s)")
            self.lblAOP.SetLabel("vRA (degrees)")
            self.lblMA.SetLabel("vDEC (degrees)")
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_state_representation == 2:
            self.lblSMA.SetLabel("r (km)")
            self.lblECC.SetLabel("RA (degrees)")
            self.lblINC.SetLabel("DEC (degrees)")
            self.lblRAAN.SetLabel("v (km/s)")
            self.lblAOP.SetLabel("AZ (degrees)")
            self.lblMA.SetLabel("FPA (degrees)")
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_state_representation == 3:
            self.lblSMA.SetLabel("SMA (km)")
            self.lblECC.SetLabel("ECC")
            self.lblINC.SetLabel("INC (degrees)")
            self.lblRAAN.SetLabel("RAAN (degrees)")
            self.lblAOP.SetLabel("AOP (degrees)")
            self.lblMA.SetLabel("TA (degrees)")
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].arrival_elements_state_representation == 4:
            self.lblSMA.SetLabel("P (km)")
            self.lblECC.SetLabel("F")
            self.lblINC.SetLabel("G")
            self.lblRAAN.SetLabel("H")
            self.lblAOP.SetLabel("K")
            self.lblMA.SetLabel("L (degrees)")
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_state_representation in [5, 6]:
            self.lblSMA.SetLabel("vinf (km/s)")
            self.lblECC.SetLabel("RHA (degrees)")
            self.lblINC.SetLabel("DHA (degrees)")
            self.lblRAAN.SetLabel("bradius (km)")
            self.lblAOP.SetLabel("btheta (degrees)")
            self.lblMA.SetLabel("TA (degrees)")
        elif self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_state_representation in [7, 8]:
            self.lblSMA.SetLabel("vinf (km/s)")
            self.lblECC.SetLabel("RHA (degrees)")
            self.lblINC.SetLabel("DHA (degrees)")
            self.lblRAAN.SetLabel("Rp (km)")
            self.lblAOP.SetLabel("btheta (degrees)")
            self.lblMA.SetLabel("TA (degrees)")

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[0] == 0:
            self.txtSMA.Enable()
            self.txtSMA0.Disable()
            self.txtSMA1.Disable()
        else:
            self.txtSMA.Disable()
            self.txtSMA0.Enable()
            self.txtSMA1.Enable()

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[1] == 0:
            self.txtECC.Enable()
            self.txtECC0.Disable()
            self.txtECC1.Disable()
        else:
            self.txtECC.Disable()
            self.txtECC0.Enable()
            self.txtECC1.Enable()

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[2] == 0:
            self.txtINC.Enable()
            self.txtINC0.Disable()
            self.txtINC1.Disable()
        else:
            self.txtINC.Disable()
            self.txtINC0.Enable()
            self.txtINC1.Enable()

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[3] == 0:
            self.txtRAAN.Enable()
            self.txtRAAN0.Disable()
            self.txtRAAN1.Disable()
        else:
            self.txtRAAN.Disable()
            self.txtRAAN0.Enable()
            self.txtRAAN1.Enable()
            
        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[4] == 0:
            self.txtAOP.Enable()
            self.txtAOP0.Disable()
            self.txtAOP1.Disable()
        else:
            self.txtAOP.Disable()
            self.txtAOP0.Enable()
            self.txtAOP1.Enable()

        if self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[5] == 0:
            self.txtMA.Enable()
            self.txtMA0.Disable()
            self.txtMA1.Disable()
        else:
            self.txtMA.Disable()
            self.txtMA0.Enable()
            self.txtMA1.Enable()
        
        
        self.Layout()
        if platform.system() == 'Windows':
            self.parent.SetupScrolling(scrollToTop=False)

    #event handlers for journey options
    def Changedeparture_elements_representation(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_state_representation = self.cmbelements_representation.GetSelection()
        self.update()
        self.missionoptions.DisassembleMasterDecisionVector()
        self.missionoptions.ConvertDecisionVector()
        self.missionoptions.AssembleMasterDecisionVector()
        e.Skip()
    
    def Changedeparture_elements_frame(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_frame = self.cmbelements_frame.GetSelection()
        self.parent.update()

    def Changedeparture_elements_reference_epoch(self, e):
        e.Skip()

        dateString = self.txtelements_reference_epoch.GetValue()

        from timeUtilities import stringToJD

        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_reference_epoch = stringToJD(dateString, self.missionoptions.universe_folder)

        self.update()

    def ChangeAllowJourneyFreePointDepartureToPropagate(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].AllowJourneyFreePointDepartureToPropagate = int(self.chkAllowJourneyFreePointToPropagate.GetValue())
        self.update()

    def ChangevarySMA_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[0] = int(self.chkSMA.GetValue())
        self.update()

    def ChangevaryECC_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[1] = int(self.chkECC.GetValue())
        self.update()

    def ChangevaryINC_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[2] = int(self.chkINC.GetValue())
        self.update()

    def ChangevaryRAAN_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[3] = int(self.chkRAAN.GetValue())
        self.update()

    def ChangevaryAOP_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[4] = int(self.chkAOP.GetValue())
        self.update()

    def ChangevaryMA_departure(self, e):
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_vary_flag[5] = int(self.chkMA.GetValue())
        self.update()

    def ChangeSMA_departure(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[0] = eval(self.txtSMA.GetValue())
        self.update()

    def ChangeECC_departure(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[1] = eval(self.txtECC.GetValue())
        self.update()

    def ChangeINC_departure(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[2] = eval(self.txtINC.GetValue())
        self.update()

    def ChangeRAAN_departure(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[3] = eval(self.txtRAAN.GetValue())
        self.update()

    def ChangeAOP_departure(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[4] = eval(self.txtAOP.GetValue())
        self.update()

    def ChangeMA_departure(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements[5] = eval(self.txtMA.GetValue())
        self.update()

    def ChangeSMA0(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[0] = eval(self.txtSMA0.GetValue())
        self.update()

    def ChangeECC0(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[2] = eval(self.txtECC0.GetValue())
        self.update()

    def ChangeINC0(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[4] = eval(self.txtINC0.GetValue())
        self.update()

    def ChangeRAAN0(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[6] = eval(self.txtRAAN0.GetValue())
        self.update()

    def ChangeAOP0(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[8] = eval(self.txtAOP0.GetValue())
        self.update()

    def ChangeMA0(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[10] = eval(self.txtMA0.GetValue())
        self.update()

    def ChangeSMA1(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[1] = eval(self.txtSMA1.GetValue())
        self.update()

    def ChangeECC1(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[3] = eval(self.txtECC1.GetValue())
        self.update()

    def ChangeINC1(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[5] = eval(self.txtINC1.GetValue())
        self.update()

    def ChangeRAAN1(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[7] = eval(self.txtRAAN1.GetValue())
        self.update()

    def ChangeAOP1(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[9] = eval(self.txtAOP1.GetValue())
        self.update()

    def ChangeMA1(self, e):
        e.Skip()
        self.missionoptions.Journeys[self.missionoptions.ActiveJourney].departure_elements_bounds[11] = eval(self.txtMA1.GetValue())
        self.update()