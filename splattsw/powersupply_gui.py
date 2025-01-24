import time
import threading
import os.path as op
from splattsw.devices.powersupplier import PowerSupplier 
from guietta import Gui, _, III, L, G, HB # type: ignore
from PyQt5.QtCore import Qt


class PowerSupplyGui:

    def __init__(self):
        self.power_supply = PowerSupplier()
        self.ch = [False, False, False]
        self.isRunning = False
        self.monitoring_thread = False

    def startGui(self):
        """
        Initializes and starts the graphical user interface (GUI) for controlling the power supply.
        The GUI consists of three main channels (CH1, CH2, CH3) and three control panels for each channel.
        It provides functionalities to start all channels, connect to the power supply, and perform an emergency stop.
        The GUI includes the following components:
        - Channel indicators (CH1, CH2, CH3)
        - Control panels for each channel (control1_gui, control2_gui, control3_gui)
        - Buttons for starting all channels, connecting, and emergency stop
        - LED indicators for channel status (on/off)
        - Real-time monitoring of voltage and current for each channel
        Methods:
        - CH1(gui, *args): Toggles the power state of channel 1.
        - CH2(gui, *args): Toggles the power state of channel 2.
        - CH3(gui, *args): Toggles the power state of channel 3.
        - Start_All(gui, *args): Turns on all channels.
        - EMERGENCY_STOP(gui, *args): Turns off all channels and stops the monitoring thread.
        - connect(gui, *args): Connects to the power supply and starts monitoring.
        - _check_state(): Checks the current state of each channel.
        - _update_graphic(): Updates the GUI to reflect the current state of each channel.
        - _monitoring(gui, *args): Monitors the voltage and current of each channel.
        - live_monitornig(gui, *args): Continuously monitors the power supply status in a separate thread.
        The GUI is titled 'Rigol DP832A' and is resized to 340x300 pixels.
        """
        LBS = ["SIGGen", "BCU+Logic", "Drivers"]
        # Led indicators. Max 50x30 (for now)
        images_dir = op.dirname(__file__)+'/_images/'
        off = [images_dir+f"{n}off.png" for n in [1,2,3]]
        on = [images_dir+f"{n}on.png" for n in [1,2,3]]
        gui = Gui(
            [  ['CH1'] ,  ['CH2']  ,  ['CH3']  ],
            [  L("s1") ,  L("s2")  ,  L("s3")  ],
            [ G('c1'), G('c2') , G('c3') ],
            [    III   ,    III    ,    III    ],
            [   ['SA'] ,   ['CN']  ,  ['ES']   ]
        )
        control1_gui = Gui(
            [ 'V' ],
            [ 'A' ]
        )
        control2_gui = Gui(
            [ 'V' ],
            [ 'A' ]
        )
        control3_gui = Gui(
            [ 'V' ],
            [ 'A' ]
        )
        self.gui = gui
        cguis = [control1_gui, control2_gui, control3_gui]

        def _formatGraphics(gui, *args):
            for ch in range (3):
                cguis[ch].title(LBS[ch])
                gui.widgets['c'+str(ch+1)].setTitle(LBS[ch])
                gui.widgets['c'+str(ch+1)].setAlignment(Qt.AlignHCenter)
                gui.widgets['c'+str(ch+1)].setStyleSheet(
                    "QGroupBox { font-size: 9pt; font-weight: bold; }"
                    )
                gui.widgets['c'+str(ch+1)].setFixedSize(100, 140)
                cguis[ch].V = 'OFF'
                cguis[ch].A = 'OFF'
                cguis[ch].widgets['V'].setContentsMargins(0,0,5,0)
                cguis[ch].widgets['A'].setContentsMargins(0,0,5,0)
                cguis[ch].widgets['V'].setStyleSheet(
                    "QLabel { font-size: 10pt; font-weight: semibold; }"
                    )
                cguis[ch].widgets['A'].setStyleSheet(
                    "QLabel { font-size: 10pt; font-weight: semibold; }"
                    )
                gui.widgets[f"s{ch+1}"].setContentsMargins(25,0,20,0)
            gui.c1 = control1_gui
            gui.c2 = control2_gui
            gui.c3 = control3_gui

###################################################################
 ### Functions for testing graphic updates and event handling ####
###################################################################

        def CH1(gui, *args):
            if self.ch[0]:
                self.ch[0] = False
            else:
                self.ch[0] = True
            _update_graphic()

        def CH2(gui, *args):
            if self.ch[1]:
                self.ch[1] = False
            else:
                self.ch[1] = True
            _update_graphic()


        def CH3(gui, *args):
            if self.ch[2]:
                self.ch[2] = False
            else:
                self.ch[2] = True
            _update_graphic()

        def StartAll(gui, *args):
            for n,ch in enumerate(self.ch):
                if ch is False:
                    self.ch[n] = True
            _update_graphic()

        def EmergencySTOP(gui, *args):
            self.ch = [False, False, False]
            _update_graphic()
            self.isRunning = False

        def connect(gui, *args):
            if not self.isRunning:
                self.isRunning = True
                gui.CN.setText('STOP\nMONITORING')
                self.monitoring_thread = threading.Thread(target=_live_monitoring, args=(gui,))
                self.monitoring_thread.start()
                # change to HeartBeat class
                # gui.s1 = HB(on[0], off[0]).create(gui)
                # gui.s2 = HB(on[1], off[1]).create(gui)
                # gui.s3 = HB(on[2], off[2]).create(gui)
            else:
                self.isRunning = False
                gui.CN.setText('START\nMONITORING')
                self.monitoring_thread.join()
                for ch in range(3):
                    voltage = 'null'
                    current = 'null'
                    cguis[ch].V = voltage
                    cguis[ch].A = current

        def _live_monitoring(gui, *args):
            import numpy as np
            while self.isRunning:
                for ch in range(3):
                    if self.ch[ch]:
                        voltage = f"{np.random.uniform(0, 5.5):.3f} V"
                        current = f"{np.random.uniform(0, 2.5):.3f} A"
                        cguis[ch].V = voltage
                        cguis[ch].A = current
                    else:
                        voltage = 'null'
                        current = 'null'
                        cguis[ch].V = voltage
                        cguis[ch].A = current
                time.sleep(0.6)


###################################################################
  ####                        END                           #####
###################################################################

        # def CH1(gui, *args):
        #     if self.ch[0]:
        #         self.power_supply.switch_off(1)
        #         self.ch[0] = False
        #     else:
        #         self.power_supply.switch_on(1)
        #         self.ch[0] = True
        #     _update_graphic()

        # def CH2(gui, *args):
        #     if self.ch[1]:
        #         self.power_supply.switch_off(2)
        #         self.ch[1] = False
        #     else:
        #         self.power_supply.switch_on(2)
        #         self.ch[1] = True
        #     _update_graphic()

        # def CH3(gui, *args):
        #     if self.ch[2]:
        #         self.power_supply.switch_off(3)
        #         self.ch[2] = False
        #     else:
        #         self.power_supply.switch_on(3)
        #         self.ch[2] = True
        #     _update_graphic()

        # def StartAll(gui, *args):
        #     for n,ch in enumerate(self.ch):
        #         if ch is False:
        #             self.power_supply.switch_on(n+1)
        #             self.ch[n] = True
        #     _update_graphic()

        # def EmergencySTOP(gui, *args):
        #     self.power_supply.switch_off(3)
        #     self.power_supply.switch_off(2)
        #     self.power_supply.switch_off(1)
        #     self.ch = [False, False, False]
        #     _update_graphic()
        #     self.isRunning = False

        # def connect(gui, *args):
        #     if not self.isRunning:
        #         self.isRunning = True
        #         gui.CN.setText('STOP\nMONITORING')
        #         self.monitoring_thread = threading.Thread(target=live_monitoring, args=(gui,))
        #         self.monitoring_thread.start()
        #     else:
        #         self.isRunning = False
        #         gui.CN.setText('START\nMONITORING')
        #         if self.monitoring_thread:
        #             self.monitoring_thread.join()
                
        # def live_monitoring(gui, *args):
        #     while self.isRunning:
        #         _check_state()
        #         _monitoring(gui)
        #         time.sleep(0.6)

        def _check_state():
            for n,_ in enumerate(self.ch):
                state = self.power_supply.switch_state(n+1)
                self.ch[n] = True if state == 'ON' else False

        def _update_graphic():
            """update the graphic for the channel status."""
            gui.s1 = on[0] if self.ch[0] else off[0]
            gui.s2 = on[1] if self.ch[1] else off[1]
            gui.s3 = on[2] if self.ch[2] else off[2]

        def _monitoring(gui, *args):
            """Monitor the power supply status."""
            for ch in range(1,4):
                if self.ch[ch-1]:
                    voltage = f"{self.power_supply.read_voltage(ch):.3f} V"
                    current = f"{self.power_supply.read_current(ch):.3f} A"
                else:
                    voltage = 'OFF'
                    current = 'OFF'
                cguis[ch-1].V = voltage
                cguis[ch-1].A = current

        def _closeEvent(event):
            self.isRunning = False
            if self.monitoring_thread:
                self.monitoring_thread.join()
            event.accept()

        
        gui.events(
            [CH1, CH2, CH3],
            [_  , _  , _  ],
            [_  , _  , _  ],
            [_  , _  , _  ],
            [StartAll,connect,EmergencySTOP]
        )

        gui.title('Rigol DP832A')
        gui.SA.setText('START\n ALL')
        gui.ES.setText('EMERGENCY\n  STOP')
        gui.CN.setText('START\nMONITORING')
        _formatGraphics(gui)
        #_check_state()
        _update_graphic()
        gui.window()
        gui.window().resize(340, 300)
        gui.window().closeEvent = _closeEvent
        gui.run()

if __name__ == '__main__':
    psg = PowerSupplyGui()
    psg.startGui()
