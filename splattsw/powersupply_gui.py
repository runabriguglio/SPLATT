import os.path
import threading
from splattsw.devices import powersupplier as ps # type: ignore
from guietta import Gui, _, ___, III, L, G

class PowerSupplyGui:

    def __init__(self):
        self.power_supply = ps
        self.ch = [False, False, False]
        self.isRunning = False

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
        LBS = ["         SIGGen", "        BCU+Logic", "        Driver+EMStop"]
        gui = Gui(
            [  ['CH1'] ,  ['CH2'] ,  ['CH3']  ],
            [  L("s1") ,  L("s2") ,  L("s3")  ],
            [  G(LBS[0]) ,  G(LBS[1]) ,  G(LBS[2])  ],
            [    III   ,    III   ,    III    ],
            [['START\n  ALL'],['CONNECT'],['EMERGENCY\n  STOP']],
            images_dir = os.path.dirname(__file__)
        )
        control1_gui = Gui(
            [ 'V1' ],
            [ 'A1' ]
        )
        control2_gui = Gui(
            [ 'V2' ],
            [ 'A2' ]
        )
        control3_gui = Gui(
            [ 'V3' ],
            [ 'A3' ]
        )
        control1_gui.title(LBS[0])
        control2_gui.title(LBS[1])
        control3_gui.title(LBS[2])
        # Led indicators. Max 50x30 (for now)
        off = 'coff.png'
        on = 'con.png'
        gui.C1 = control1_gui
        gui.C2 = control2_gui
        gui.C3 = control3_gui

        def CH1(gui, *args):
            if self.ch[0]:
                self.power_supply.switch_off(1)
                self.ch[0] = False
            else:
                self.power_supply.switch_on(1)
                self.ch[0] = True
            _update_graphic()

        def CH2(gui, *args):
            if self.ch[1]:
                self.power_supply.switch_off(2)
                self.ch[1] = False
            else:
                self.power_supply.switch_on(2)
                self.ch[1] = True
            _update_graphic()

        def CH3(gui, *args):
            if self.ch[2]:
                self.power_supply.switch_off(3)
                self.ch[2] = False
            else:
                self.power_supply.switch_on(3)
                self.ch[2] = True
            _update_graphic()

        def Start_All(gui, *args):
            for n,ch in enumerate(self.ch):
                if ch is False:
                    self.power_supply.switch_on(n+1)
                    self.ch[n] = True
            gui.widgets['Start_All'].setEnabled(False)
            _update_graphic()

        def EMERGENCY_STOP(gui, *args):
            self.power_supply.switch_off(3)
            self.power_supply.switch_off(2)
            self.power_supply.switch_off(1)
            self.ch = [False, False, False]
            gui.widgets['Start_All'].setEnabled(True)
            _update_graphic()
            self.isRunning = False

        def connect(gui, *args):
            self.isRunning = True
            _check_state()
            _update_graphic()
            live_monitornig(gui)

        def _check_state():
            for n,_ in enumerate(self.ch):
                state = self.power_supply.switch_state(n+1)
                self.ch[n] = True if state == 'ON' else False

        def _update_graphic():
            """update the graphic for the channel status."""
            gui.s1 = on if self.ch[0] else off
            gui.s2 = on if self.ch[1] else off
            gui.s3 = on if self.ch[2] else off

        def _monitoring(gui, *args):
            """Monitor the power supply status."""
            for ch in (3):
                voltage = self.power_supply.read_voltage(ch)
                current = self.power_supply.read_current(ch)
                gui.widgets['V'+str(ch+1)].setText(str(voltage))
                gui.widgets['A'+str(ch+1)].setText(str(current))

        def live_monitornig(gui, *args):
            while self.isRunning:
                self.monitoring_thread = threading.Thread(target=_monitoring, args=(gui,))
                self.monitoring_thread.start()
        
        gui.events(
            [CH1, CH2, CH3],
            [_  , _  , _  ],
            [Start_All,connect, EMERGENCY_STOP]
        )

        gui.title('Rigol DP832A')
        gui.window()
        _update_graphic()
        gui.window().resize(340, 300)
        gui.run()

if __name__ == '__main__':
    psg = PowerSupplyGui()
    psg.startGui()
