import os.path
from splattsw.devices import powersupplier as ps # type: ignore
from guietta import Gui, _, ___, III, L

class PowerSupplyGui:

    def __init__(self):
        self.power_supply = ps
        self.ch = [False, False, False]
        self.isRunning = False

    def startGui(self):
        gui = Gui(
            [ ['CH1'] , ['CH2'] , ['CH3'] ],
            [ L("s1") , L("s2") , L("s3") ],
            [ ['Start_All'],_,['EMERGENCY\n  STOP']],
            images_dir = os.path.dirname(__file__)
        )
        off = 'switch_off.png'
        on = 'switch_on.png'

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

        def _check_state():
            for n,_ in enumerate(self.ch):
                state = self.power_supply.switch_state(n+1)
                self.ch[n] = True if state == 'ON' else False

        def _update_graphic():
            """update the graphic for the channel status."""
            gui.s1 = on if self.ch[0] else off
            gui.s2 = on if self.ch[1] else off
            gui.s3 = on if self.ch[2] else off

        gui.events(
            [CH1, CH2, CH3],
            [_  , _  , _  ],
            [Start_All,_, EMERGENCY_STOP]
        )
        _check_state()
        _update_graphic()
        gui.title('SPLATT Power Supplier')
        gui.window()
        gui.run()

if __name__ == '__main__':
    psg = PowerSupplyGui()
    psg.startGui()
