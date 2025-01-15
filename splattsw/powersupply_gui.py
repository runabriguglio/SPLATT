from splattsw.devices import powersupplier as ps
from guietta import Gui, _, ___, III

class PowerSupplyGui:

    def __init__(self):
        self.power_supply = ps
        self.ch = [False, False, False]
        self.isRunning = False

    def startGui(self):
        gui = Gui(
            [ ['CH1'] , ['CH2'] , ['CH3']],
            [   III   ,   III   ,   III  ],
            [['Start_All'] , ___ , ['EMERGENCY_STOP'] ],
        )

        pass

        def CH1(gui, *args):
            if self.ch[0]:
                self.ch[0] = False
                self.power_supply.switch_off(1)
            else:
                self.ch[0] = True
                self.power_supply.switch_on(1)

        def CH2(gui, *args):
            if self.ch[1]:
                self.ch[1] = False
                self.power_supply.switch_off(2)
            else:
                self.ch[1] = True
                self.power_supply.switch_on(2)

        def CH3(gui, *args):
            if self.ch[2]:
                self.ch[2] = False
                self.power_supply.switch_off(3)
            else:
                self.ch[2] = True
                self.power_supply.switch_on(3)

        def Start_All(gui, *args):
            for n,ch in enumerate(self.ch):
                if ch is False:
                    self.power_supply.switch_on(n+1)
                    self.ch[n] = True
            gui.widgets['Start_All'].setEnabled(False)

        def EMERGENCY_STOP(gui, *args):
            self.power_supply.switch_off(3)
            self.power_supply.switch_off(2)
            self.power_supply.switch_off(1)
            self.ch = [False, False, False]
            gui.widgets['Start_All'].setEnabled(True)

        gui.title('SPLATT Power Supplier')
        gui.window()
        gui.run()

if __name__ == '__main__':
    psg = PowerSupplyGui()
    psg.startGui()
