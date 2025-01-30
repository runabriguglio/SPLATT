import time
import threading
import numpy as np
from PyQt5.QtCore import Qt
from matplotlib import pyplot as plt
from matplotlib.pyplot import tight_layout
from splattsw.devices import wavegenerators as wg
from guietta import Gui, _, ___, III, B, M, G, HValueSlider as hvs, execute_in_main_thread

class GeneratorGui():

    def __init__(self):
        self.generator = wg.WaveGenerator()
        self.sent = False
        self.waveType=None
        # DEFINIRE ALTRI PARAMETRI UTILI



    def startGui(self):
        """
        Questa sarà la funzione che avvierà la GUI, e in cui saranno
        definite tutte le funzioni che faranno accadere cose
        """

        fslider = hvs('Frequency', unit='   [Hz]', myrange=range(1,120,1), default=1)
        aslider = hvs('Amplitude', unit='   [mV]', myrange=range(0,2000,100), default=1)
        pslider = hvs('Phase'    , unit='  [deg]', myrange=range(0,359,1), default=0)
        gui = Gui(
            [M('plot'),___,  ___  ],
            [ III  , III  ,  III  ],
            [ III  , III  ,  III  ],
            [G('CONTROLS'),___,___],
            [ III  , III  ,  III ,],
            [ III  , III  ,  III ,],
            [ III  , III  ,  III ,],
            [ III  , III  ,  III ,]
        )

        gui_panel = Gui(
            [     fslider   ,    ___    , 'Frequency'  ],
            [     aslider   ,    ___    , 'Amplitude'  ],
            [     pslider   ,    ___    , 'Wave Phase' ],
            [     ['Sine']  , ['Square'],   ['SEND']   ],
            [     ['Pulse'] , ['Sweep'] ,   ['Close']  ]
        )
        gui.CONTROLS = gui_panel
        gui_panel.main_gui = gui
        gui.widgets['CONTROLS'].setTitle("Wave Function Generator Control Panel")
        gui.widgets['CONTROLS'].setStyleSheet("QGroupBox { font-size: 9pt; font-weight: bold; }")
        gui.widgets['CONTROLS'].setAlignment(Qt.AlignHCenter)
        gui_panel.widgets['PLOT'].setStyleSheet("background-color: green; font: bold")

        def sin_wave(gui, *arg):
            """

            """

            # SCRIVERE LA FUNZIONE CHE, PRESI I PARAMETRI GIUSTI, MANDA IL COMANDO
            # DI GENERAZIONE D'ONDA

        def square_wave(gui, *arg):
            """

            """

            # SCRIVERE LA FUNZIONE CHE, PRESI I PARAMETRI GIUSTI, MANDA IL COMANDO
            # DI GENERAZIONE D'ONDA

        def pulse_wave(gui, *arg):
            """

            """
            # SCRIVERE LA FUNZIONE CHE, PRESI I PARAMETRI GIUSTI, MANDA IL COMANDO
            # DI GENERAZIONE D'ONDA

        def sweep_wave(gui, *args):
            """
            
            """

            # SCRIVERE LA FUNZIONE CHE, PRESI I PARAMETRI GIUSTI, MANDA IL COMANDO
            # DI GENERAZIONE D'ONDA


        def send_wave(gui, *args):
            """
            
            """
            if self.sent is False:
                gui.SEND.setStyleSheet("background-color: red; font: bold")
                _start_loop()
            else:
                self.sent = False
                gui.SEND.setStyleSheet("background-color: green; font: bold")
                if self.update_thread:
                    self.update_thread.join()

        def _start_loop():
            if not self.sent:
                self.sent = True
                self.update_thread = threading.Thread(target=_looping_function, args=(gui,))
                self.update_thread.start()

        def _looping_function(gui, *args):
            while self.sent:
                plot_wave()
                time.sleep(0.05)

        @execute_in_main_thread(gui)
        def plot_wave(*arg):
            """

            """
            ax = gui.plot.ax
            ax.clear()
            f,x = self._waveFunction()
            lw = 3 if self.freq < 20 else 1
            ax.plot(x ,f, linewidth=lw)
            ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
            ax.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
            ax.grid(True)
            ax.set_title("Wave Form")
            ax.set_in_layout(tight_layout)
            ax.figure.canvas.draw()

        def close_gui(sub_gui, *arg):
            """

            """
            if self.sent:
                self.sent = False
            if self.update_thread:
                self.update_thread.join()
            sub_gui.main_gui.close()

        def _close_event(event):
            if self.sent:
                self.sent = False
            if self.update_thread:
                self.update_thread.join()
            event.accept()

        def _freq(gui, value):
            self.freq = value

        def _ampl(gui, value):
            self.amp = value/10

        def _phase(gui, value):
            self.phase = (value*np.pi)/180


        gui.events(
            [plot_wave ,     _      ,     _    ],
            [    _     ,     _      ,     _    ],
            [    _     ,     _      ,     _    ]
        )
        
        gui_panel.events(
            [  _freq   ,     _      ,     _    ],
            [  _ampl   ,     _      ,     _    ],
            [ _phase   ,     _      ,     _    ],
            [sin_wave  , square_wave, send_wave],
            [pulse_wave,     _      , close_gui]
        )

        gui.window()
        gui.window().closeEvent = _close_event
        gui.run()
        

    def _waveFunction(self):
        x = np.arange(0,2*np.pi,np.pi/1000)
        if self.waveType == 'sine':
            f = self.amp*np.sin(self.freq*x + self.phase)
        elif self.waveType == 'square':
            f = self.amp * np.sign(np.sin(self.freq * x + self.phase))
        elif self.waveType == 'pulse':
            f = self.amp * (np.sin(self.freq * x + self.phase) > 0).astype(float)
        else:
            f = np.zeros(len(x))
        return f, x

if __name__=='__main__':
    gui = GeneratorGui()
    gui.startGui()
