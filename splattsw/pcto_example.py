# In questo file c'è la bozza per costruire una GUI (Graphic User Interface) 
# per il generatore di funzioni d'onda Rigol DG1000Z. L'obiettico è quello di
# procurarsi i comandi giusti per poter applicare una corretta funzione d'onda
# al generatore, e di poter visualizzare i parametri della funzione d'onda nella GUI.
# NOTA:
#     Non eliminare nulla di ciò che è già scritto, in quanto è necessario al
#     corretto funzionamento finale dello script. Inserire il codice necessario
#     negli spazi indicati con i commenti.
#
#     Se si usa VSCode, con ctrl+click si può andare alla definizione di una funzione
#     o di una variabile, e con ctrl+spazio si possono vedere le possibili funzioni
#     o variabili da utilizzare.

from guietta import Gui, _, B
power = 0
#######################################################################################
# Questo è un esempio di gui molto semplice da cui partire, per capire come collegare  #
# i vari elementi tra loro. C'è una variabile `gui`, in cui si definisce la "grafica"  #
# della GUI; nel mezzo ci sono le varie funzioni che di collegano agli elementi grafi- #
# ci definiti sopra; infine, si collegano i due con il `gui.events`, in cui si mettono #
# i nomi delle funzioni definite in corrispondenza degli elementi grafici. L'ultimo    #
# comando, `gui.run()`, avvia la GUI.                                                  #
#######################################################################################
gui = Gui(
    [      _      ,"Generatore di\nfunzioni d'onda",     _      ],
    [ B('Genera') ,            '__p__'             ,    'P'     ],
    [      _      ,               _                , B('chiudi')]
)
gui.P = power

def genera(gui, *arg):
    """
    Questa funzione sarà chiamata quando si preme il bottone "Genera"
    """
    power = float(gui.p)
    gui.P = f"{power:.2f} W"

def chiudi(gui, *arg):
    """
    Questa funzione sarà chiamata quando si preme il bottone "chiudi"
    """
    gui.close()

gui.events(
    [   _   , _ ,    _   ],
    [genera , _ ,    _   ],
    [   _   , _ , chiudi ]
)

gui.window()
gui.window.title('')
gui.run()
########################################################################################
# Fine esempio. Per visualizzare il risultato di questa GUI di esempio, per prima cosa  #
# rimuovere tutti i commenti (ovvero i `#`) alle linee di codice, e dopo andare su ter- #
# minale, assicurarsi di essere nella stessa cartella in cui si trova questo file, e    #
# mandare `python3 PCTO_GUI.py`.                                                        #
########################################################################################