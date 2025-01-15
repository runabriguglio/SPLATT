import socket
import time

ip = "192.168.0.101"
port = 5555
delimiter = '\r\n'
rigol = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

def open_read_ch1():
    try:
        print("Connetto")
        RDP = rigol.connect((ip,port))
        RDP.settimeout(2)
        print("Accendo il canale 1...")
        RDP.send((":OUTP CH1,ON" + delimiter).encode('utf-8'))
        time.sleep(0.2)
        print("Leggo lo stato del canale 1...")
        RDP.send((":OUTP? CH1" + delimiter).encode('utf-8'))
        time.sleep(0.2)
        response = RDP.recv(4096).decode('utf-8').strip()
        print(f"Stato canale 1: {response}")
        print("Controllo errori...")
        RDP.send((":SYST:ERR?" + delimiter).encode('utf-8'))
        error_response = RDP.recv(4096).decode('utf-8').strip()
        print(f"Errore: {error_response}")
    finally:
        RDP.close()

def read(cs=4096):
    """Receive text string and return it after removing the delimiter."""
    RDP = rigol.connect((ip, port))
    RDP.settimeout(2)
    msg = ''
    while True:
        chunk = RDP.recv(cs).decode('utf-8')
        print(f"Ricevuto chunk: {chunk!r}")  # Log del chunk
        msg += chunk
        if msg.endswith(delimiter):
            break
    print(f"Messaggio completo: {msg!r}")  # Log del messaggio completo
    return msg[: -len(delimiter)]
