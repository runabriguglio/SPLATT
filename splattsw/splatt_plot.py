import numpy as np
import matplotlib.pyplot as plt


def plot_splatt_data(values):
    # These variables, togheter with folders, should be initialized differently...
    nActs = 19
    coordAct = np.array([[0,0],
    [0.0000,-0.0925],
    [-0.0801,-0.0463],
    [-0.0801,0.0462],
    [-0.0000,0.0925],
    [0.0801,0.0463],
    [0.0801,-0.0462],
    [0.0000,-0.1855],
    [-0.0927,-0.1606],
    [-0.1606,-0.0928],
    [-0.1855,-0.0000],
    [-0.1606,0.0927],
    [-0.0927,0.1606],
    [-0.0000,0.1855],
    [0.0927,0.1606],
    [0.1606,0.0928],
    [0.1855,0.0000],
    [0.1606,-0.0927],
    [0.0928,-0.1606]])

    # Perform matrix rotation to align with 'gravity'
    phi = 60./180*np.pi
    c=np.cos(phi)
    s=np.sin(phi)
    MatRot=[[c,-s],[s,c]]
    coordAct = MatRot@coordAct.T
    coordAct = coordAct.T

    # Set scatter plot variables
    Margin = 0.03
    markerSize = 800
    x = coordAct[:,0]
    y = coordAct[:,1]
    indices = np.arange(nActs)+1

    # Plot
    plt.figure()
    ax = plt.axes()
    ax.set_xlim(min(x)-Margin,max(x)+Margin)
    ax.set_ylim(min(y)-Margin,max(y)+Margin)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.scatter(x, y, c=values, s=markerSize, edgecolor='k', cmap='hot')
    plt.colorbar()

    # Write 'G' reference and actuator indices
    for i in range(nActs):
        plt.text(x[i]*2/3, y[i]+Margin*2/3, str(indices[i]))
    plt.text(x[15],y[15]*1.3,'G')
    plt.show()