import PySimpleGUI as sg

from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
matplotlib.use('TkAgg')


sg.theme('DarkPurple2')


#fig = matplotlib.figure.Figure(figsize=(5, 4), dpi=100)
# ------------------------------- Beginning of Matplotlib helper code -----------------------


#Quisiera quitar esta y quedarme con las demás.
def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg



def pack_figure(graph, figure):
    canvas = FigureCanvasTkAgg(figure, graph.Widget)
    plot_widget = canvas.get_tk_widget()
    plot_widget.pack(side='top', fill='both', expand=1)
    return plot_widget

def plot_figure(index, theta):
    fig = plt.figure(index)         # Active an existing figure
    ax = plt.gca()                  # Get the current axes
    x = [degree for degree in range(1080)]
    y = [math.sin((degree+theta)/180*math.pi) for degree in range(1080)]
    ax.cla()                        # Clear the current axes
    ax.set_title(f"Sensor Data {index}")
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.set_xscale('log')
    ax.grid()
    plt.plot(x, y)                  # Plot y versus x as lines and/or markers
    fig.canvas.draw() 


# ----------------------- Funciones auxiliares ---------------------------------

def are_inputs_valid(n, k):
    """
    Función que determina si los input n y k (dimension y grado) dados por el usuario son
    válidos.
    Se checa que ambos sean enteros, que n esté entre 2 y 69 (inclusivo) y que k esté
    entre 0 y n-1. 

    Si los input son válidos, la función regresa True. En caso contrario, regresa False.
    """
    #TODO tal vez debería quitar este int, pues estoy convirtiendo los 
    #valores ingresados 
    if (isinstance(n, int) and isinstance(k, int)) == False:
        return False
    if (2 <= n and n <= 69) == False:
        return False
    if (0 <= k and k <= n-1) == False:
        return False
    return True

# ------------------------ Template de la GUI -----------------------------------

file_header = [
    [
        sg.Text("Análisis espectral de los Polinomios discretos de Legendre", justification = "center")
        ],
    [
        sg.Text("Introduzca la dimensión (entero entre 2 y 69, inclusivo) y el grado (entero no negativo y no mayor a la dimensión) del PDL a analizar")
        ],
        [sg.Text(size=(40, 1), key="-WARNING-")],
        ]

inputs =[
    [
        sg.Text("Dimensión"), sg.InputText(key = "-DIM-")
        ],
    [
        sg.Text("Grado"), sg.InputText(key = "-GRADO-")
        ],
    [
        sg.Button("OK", key = "-OK-")
        ],
    [
        sg.Button("Cancelar", key = "-CANCEL-")
        ]
        ]

layout = [
    [
        sg.Column(file_header),
        #sg.VSeparator(),
        sg.Column(inputs),
        sg.Canvas(size=(700,500), background_color='grey', key='-CANVAS-')
        ]
        ]

window = sg.Window("PDL", layout )


# ------------------------ Lógica de la GUI --------------------------------------

t = np.arange(0, 3, .01)

while True:
    event, values = window.read()

    if event == sg.WIN_CLOSED or event == "-CANCEL-":
        break
    if event == "-OK-":
        try:
            #hay que convertir a int, pues el tipo original es str.
            n, k = int(values["-DIM-"]), int(values["-GRADO-"])

            inputs_validez = are_inputs_valid(n, k)
            #sg.popup(str(inputs_validez) + str(type(n)) + str(k))
            #TODO agregar un botón de clear.
            if inputs_validez  == False:
                #TODO tal vez sea mejor poner un sg.popup!
                window["-WARNING-"].update("INPUT INVÁLIDO", background_color = "black", text_color = "white")
            else:
                window["-WARNING-"].update(" ")

                fig= matplotlib.figure.Figure(dpi=100)
                
                fig.add_subplot(111).plot(t, n*t + k)

                #TODO i need to clear this window first...
                fig_canvas_agg = draw_figure(window['-CANVAS-'].TKCanvas, fig)


        except ValueError:
            window["-WARNING-"].update("INPUT INVÁLIDO", background_color = "black", text_color = "white")

    #TODO poner la grafica del espectro en un sg.popup!

window.close()

