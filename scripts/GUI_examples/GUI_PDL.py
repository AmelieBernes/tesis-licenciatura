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

def draw_figure(canvas, figure):
    figure_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    figure_canvas_agg.draw()
    figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
    return figure_canvas_agg



# ----------------------- Funciones auxiliares ---------------------------------

def are_inputs_valid(n, k):
    """
    Función que determina si los input n y k (dimension y grado) dados por el usuario son
    válidos.
    Se checa que ambos sean enteros, que n esté entre 2 y 69 (inclusivo) y que k esté
    entre 0 y n-1. 

    Si los input son válidos, la función regresa True. En caso contrario, regresa False.
    """
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
        ]
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
        sg.VSeparator(),
        sg.Column(inputs),
        ]
        ]

window = sg.Window("PDL", layout )

while True:
    event, values = window.read()
    if event == sg.WIN_CLOSED or event == "-CANCEL-":
        break
    if event == "-OK-":
        try:
            inputs_validez = are_inputs_valid("-DIM-", "-GRADO-")
        except:
            pass


window.close()
