import PySimpleGUI as sg

from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib
import ae_para_GUI as ae_GUI



matplotlib.use('TkAgg')

#sg.theme('DarkPurple2')
sg.theme('DarkPurple1')

# ------------------------------- Beginning of Matplotlib helper code -----------------------



def draw_figure(canvas, figure, clear = False):
    if clear == True:
        if canvas.children:
            for child in canvas.winfo_children():
                child.destroy()

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
        sg.Button("OK", key = "-OK-"), sg.Button("Cancelar", key = "-CANCEL-")
        ],
        ]

layout = [
    [
        #sg.Column(file_header),
        #sg.VSeparator(),
        file_header,
        inputs, 
        [sg.Text("TDF: ", key = "-TDF-")], 
        [sg.Text("Esp. Monof: ", key = "-EM-")], 
        #sg.Canvas(size=(700,761), background_color='grey', key='-CANVAS-')
        sg.Canvas(size=(660,717), background_color='grey', key='-CANVAS-')
        ]
        ]



# ------------------------ Elementos para la gráfica --------------------------------------


t = np.arange(0, 3, .01)


#Creamos la imagen en la que se mostrarán los espectros pedidos por el usuario
#antes del loop
plt.figure(1, figsize = (6.875, 7.46))
fig = plt.gcf() 
#fig.tight_layout() #does not work.
gs = fig.add_gridspec(2,2)
ax0 = fig.add_subplot(gs[1, :1])
ax1 = fig.add_subplot(gs[1, 1:])
ax2 = fig.add_subplot(gs[0, :2])


# ------------------------ Lógica de la GUI --------------------------------------

window = sg.Window("PDL", layout)

while True:
    event, values = window.read()

    if event == sg.WIN_CLOSED or event == "-CANCEL-":
        break
    if event == "-OK-":
        try:
            #hay que convertir a int, pues el tipo original es str.
            n, k = int(values["-DIM-"]), int(values["-GRADO-"])

            inputs_validez = are_inputs_valid(n, k)

            if inputs_validez  == False:
                sg.popup("INPUT INVÁLIDO")
            else:
                
                label_TDF, label_EM = ae_GUI.grafica_analisisEspectrales_PDL_GUI(fig, n,k)
                draw_figure(window['-CANVAS-'].TKCanvas, fig, clear = True)
                window['-TDF-'].update("TDF: " + label_TDF)
                window['-EM-'].update("Esp. Monof: " + label_EM) #por qué no funciona?? sólo funciona en los casos extremos


        except ValueError:
            sg.popup("INPUT INVÁLIDO")

    #TODO poner la grafica del espectro en un sg.popup!

window.close()

