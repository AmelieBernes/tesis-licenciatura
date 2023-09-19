import PySimpleGUI as sg
import os.path



#FileBrowse() and FolderBrowse() are different widgets.


layout = [
    [sg.Text("Análisis morfológico y espectral usando los PDL y espacios monofrecuenciales")],
    [
        sg.Text("Archivo .txt con las mediciones:"),
        sg.In(enable_events=True, key="-MEDICIONES-"),
        sg.FileBrowse(),
    ],
    [sg.Text("Tolerancia para el cálculo del grado:"), sg.Input(enable_events = True, key = "-EPSILON-")],
    [sg.Button("Salir", key = "-SALIR-")]
]

window = sg.Window("Análisis con PDL", layout)

while True:
    event, values = window.read() 
    #The 'event' will be the key string of whichever element the user interacts with
    #The 'values' variable contains a Python dictionary that maps the element key to a value.

    #The conditional statements are used to control what happens.
    if event == "-SALIR-" or event == sg.WIN_CLOSED:
        break
    if event == "-MEDICIONES-":
        archivo = values["-MEDICIONES-"] #Tipo string, es la ruta al archivo seleccionado por el usuario.
        if archivo.lower().endswith((".txt")) == False: #si el archivo seleccionado no es un txt,
                                                        #lanzamos un error
            sg.popup("INPUT INVÁLIDO: Debe de seleccionar un archivo de extensión .txt")
            window["-MEDICIONES-"].update("")


window.close()