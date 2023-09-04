import PySimpleGUI as sg
import os.path


layout = [
    [sg.Text("Análisis morfológico y espectral usando los PDL y espacios monofrecuenciales")],
    [
        sg.Text("Folder con el archivo .txt con las mediciones:"),
        sg.In(enable_events=True, key="-FOLDER-"),
        sg.FolderBrowse(),
    ],
    [sg.Text("Tolerancia para el cálculo del grado:"), sg.Input(enable_events = True, key = "-EPSILON-")],
    [
        sg.Listbox(
            values=[], enable_events=True, size=(40, 20), key="-FILE LIST-"
        )
    ],
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
    if event == "-FOLDER-":
        folder = values["-FOLDER-"]
        try:
            # Get list of files in folder
            file_list = os.listdir(folder)
        except:
            file_list = []

        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(folder, f))
            and f.lower().endswith((".txt"))
        ]
        window["-FILE LIST-"].update(fnames)

window.close()