
import PySimpleGUI as sg

#We save the elements of the GUI on the layout array.
layout = [[sg.Text("Polinomios discretos de Legendre")], [sg.Button("OK")]]

#Creating the window
window = sg.Window("Demo", layout)

#Create an event loop: infinite while loop that reads events from the 'window' object.
while True:
    event, values = window.read()
    #End program if user closes the window or presses the OK button
    if event == "OK" or event == sg.WIN_CLOSED:
        break

window.close()



#Some theory facts (for myself!)
# 1.- PySimpleGUI uses nested Python lists to lay out its elements.
