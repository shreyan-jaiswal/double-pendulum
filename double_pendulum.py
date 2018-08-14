"""
Shreyan Jaiswal
Inspired by http://www.physics.usyd.edu.au/~wheat/dpend_html/
"""

from numpy import square, sin, cos, divide, round, pi
import numpy as np
from tkinter import *
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D
from scipy import integrate

L1, L2, M1, M2, G  = 1.00, 1.00, 2.00, 2.00, 9.80
state = [pi, .01, 2*pi/3, 0]
dt = 1./60

def derivative(state, t):
    global L1, L2, M1, M2, G
    theta1, omega1, theta2, omega2 = state
    delta = theta2 - theta1
    dydt = np.zeros_like(state)

    dydt[0] = omega1 #d/dt(theta1) = omega1
    dydt[2] = omega2 #d/dt(theta2) = omega2

    dydt[1] = divide(M2*L1*square(omega1)*sin(delta)*cos(delta) +
                          M2*G*sin(theta2)*cos(delta) + M2*L2*square(omega2)*sin(delta) -
                          (M2+M1)*G*sin(theta1),
                          (M2+M1)*L1-M2*L1*square(cos(delta))
              )
    dydt[3] = divide(-M2*L2*square(omega2)*sin(delta)*cos(delta) +
                        (M1+M2)*(G*sin(theta1)*cos(delta)-L1*square(omega1)*sin(delta)-G*sin(theta2)),
                        (M1+M2)*L2 - M2*L2*square(cos(delta))
                        )

    return dydt

def step(state):
    global dt
    new_state = integrate.odeint(derivative, state, [0, dt])[1]
    return new_state

def xy(state):
    global L1, L2
    theta1, omega1, theta2, omega2 = state
    return L1*sin(theta1), -L1*cos(theta1), L1*sin(theta1)+L2*sin(theta2), -L1*cos(theta1)-L2*cos(theta2)


class GUIFrame(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)

        self.gui = {}
        for var in ['L1', 'L2', 'M1', 'M2', 'G']: self.gui[var] = {}

        self.gui['L1']['Label'] = Label(self, text="Length of Rod 1")
        self.gui['L1']['Scale'] = Scale(self, from_=0.1, to=2, orient=HORIZONTAL, resolution=0.01, command = lambda value: self.onScaleInput('L1', value))
        self.gui['L1']['Scale'].set(1.0)
        self.gui['L1']['Button'] = Button(self, text='Set From Box', command=lambda: self.onButtonPress('L1'))

        self.gui['L2']['Label'] = Label(self, text="Length of Rod 2")
        self.gui['L2']['Scale'] = Scale(self, from_=0.1, to=2, orient=HORIZONTAL, resolution=0.01,  command = lambda value: self.onScaleInput('L2', value))
        self.gui['L2']['Scale'].set(1.0)
        self.gui['L2']['Button'] = Button(self, text='Set From Box', command=lambda: self.onButtonPress('L2'))

        self.gui['M1']['Label'] = Label(self, text="Mass 1")
        self.gui['M1']['Scale'] = Scale(self, from_=0.01, to=20, orient=HORIZONTAL, resolution=0.01, command = lambda value: self.onScaleInput('M1', value))
        self.gui['M1']['Scale'].set(2.0)
        self.gui['M1']['Button'] = Button(self, text='Set From Box', command=lambda: self.onButtonPress('M1'))

        self.gui['M2']['Label'] = Label(self, text="Mass 2")
        self.gui['M2']['Scale'] = Scale(self, from_=0.01, to=20, orient=HORIZONTAL, resolution=0.01, command = lambda value: self.onScaleInput('M2', value))
        self.gui['M2']['Scale'].set(2.0)
        self.gui['M2']['Button'] = Button(self, text='Set From Box', command=lambda: self.onButtonPress('M2'))

        self.gui['G']['Label'] = Label(self, text="Gravity")
        self.gui['G']['Scale'] = Scale(self, from_=0.01, to=100, orient=HORIZONTAL, resolution=0.01, command = lambda value: self.onScaleInput('G', value))
        self.gui['G']['Scale'].set(9.80)
        self.gui['G']['Button'] = Button(self, text='Set From Box', command=lambda: self.onButtonPress('G'))

        for var in ['L1', 'L2', 'M1', 'M2', 'G']:
            self.gui[var]['StringVar'] = StringVar()
            self.gui[var]['Box'] = Entry(self, width=8, textvariable=self.gui[var]['StringVar'])

        for r, var in enumerate(['L1', 'L2', 'M1', 'M2', 'G']):
            self.gui[var]['Label'].grid(row = r, column = 0)
            self.gui[var]['Scale'].grid(row=r, column=1)
            self.gui[var]['Box'].grid(row=r, column=2)
            self.gui[var]['Button'].grid(row=r, column=3, padx = (0, 10))

    def onScaleInput(self, var, value):
        self.gui[var]['Box'].delete(0, END)
        self.gui[var]['Box'].insert(0, str(value))
        self.setVar(var, float(value))

    def onButtonPress(self, var):
        try:
            value = round(float(self.gui[var]['StringVar'].get()), 2)
        except:
            value = self.gui[var]['Scale'].get()
        self.gui[var]['Box'].delete(0, END)
        self.gui[var]['Box'].insert(0, str(value))
        self.gui[var]['Scale'].set(value)
        self.setVar(var, float(value))

    def setVar(self, var, value):
        global L1, L2, M1, M2, G
        if var == 'L1':L1 = value
        elif var=='L2':L2 = value
        elif var=='M1':M1 = value
        elif var=='M2':M2 = value
        elif var== 'G': G = value


root = Tk()
root.title("Double Pendulum")

gui_frame = GUIFrame(root)
gui_frame.pack(side = RIGHT)

animation_frame = Frame(root)
animation_frame.pack(side = LEFT)

fig = plt.figure()
subplot1 = plt.subplot()
subplot1.set_xlim([-2, 2])
subplot1.set_ylim([-2, 2])
x1_init, y1_init, x2_init, y2_init = state
line, = subplot1.plot([0, x1_init, x2_init], [0, y1_init, y2_init], 'o-', lw = 2, ms = 20)

canvas = FigureCanvasTkAgg(fig, animation_frame)
canvas.get_tk_widget().pack()

def animate(i):
    global state
    state = step(state)
    x1, y1, x2, y2 = xy(state)
    line.set_data([0, x1, x2], [0, y1, y2])
    return line,


from time import time
t0 = time()
animate(0)
t1 = time()
interval = 1000 * dt - (t1 - t0)
ani = animation.FuncAnimation(fig, animate, interval=interval, blit = True)

root.mainloop()



