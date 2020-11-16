# MAT231 Project B
# TUT0103 Group 9
# Modelling the CoViD-19 Pandemic in Canada
# This basic program asks the user to input the date, population size, number of recovered cases, number of active cases, and the CoViD-19 basic reproductive number, and whether there are vaccines, all on that date, and then creates a model of the pandemic with a graph that has separate curves for individuals who are susceptible, infected, and recovered
# Please note, the code below is based on "The SIR Epidemic Model" code in "Learning Scientific Programming with Python" (cited in the next line in IEEE)
# C. Hill, “The SIR Epidemic Model,” Learning Scientific Programming with Python, 2016. [Online]. Available: https://scipython.com/book/chapter-8-scipy/additional-examples/the-sir-epidemic-model/. [Accessed: 10-Nov-2020]. 

import numpy
from scipy.integrate import odeint
import matplotlib.pyplot as plt


## GETTING USER INPUT

# print the title of the program
print("\nInteractively Modelling the Spread of CoViD-19 in Canada\n")

# Asks user for date
date = input("Enter today's date (DDMMYYYY): ")
new_date = date[0:2] + "/" + date[2:4] + "/" + date[4: ]

# Asks user what region they're inputting data for
region = input("Enter the name of the region for which you are modelling the spread of CoViD-19: ")

# Asks user over how many days they want to model the pandemic
days = int(float(input("Enter how many days you want the spread of CoViD-19 to be modelled over: ")))

# Asks user for the size of the population, N
N = int(float(input("Enter the population size: ")))

# Asks user for iniitial number of active cases, I0
I0 = int(float(input("Enter the number of active CoViD-19 cases today: ")))

# Asks user for initial number of recovered cases, R0
Recovered0 = int(float(input("Enter the current number of recovered CoViD-19 cases: ")))

# Asks user for initial number of deaths
Deaths0 = int(float(input("Enter the current number of deaths due to CoViD-19: ")))

# Asks user for the basic reproduction number of covid
B = float(input("Enter the CoViD-19 basic reproduction number: "))

# Asks user if there is a vaccine, and if so how many are administered each day
A = input("Is there currently a vaccine? (y/n): ")
if A == "Y" or A == "y" :
    vaccine_rate = int(input("Enter number of vaccines administered daily: "))
if A == "N" or A == "n" :
    vaccine_rate = 0


## DEFINING VARIABLES AND EQUATIONS

# Initial number of removed individuals is the sum of the recovered cases and deaths from CoViD-19, as removed individuals have either recovered or died from the disease
R0 = Recovered0 + Deaths0

# Everyone else in the population, S0, is susceptible to the infection initially.
S0 = N - I0 - R0

# Contact rate (beta), and mean recovery rate (gamma), (in 1/days).
beta, gamma = B/14, 1./14

# A grid of time points (in days)
t = numpy.linspace(0, days, days)

# The SIR model differential equations
def deriv(y, t, N, beta, gamma):
    S, I, R = y
    S_prime = ( -beta * S * I / N ) - vaccine_rate
    I_prime = ( beta * S * I / N ) - ( gamma * I )
    R_prime = ( gamma * I ) + vaccine_rate
    return S_prime, I_prime, R_prime

# Initial conditions
y0 = S0, I0, R0
# Integrating the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, gamma))
S, I, R = ret.T


## PLOTTING DATA

# Get the scale factor
if N < 101 :
    scale_factor = 1
    str_scale = ""
elif N < 1001 :
    scale_factor = 100
    str_scale = " (100s)"
elif N < 10001 :
    scale_factor = 1000
    str_scale = " (1,000s)"
elif N < 100001 :
    scale_factor = 10000
    str_scale = " (10,000s)"
elif N < 1000001 :
    scale_factor = 100000
    str_scale = " (100,000s)"
elif N < 100000001 :
    scale_factor = 1000000
    str_scale = " (1,000,000s)"
elif N < 10000000001 :
    scale_factor = 100000000
    str_scale = " (100,000,000s)"   

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor="w")
ax = fig.add_subplot(111, axisbelow=True)
# plot S(t) curve
ax.plot(t, S/scale_factor, 'b', alpha=0.5, lw=2, label='Susceptible Individuals')
# plot I(t) curve
ax.plot(t, I/scale_factor, 'r', alpha=0.5, lw=2, label='Infected Individuals')
# plot R(t) curve
ax.plot(t, R/scale_factor, 'g', alpha=0.5, lw=2, label='Removed Individuals')
# add graph title
fig.suptitle("SIR Model of Coronavirus Pandemic in " + region)
# add axis titles
ax.set_xlabel('Days since ' + new_date)
ax.set_ylabel('Number of People' + str_scale)
# set bounds for x and y axes
ax.set_ylim(0, N / scale_factor )
ax.set_xlim(0,days)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(True)
# display graph
plt.show()
