"""
Linear Combination of Debye and Einstein Functions for Heat Capacity 
by Jon A Newnham
"""



"""Data Set"""

StartT= 1.8 #K   <- Model Starting Temperature
EndT = 300 #K    <- Model Ending Temperature
increment = 0.1 #K <- Model Temperature Steps

data_location = "example_data.csv"   # <- two column .CSV file with no headers 


""""Refined Perameters"""

#Sum of each the prefactors (Td_p + Te_p) must equal 1.
#Turn off a component by setting the prefactor value to 0.


#Debye Components   

#1st Debye Temperature
Td1 = 300 #K
Td1p = 0.40

#2nd Debye Temperature
Td2 = 840 #K
Td2p = 0.5374

#3rd Debye Temperature
Td3 = 360 #K
Td3p = 0.0



#Einstine Components

#1st Einstein Temperature
Te1 = 69 #K
Te1p = 0.06

#2nd Einstein Temperature
Te2 = 30 #K
Te2p = 0.0026

#3rd Einstein Temperature
Te3 = 12 #K
Te3p = 0



#Linear Component
gamma = 0.0018 #J mol-1 K-2


#Schottkey Anomaly
delta = 5.0 #K
defect_conc = 0.0011 #Set this to zero to turn off Schottkey Anomaly



"""Plotting"""

Data_view = 0

#Data_view =
# 0 for Cp(T)
# 1 for Cp/T^3(T)
# 2 for Cp/T(T^2) full range
# 3 for Cp/T(T^2) low T range
# 4 for all plots togeater


export_graph = 0 # set to 1 to export graph as a .png
export_model = 0 # set to 1 to export your model as a .csv



#Do not edit below here


"""Header"""
from math import e
from scipy.constants import R
import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib.pyplot as plt



"""Einstein Function"""
def Einstein(T, Te):
    return ((3*R*((Te/T)**2))*((e**(Te/T))/((e**(Te/T)-1)**2)))


"""Debye Function"""
def integrate(x):
    return (np.power(x,4)*np.power(e,x))/(np.power((np.power(e,x)-1),2))

def Debye(T, Td):
    return (9*R*((T/Td)**3))*quad(integrate,0,Td/T)[0]


"""Linear Function"""
def y(gamma, T):
    return gamma*T


""""Schottky Anomaly"""
def Schottky(T):
    return defect_conc * R * ((delta/T)**2) *( (e**(delta/T))   / (  (1+(e**(delta/T)))**2  ) )


"""Check sum to 1"""
einstein_comps = [(Te1, Te1p), (Te2, Te2p), (Te3, Te3p)]
debye_comps = [(Td1, Td1p), (Td2, Td2p), (Td3, Td3p)]

if Td1p + Td2p + Td3p + Te1p + Te2p + Te3p != 1:
    print("\nDebye and Einstein Prefactors Do Not Sum to 1!\n")



"""import data"""
data = pd.read_csv(data_location, header=None)


"""plotting"""

xaxis = []

debye_ys = [[], [], []]
einstein_ys = [[], [], []]
linear_ys = []
schottky_ys = []
temps = np.arange(StartT, EndT, increment)

for T in temps:
    for ys, (Td, Tdp) in zip(debye_ys, debye_comps):
        if Tdp:
            ys.append((Tdp * Debye(T, Td)))
    for ys, (Te, Tep) in zip(einstein_ys ,einstein_comps):
        if Tep:
            ys.append((Tep * Einstein(T, Te)))
    if gamma:
        linear_ys.append((y(gamma, T)))  
    if delta:
        schottky_ys.append(Schottky(T))       

all_ys = zip(*[x for x in debye_ys + einstein_ys + [linear_ys] + [schottky_ys] if len(x)>0])
totaly = list(map(sum,all_ys))

data_dict = {"T":temps, "Total":totaly}


#I'm a bit embarrassed by the code in this section, but hey, it works


for i, ys in enumerate(debye_ys):
    if len(ys)>0:
        data_dict[f"D{i+1}"] = ys
        if Data_view == 0:
            plt.plot(temps, ys, label =f"D{i+1}")
        if Data_view == 1:
            plt.plot(temps, ys/temps**3, label =f"D{i+1}")
        if Data_view == 2 or Data_view == 3:
            plt.plot(temps**2, ys/temps, label =f"D{i+1}")
            
for i, ys in enumerate(einstein_ys):
    if len(ys)>0:
        data_dict[f"E{i+1}"] = ys
        if Data_view == 0:
            plt.plot(temps, ys, label =f"E{i+1}")
        if Data_view == 1:
            plt.plot(temps, ys/temps**3, label =f"E{i+1}")
        if Data_view == 2 or Data_view == 3:
            plt.plot(temps**2, ys/temps, label =f"E{i+1}")  
            
if gamma:
    data_dict["Linear"] = linear_ys
    if Data_view == 0:
        plt.plot(temps, linear_ys, label="Linear")  
    if Data_view == 1:
        plt.plot(temps, linear_ys/temps**3, label="Linear") 
    if Data_view == 2 or Data_view == 3:
        plt.plot(temps**2, linear_ys/temps, label="Linear")
        
if delta:
     data_dict["Schottky"] = schottky_ys
     if Data_view == 0:
         plt.plot(temps, schottky_ys, label="Schottky")
     if Data_view == 1:
         plt.plot(temps, schottky_ys/temps**3, label="Schottky")
     if Data_view == 2 or Data_view == 3:
        plt.plot(temps**2, schottky_ys/temps, label="Schottky")

if Data_view == 0:
    plt.scatter(data[0],data[1], c="r", s=1, label="data")
    plt.plot(temps, totaly, c="black", label="total")
    plt.xlim([StartT, EndT])
    plt.ylim([0, None])
    plt.ylabel("C$_p$ (Jmol$^{-1}$K$^{-1}$)")
    plt.xlabel("Temperature (K)")
    plt.legend()

if Data_view == 1:
    plt.scatter(data[0],data[1]/data[0]**3, c="r", s=1, label="data")
    plt.plot(temps, totaly/temps**3, c="black", label="total")
    plt.xlim([StartT, EndT])
    plt.xscale("log")
    plt.ylim([0, None])
    plt.xlabel("Temperature (K)")
    plt.ylabel("C$_p$/T$^3$ (Jmol$^{-1}$K$^{-4}$)")
    plt.legend()

if Data_view == 2:
    plt.scatter(data[0]**2,data[1]/data[0], c="r", s=1, label="data")
    plt.plot(temps**2, totaly/temps, c="black", label="total")
    plt.xlim([StartT**2, EndT**2])
    plt.ylim([0, None])
    plt.xlabel("Temperature$^2$ (K$^2$)")
    plt.ylabel("C$_p$/T (Jmol$^{-1}$K$^{-2}$)")
    plt.legend()

if Data_view == 3:
    plt.scatter(data[0]**2,data[1]/data[0], c="r", s=1, label="data")
    plt.plot(temps**2, totaly/temps, c="black", label="total")
    plt.xlim([StartT**2, 40])
    plt.ylim([0, 0.007])
    plt.xlabel("Temperature$^2$ (K$^2$)")
    plt.ylabel("C$_p$/T (Jmol$^{-1}$K$^{-2}$)")
    plt.legend()

    
if Data_view == 4:
    fig, axs = plt.subplots(2,2)
    
    #top left
    axs[0,0].scatter(data[0],data[1], c="r", s=1, label="data")
    axs[0,0].plot(temps, totaly, c="black", label="total")
    for i, ys in enumerate(debye_ys):
        if len(ys)>0:
            data_dict[f"D{i+1}"] = ys
            axs[0,0].plot(temps, ys, label =f"D{i+1}")
    for i, ys in enumerate(einstein_ys):
        if len(ys)>0:
            data_dict[f"E{i+1}"] = ys
            axs[0,0].plot(temps, ys, label =f"E{i+1}")
    if gamma:
        data_dict["Linear"] = linear_ys
        axs[0,0].plot(temps, linear_ys, label="Linear")  
    if delta:
         data_dict["Schottky"] = schottky_ys
         axs[0,0].plot(temps, schottky_ys, label="Schottky")
    axs[0,0].set_xlim([StartT, EndT])
    axs[0,0].set_ylim([0, None])
    axs[0,0].set_ylabel("C$_p$ (Jmol$^{-1}$K$^{-1}$)")
    axs[0,0].set_xlabel("Temperature (K)")


    #top right
    axs[0,1].scatter(data[0],data[1]/data[0]**3, c="r", s=1, label="data")
    axs[0,1].plot(temps, totaly/temps**3, c="black", label="total")
    for i, ys in enumerate(debye_ys):
        if len(ys)>0:
            axs[0,1].plot(temps, ys/temps**3, label =f"D{i+1}")
    for i, ys in enumerate(einstein_ys):
        if len(ys)>0:
            axs[0,1].plot(temps, ys/temps**3, label =f"E{i+1}")
    if gamma:
        axs[0,1].plot(temps, linear_ys/temps**3, label="Linear")  
    if delta:
         axs[0,1].plot(temps, schottky_ys/temps**3, label="Schottky")
    axs[0,1].set_ylim([0, None])    
    axs[0,1].set_xscale('log')
    axs[0,1].set_ylabel("C$_p$/T$^3$ (Jmol$^{-1}$K$^{-4}$)")
    axs[0,1].set_xlabel("Temperature (K)")
    axs[0,1].set_xlim([StartT, EndT])

    #bottom left
    axs[1,0].scatter(data[0]**2,data[1]/data[0], c="r", s=1, label="data")
    axs[1,0].plot(temps**2, totaly/temps, c="black", label="total")
    for i, ys in enumerate(debye_ys):
        if len(ys)>0:
            axs[1,0].plot(temps**2, ys/temps, label =f"D{i+1}")
    for i, ys in enumerate(einstein_ys):
        if len(ys)>0:
            axs[1,0].plot(temps**2, ys/temps, label =f"E{i+1}")
    if gamma:
        axs[1,0].plot(temps**2, linear_ys/temps, label="Linear")  
    if delta:
         axs[1,0].plot(temps**2, schottky_ys/temps, label="Schottky")
    axs[1,0].set_xlim([StartT**2, 40])
    axs[1,0].set_ylim([0, 0.007])
    axs[1,0].set_xlabel("Temperature$^2$ (K$^2$)")
    axs[1,0].set_ylabel("C$_p$/T (Jmol$^{-1}$K$^{-2}$)")

    #bottom right
    axs[1,1].set_facecolor("None")
    axs[1,1].scatter(data[0]**2,data[1]/data[0], c="r", s=1, label="Data")
    axs[1,1].plot(temps**2, totaly/temps, c="black", label="Total")
    for i, ys in enumerate(debye_ys):
        if len(ys)>0:
            axs[1,1].plot(temps**2, ys/temps, label =f"D{i+1}")
    for i, ys in enumerate(einstein_ys):
        if len(ys)>0:
            axs[1,1].plot(temps**2, ys/temps, label =f"E{i+1}")
    if gamma:
        axs[1,1].plot(temps**2, linear_ys/temps, label="Linear")  
    if delta:
         axs[1,1].plot(temps**2, schottky_ys/temps, label="Schottky")
    axs[1,1].set_xlim([StartT**2, EndT**2])
    axs[1,1].set_ylim([0, None])
    axs[1,1].set_xlabel("Temperature$^2$ (K$^2$)")
    axs[1,1].set_ylabel("C$_p$/T (Jmol$^{-1}$K$^{-2}$)")
    fig.tight_layout()



"""export data"""

if export_graph:
    plt.savefig('Fit.png', dpi=300)
if export_model:
    df = pd.DataFrame(data_dict)
    df.to_csv("output.csv")

