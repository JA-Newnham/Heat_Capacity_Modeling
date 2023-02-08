Modelling Heat Capacity by Jon A. Newnham

About the code:

This code allows for the modelling of Heat Capacity (Cp) data as a function of temperature using a linear combination of Debye and Einstein functions.

Heat capacity data can be fit with up to 3 Debye temperatures and 3 Einstein temperatures along with their respective refactors, as well as one linear function and one Schottky anomaly function if needed. The model and graphs produced can be outputted as .csv and .png files, respectively.

Plots can be viewed as:
Cp vs. T
Cp/T^3 vs. T
Cp/T vs T^2
or all 3 simultaneously.


How to use:
1.	Under “““Data Set”””, enter the temperature range and temperature steps to model your data with.
2.	Under “““Refined Parameters ”””, enter some initial starting parameters (i.e. Debye (Td) and Einstein (Te) temperatures) to start your modelling your data with. Parameters can be turned off entirely by setting the respective pre-factor to 0. The sum of each the Debye and Einstein prefactors (Td_p + Te_p) must equal 1.
3.	If needed, add a linear component or Schottky Anomaly term.
4.	Under “““Plotting”””, select how you would like to view your data, either as: Cp vs. T, Cp/T^3 vs. T, Cp/T vs T^2, or all 3 simultaneously.
5.	Run the code and iteratively refine your model until you’ve fit your data.
6.	Choose whether you would like to export the a .png of the graphs and/or a .csv of the model by setting the respective variable to 1 or 0.

GLHF!
