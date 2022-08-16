# MTL_Simulator
An object-oriented python module for simulating non-uniform cables using many-conductor transmission line theory.

## Introduction
This program allows the simulation of cables consisting of uncoated or coated wires. It uses an expansion of classic many-conductor transmission line theory as [published](https://doi.org/10.1063/5.0059393 'Surface wave transmission line theory for single and many wire systems') as well as unpublished results which can be found in my thesis (to be uploaded to the Cambridge University repository soon).  

## Dependencies
The program was written in python 3.7.7 and uses python's numpy (tested with v.1.21.2), matplotlib (v.3.4.3) and math (v.1.2.1) modules. Additionally, it uses the  scikit-rf toolbox (v.0.15.4).

## Running a simulation
The program has a few examples for running simulations on non-uniform twisted pairs. Just copy one of the example files into the main folder to run. To write your own simulation, the following steps must be made. 
1. import the em_env module and define the electromagnetic environment i.e. the frequencies, relative permittivity and loss tangent of the surrounding medium. Finally create an environment object with the EmEnv command. 
2. import the geometry module and create the geometry of your cable by defining cable objects and adding them to a Geometry object.
3. import the simulation module and initialise with the geometry and environment objects. Define a number points along the cable which will be analysed.
4. Run the simulation to get the abcd parameters. Use the s-parameter class to convert them to S-parameters and scikit-rf for furhter post-processing. 

## Notes
I do not take any responsability for the accuracy or inaccuracy of simulation results. 
Long cables can take a long time to simulate. It may prove easier to simulate a small section and then concatenate the S-parameters using scikit-rf. 
Note that the excitation of the cables needs to be explicitly specified later on. In general, converting to S-parameters will give the results for exciting individual wires with 50 Ohms, usually a pretty poor match!
