# Particle-Charge-Distribution-Calculator
## Overview
This code calculate the variation of particle charge distribution.  It is assumed that the particle charge distribution is changed by traveling the ion existing space (e.g., ionizer but this calcualtion is not only applicable to high ion concentration cases).  The initial charge distribution of $z$ charged particle is given as $N_z(d_p,0)$ and it is changed to $N_z(d_p,t)$ by colliding with ions for residence time in the charging area ($t$).  

## Theory
### Population balance equaiton
The elementary reaction of the charging process is the second order reaction as:  
$$A^z+B^{\pm}{\rightarrow}A^{z\pm1}$$  
where $A^z$ is $z$ charged particle and $B^{\pm}$ is the positive or  negative ions.  In fact, exact right hand is $AB^{z\pm1}$ but ion $B$ is generally omitted since it negligibly smaller than the particle $A$. Based on this reaction, the populaiton balance equaiton is given as:  
$${dN_{z}\over{dt}}=-\beta_{z}^{+}N_{+}N_{z}-\beta_{z}^{-}N_{-}N_{z}+\beta_{z+1}^{-}N_{-}N_{z+1}+\beta_{z-1}^{+}N_{+}N_{z-1}$$  
$N_z$ is the concentraiton of charged particle $A^z$ and $N_{+}/N_{-}$ are the positive/negative charged ion concentrations.  $\beta_{z}^{\pm}$ is the collision rate coefficient of $z$ charged particle and ion ($A^z$ and $B^{\pm}$ ).  While concentrations are experimental conditions, collision rate coefficient is theoretically estimatable and that calculation process is shown in following section.
### Collision rate coefficient calculation ($\beta$)
#### Limiting sphere theory (Fuchs)
#### Modified limiting sphere theory (Hoppel & Frick)
Under preparation
#### LD based equation
### Time evolution

## Usage
### Use execute file
### Run on the Python

# Author
* Tomoya Tamadate
* University of Minnesota
* tamalab0109[at]gmail.com

