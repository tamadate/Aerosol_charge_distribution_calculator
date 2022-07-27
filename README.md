# Aerosol charge distribution calculator
## Overview
This code calculate the variation of particle charge distribution.  It is assumed that the particle charge distribution is changed by colliding with positive or negative ions during traveling the ion existing space.  For example, charging process via ion generator (ionizer) as shown in follwing figure is typical calculation case but it should be noted that applicable ion concentration of this code is not limited to high concentration.  Also, this code can be applied to the any positive/negative ion balance (i.e., both of the unipolar and bipolar is calculatable).
<br>
<br>  
![IonizationProcess](https://user-images.githubusercontent.com/75816343/181148783-0b59f3d4-7f05-4a97-83fd-7ecc74bad322.png)
Figure 1. Schematic diagram of charging process.  Particle have a arbitrary charge distribution before passing through the charger (left) and the it is varied at after passing the charger (right).  In this schematic, the after distribution is shifted to positive side because the ion polarity in the charger is positive biased.
<br>
<br>
## Theory
### Population balance equaiton
The elementary reaction of the charging process is the second order reaction as:  
$$A^z+B^{\pm}{\rightarrow}A^{z\pm1}$$  
where $A^z$ is $z$ charged particle and $B^{\pm}$ is the positive or  negative ions.  In fact, exact right hand is $AB^{z\pm1}$ but ion $B$ is generally omitted since it negligibly smaller than the particle $A$. Based on this reaction, the populaiton balance equaiton is given as:  
$${dN_{z}\over{dt}}=-\beta_{z}^{+}N_{+}N_{z}-\beta_{z}^{-}N_{-}N_{z}+\beta_{z+1}^{-}N_{-}N_{z+1}+\beta_{z-1}^{+}N_{+}N_{z-1}$$  
$N_z$ is the concentraiton of charged particle $A^z$ and $N_{+}/N_{-}$ are the positive/negative charged ion concentrations.  $\beta_{z}^{\pm}$ is the collision rate coefficient of $z$ charged particle and negatively or positively charged ion ($A^z$ and $B^{\pm}$ ).  To simplify the calculation, in this calculation code, the particle concentration is normalized by the total particle concentration $N_0=\Sigma^\infty_{i=-\infty}N_i$.
In this equaiton, the concentrations are experimental conditions and collision rate coefficients are theoretically estimatable parameter which  calculation method is shown in following part.
### Collision rate coefficient ($\beta$)
Currently, this code support three different collision rate coefficient calculaiton method (Fuchs, Hoppel & Frick, and LD based).
* Limiting sphere model (Fuchs)
Fuchs proposed follwing equaiton which called limiting sphere theory:
$$\beta_z^{\pm} = {{{\pi}c_{\pm}p{\delta}^2exp[-{\phi}({\delta})/k_bT]} \over {1+exp[-{\phi}({\delta})/k_bT]{{c_{\pm}p{\delta}^2} \over {4D_{\pm}a}}}\int_0^{a/{\delta}}exp(-{\phi}(a/x)/k_bT)dx}$$  
$${\phi}(r)={z_{ion}ze^2 \over 4{\pi}\epsilon_0r} - {\epsilon-1 \over \epsilon+1}{e^2a^3 \over 8{\pi}\epsilon_0r^2(r^2-a^2)}$$  
where, $c_{\pm}$ is the ions mean thermal velocity, $p$ is the collision probability, $\delta$ is the radius of the limiting sphere, $k_b$ is Boltzumann constant, $T$ is the temperature, $D_{\pm}$ is the ion diffusion coefficient, $a$ is the particle radius, $\phi(r)$ is the ion-particle potential function which generally defined as a sum of the electric and image potential, $z_{ion}$ is the ion's number of charge (${\pm}1$), $r$ is ion-particle distance, and $\epsilon_0$ is the dielectric constant of vacuum.  Collision probability $p$ is difined as a number of atoms per total injecting atoms to the limiting sphere and it is calculated in this theory from:
$$p={[b(min) \over \delta]}^2$$
$$b^2=\rho_m(1+{2k_bT \over 3}[\phi(\delta)-\phi(\rho_m)])
You can find more detail from this link.
* Modified limiting sphere model (Hoppel & Frick)
Under preparation

* LD based equation
LD based equation is given as:
### Time evolution
In this code, 4th order Runge-Kutta method is applied as a time evolution method.
## Usage
### Use execute file
1. Download files and open chargeDistributionCalculator.exe from downloaded directory then following window is displayed on your computer:

2. Input physical properties: (1) section is the physical properties.  General air conditions are automatically filled out when you push "Set air normal conditions".  Your values can be typed in each boxes as an input value.  It should be noted that when you change the particle diameter, the particle diffusion coefficient is needed to recalculate by "Calculate diffusion coefficient" or you need to retype in an appropriate value.
3. Set initial concentrations: (2) section is the normalized initial concentrations of each charge state $N_z/N_0$.  "Set C0" set 1 as a $z=0$ initial concentration and 0 as for the other charge statements which mean all particle is neutral at the starting time.  It can be typed in your own initial values.
4. Set collision rate coefficient: (3) section is collision rate coefficients.  Left column is for the collision with "negative" ions ($\beta_z^-$) right column is for the "positive" ions ($\beta_z^+$).  "Set Beta Fuchs" generate the values based on Fuchs' theory and the image potential effect can be controled from the check box (exclude image potential if you don't check the image force check box).  "Set Beta CG" is for the LD based equation.  You are also able to type in your own values.  
5. Run "Solve dn/dt" and the time evolution of the charge distribution is displayed (only $z=-2$ to $+2$).
6. The output file "a.dat" is created in your running directory.
### Run on the Python
1. Download files and run the chargeDistributionCalculator.py from the src.
2. Follow execute one.
## Author
* Dr. Tomoya Tamadate
* [LinkedIn](https://www.linkedin.com/in/tomoya-tamadate-953673142/)/[ResearchGate](https://www.researchgate.net/profile/Tomoya-Tamadate)/[Google Scholar](https://scholar.google.com/citations?user=XXSOgXwAAAAJ&hl=ja)
* University of Minnesota
* tamalab0109[at]gmail.com
