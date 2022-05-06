#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.patches as patches
sns.set(style="darkgrid")

mu_0=4*np.pi/10000000


# In[2]:


#EQUATIONS FOR CORE LOSSES (ALL CHECKED WITH MATHCAD)

#####################################################################
# Calculation of the c,x,y coefficient for the core losses equation
# INPUT: two pairs (P,B) for three different frequencies
#The units of each factor are in B-->mT and P-->kW/m3
# INPUT: initial guess for the coefficients c, x and y.
# OUTPUT: the values of c, x, y
#####################################################################

def mag_material(P1,B1,fsw1,P2,B2,fsw2,P3,B3,fsw3,c_guess,x_guess,y_guess):
    from scipy.optimize import fsolve
    
# The equations are defined in such a way that x[0]=c, x[1]=x, x[2]=y for the three points selected    
    def equations(x):
        return [P1-x[0]*fsw1**x[1]*B1**x[2], P2-x[0]*fsw2**x[1]*B2**x[2],P3-x[0]*fsw3**x[1]*B3**x[2]]

#The solver is called with the initial values for c, x, y
    x_guess=[c_guess,x_guess,y_guess]
    c_sol,x_sol,y_sol=fsolve(equations,x_guess)
    return c_sol,x_sol,y_sol
    
    
####################################################################    
# Core losses equation per unit of volume
# INPUT: losses coefficients c,x,y
# The units as provided by mag_material function
# INPUT: frequency and Bac (peak, not peak to peak)
# The units are fsw-->Hz and Bac-->mT
# OUTPUT: losses per unit of volume
# The units are kW/m3
####################################################################
def core_losses_per_vol(c,x,y,Bac,fsw):
    core_loss=c*fsw**x*Bac**y
    return core_loss


####################################################################    
# Core losses equation per unit of volume and total
# INPUT: core losses per unit of volume (prev equation) and volume
# The units losses_volume-->kW/m3 and Vol-->m3
# INPUT: Bac, frequency and volume
# The units are Bac-->mT, fsw-->Hz and Vol-->m3
# OUTPUT: total core losses
# The units are W
####################################################################

# Core losses equation
def core_losses(c,x,y,Bac,fsw,Vol):
    #The 1000 factor is for changing from kW to W
    return core_losses_per_vol(c,x,y,Bac,fsw)*Vol*1000


# In[3]:


factors=mag_material(40,80,100000,800,100,400000,60,60,200000,0.0000000000001,1,5)
core_losses(factors[0],factors[1],factors[2],126,150000,0.00000096)
factors


# In[4]:


# EQUATIONS FOR COPPER LOSSES (CHECKED WITH MATHCAD)
####################################################################
# Calculation of resistitivy based on temperature
# INPUT: Temperature for rho_base, operating temperature
# INPUT: Variation coefficient and base value of Rho
# The units has to be K or C and ohm·meter
# OUTPUT: resistitivity at the defined temperature
# The units are ohm·meter
#####################################################################

def rho_T(T0,T,T_coeff,rho_base):
    return rho_base+T_coeff*(T-T0)


#####################################################################
# Losses for a given wire
# INPUT: length and section of the wire, rms current and rho
# The units are SI
# OUTPUT: Copper losses
# The units are W
#####################################################################

def copper_losses(Irms,Awire,Lwire,rho):
    return Irms**2*rho*Lwire/Awire


# In[5]:


# EQUATION FOR MINIMUM NUMBER OF TURNS AND GAP (CHECKED WITH MATHCAD)
####################################################################
# Calculation of the minimum number of turns
# INPUT: Voltage, time the voltage is applied, Equivalent area
# INPUT: saturation magnetic field
# The units in SI
# OUTPUT: minimum number of turns
# Note: It is divided by 2 to take into account that the magnetic core
# will tend to a situation of minimum energy. Nonetheless, there are
# topologies in which a dc bias will have to be introduced.
####################################################################

def N_min_tr(V,T,Ae,Bsat):
    return round(V*T/Ae/Bsat/2+0.5)
def gap(N,Ae,mu_0,L,le,mu_r):
    return max(N**2*Ae*mu_0/L-le/mu_r,0)

def N_min_ind(Lind,Ipk,Bsat,Ae):
    return round(Lind*Ipk/Bsat/Ae+0.5)


# In[6]:


N_min_tr(50,0.5/150000,0.000039,0.44)
N_min_ind(0.000334,28,0.44,0.001170)


# In[7]:


# OPTIMUM NUMBER OF TURNS
# (CHECKED WITH MATHCAD 32 turns vs 31)
#####################################################################
# INPUT: equivalent area, Window Area, Window factor, primary and secondary share factor
# INPUT: resistivity, turn ratio, average length per turn, input voltage, duty,
# INPUT: rms primary and secondary currents, c, x, y for mag losses
# switching frequency and volume
# The units are SI
# OUTPUT: optimum number of turns
########################################################################

def N_opt_tr(Ae,WA,WF,WFpri,WFsec,rho,rt,Lm,Vin,D,Irms_pri,Irms_sec,c,x,y,fsw,Vol):
    N=1
    Copper_loss=100000000;
    Core_loss=100000000;
    Copper_loss_ant=Copper_loss+1;
    Core_loss_ant=Core_loss+1;
    while (Core_loss_ant+Copper_loss_ant>Core_loss+Copper_loss):
        Core_loss_ant=Core_loss
        Copper_loss_ant=Copper_loss
        Apri=WA*WF*WFpri/N
        Asec=WA*WF*WFsec/(N*rt)
        Copper_loss=copper_losses(Irms_pri,Apri,Lm*N,rho)+copper_losses(Irms_sec,Asec,Lm*N*rt,rho)
        Bac=Vin*D/Ae/N/fsw/2*1000
        Core_loss=core_losses(c,x,y,Bac,fsw,Vol)
        N=N+1
#whenever the condition is not satisfied, the number of turns is -1 because the final instruction added one
# and an additional -1 because the condition was not fulfilled so the previous value is selected.
    return N-2

def N_opt_ind(Ae,WA,WF,Ipk,Iavg,Irms,Lavg,Lind,rho,c,x,y,fsw,Vol): 
    N=1
    Copper_loss=100000000;
    Core_loss=100000000;
    Copper_loss_ant=Copper_loss+1;
    Core_loss_ant=Core_loss+1;
    while (Core_loss_ant+Copper_loss_ant>Core_loss+Copper_loss):
        Core_loss_ant=Core_loss
        Copper_loss_ant=Copper_loss
        A=WA*WF/N
        Copper_loss=copper_losses(Irms,A,Lavg*N,rho)
        Bac=Lind*(Ipk-Iavg)/N/Ae*10**3
        Core_loss=core_losses(c,x,y,Bac,fsw,Vol)
        N=N+1
#whenever the condition is not satisfied, the number of turns is -1 because the final instruction added one
# and an additional -1 because the condition was not fulfilled so the previous value is selected.
    return N-2


# In[8]:


rho=0.0000000171
N_opt_tr(0.000039,0.000020,0.3,0.5,0.5,0.0000000171,0.124,0.03800,50,0.5,0.314,3.7,factors[0],factors[1],factors[2],150000,0.00000096)
N_opt_ind(0.001170,0.000928,0.3,28,500/18,27.7,0.49500,0.000334,rho,6.198*10**-13,1.721,2.73,30000,0.0002154000)


# In[9]:


# Cores data
cores=pd.read_excel(r'D:\Trabajo\HC MathCad\Magnetics\trafos.xls',sheet_name="Cores")
cores.head()


# In[10]:


#OPTIMUM DESIGN OF THE TRANSFORMER
# Lmag not checked
######################################################################################
# Optimum design of the transformer
######################################################################################
# INPUT: all variables for the previously defined equations
# The units are as defined in the previous functions
# Vtr-->V,  fsw-->Hz, Bsat-->T , c,x,y-->defined by mag_material function, 
# Lmag-->H (1 for not required), Irms-->Arms, rho-->ohm·meter, To, Tmax-->ºC
# OUTPUT: core part, copper losses, core losses, final Temp, Nopt, Nmin, Nsel
# The units are SI
######################################################################################

def opt_trans(Vtr,Dtr,fsw,Bsat,mu_r,c,x,y,Lmag,WF,Wpri,Wsec,Irms_pri,Irms_sec,rho,To,Tmax):
   
    Tsw=1/fsw

    for i in cores.index:
        #For consistency with MathCad, which uses end as closing indentifier. Not necessary in Python
        # and used to avoid an error.
        if cores['Part'][i]=="end":
            break
        #Calculation of the minimum number of turns (only dependent on the core and the switching frequency)    
        Nmin_i=N_min_tr(Vtr,Tsw*Dtr,cores['Ae (m2)'][i],Bsat)
        Nopt_i=N_opt_tr(cores['Ae (m2)'][i],cores['Aw (m2)'][i],WF,Wpri,Wsec,rho,rtr,cores['lm (m)'][i],Vtr,Dtr,Irms_pri,Irms_sec,c,x,y,fsw,cores['Ve (m3)'][i])
        if (Nopt_i>Nmin_i):
            Nsel_i=Nopt_i
        else:
            Nsel_i=Nmin_i
        Bac=Vtr*Dtr/cores['Ae (m2)'][i]/Nsel_i/fsw/2*1000
        core_loss=core_losses(c,x,y,Bac,fsw,cores['Ve (m3)'][i])
        copper_loss_pri=copper_losses(Irms_pri,cores['Aw (m2)'][i]*WF*Wpri/Nsel_i,cores['lm (m)'][i]*Nsel_i,rho)
        copper_loss_sec=copper_losses(Irms_sec,cores['Aw (m2)'][i]*WF*Wsec/Nsel_i/rtr,cores['lm (m)'][i]*Nsel_i*rtr,rho)
        Tend=To+53*(cores['Ve (m3)'][i]*1000000)**(-0.53)*(core_loss+copper_loss_pri+copper_loss_sec)
        g=gap(Nsel_i,cores['Ae (m2)'][i],mu_0,Lmag,cores['le (m)'][i],mu_r)
        # Quedaría calcular y devolver esto
        Lmag=Nsel_i**2*mu_0*cores['Ae (m2)'][i]/(g+cores['le (m)'][i]/mu_r)
        if (Tend<Tmax):
            break
    return (cores['Part'][i],copper_loss_sec+copper_loss_pri,core_loss,Tend,Nopt_i,Nmin_i,Nsel_i,g,Lmag)


# In[11]:


#Validation for previous function
fsw=150000
Vtr=50
Dtr=0.5
Bsat=0.44
WF=0.3
Wpri=0.5
Wsec=0.5
Irms_pri=0.314
Irms_sec=3.7
rho=0.0000000171
To=10
(c,x,y)=mag_material(40,80,100000,800,100,400000,60,60,200000,0.0000000000001,1,5)
Tsw=1/fsw
rtr=0.124
Tmax=20
mu_r=2200
opt_trans(Vtr,Dtr,fsw,Bsat, mu_r,c,x,y,1,WF,Wpri,Wsec,Irms_pri,Irms_sec,rho,To,Tmax)


# ![image.png](attachment:image.png)

# ![image.png](attachment:image.png)

# In[12]:


# OPTIMUM DESIGN OF THE INDUCTOR
#(CHECKED WITH MATHCAD)
def opt_inductor(fsw,Bsat,mu_r,c,x,y,Lind,WF,Iavg,Ipk,Irms,rho,To,Tmax):
   
    Tsw=1/fsw

    for i in cores.index:
        #For consistency with MathCad, which uses end as closing indentifier. Not necessary in Python
        # and used to avoid an error.
        if cores['Part'][i]=="end":
            break
        #Calculation of the minimum number of turns (only dependent on the core and the switching frequency)    
        Nmin_i= N_min_ind(Lind,Ipk,Bsat,cores['Ae (m2)'][i])
        Nopt_i= N_opt_ind(cores['Ae (m2)'][i],cores['Aw (m2)'][i],WF,Ipk,Iavg,Irms,cores['lm (m)'][i],Lind,rho,c,x,y,fsw,cores['Ve (m3)'][i]) 
        if (Nopt_i>Nmin_i):
            Nsel_i=Nopt_i
        else:
            Nsel_i=Nmin_i
       
        Bac=Lind*(Ipk-Iavg)/Nsel_i/cores['Ae (m2)'][i]*10**3
        core_loss=core_losses(c,x,y,Bac,fsw,cores['Ve (m3)'][i])
        copper_loss=copper_losses(Irms,cores['Aw (m2)'][i]*WF/Nsel_i,cores['lm (m)'][i]*Nsel_i,rho)
        Tend=To+53*(cores['Ve (m3)'][i]*1000000)**(-0.53)*(core_loss+copper_loss)
        g=gap(Nsel_i,cores['Ae (m2)'][i],mu_0,Lind,cores['le (m)'][i],mu_r)
        if (Tend<Tmax):
            break
    return (cores['Part'][i],copper_loss,core_loss,Tend,Nopt_i,Nmin_i,Nsel_i,g)


# In[13]:


opt_inductor(30000,0.44,2200,6.198*10**-13,1.721,2.73,0.000334,0.3,500/18,28,27.7,rho,40,100)

