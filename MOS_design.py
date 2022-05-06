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


# In[2]:


#validated
def Rdson_temp(Rdson_25,Rdson_125,Tj):
    alpha=((Rdson_125/Rdson_25)**(1/100)-1)*100
    Rdson_Tj=Rdson_25*(1+alpha/100)**(Tj-25)
    return Rdson_Tj


# In[3]:


#validated
def Coss_sw_losses(Coss,Vds,fsw):
    return Coss*Vds**2*fsw/2


# In[4]:


#validated
def gate_sw_losses(Qg,Vdr,fsw):
    return Qg*Vdr*fsw


# In[5]:


#validated
def sw_losses(Vds,Ids,fsw,Qg,Idr):
#based on fairchild application note
    return Vds*Ids/2*fsw*Qg/Idr


# In[6]:


#validated
def cond_losses(Rdson,Irms):
    return Rdson*Irms**2


# In[7]:


#validated
def Id_max_T(Id_max_25,Id_max_100,Tj):
    return Id_max_25+(Id_max_100-Id_max_25)/(100-25)*(Tj-25)


# In[8]:


# MOSFET data
MOSFETs=pd.read_excel(r'D:\Trabajo\HC MathCad\Semiconductores\Space Components.xls',sheet_name="MOSFETs")
MOSFETs.tail()


# In[58]:


#Validated

def opt_MOSFET(Vds,Id,Irms,ZVS,ZCS,Vds_sw_on,Vds_sw_off,Ids_sw_on,Ids_sw_off,Vdr,Idr,fsw,To):
    Pmos_opt=100000
    sel=666
    for i in MOSFETs.index:
     
        #For consistency with MathCad, which uses end as closing indentifier. Not necessary in Python
        # and used to avoid an error.
        if (MOSFETs['Part'][i]=="end"):
            break
        if (MOSFETs['Coss (pF)'][i]=="NP" or MOSFETs['Qg (nC)'][i]=="NP" or MOSFETs['Rth jc (oC/W)'][i]=="NP"):
            continue
        Coss=MOSFETs['Coss (pF)'][i]*10**-12
        Qg=MOSFETs['Qg (nC)'][i]*10**-9
        
        if (MOSFETs['Vds (V)'][i]<Vds):
            continue
       #Losses calculation     
        Pmos_gate=gate_sw_losses(Qg,Vdr,fsw)
        
        if (ZVS==1):
            Pmos_coss=Coss_sw_losses(Coss,Vds_sw,fsw)
        else:
            Pmos_coss=0
            
        if (ZVS==1 and ZCS==1):
            Pmos_sw=0
        elif (ZVS==0 and ZCS==0):
            Pmos_sw=sw_losses(Vds_sw_on,Ids_sw_on,fsw,Qg,Idr)+sw_losses(Vds_sw_off,Ids_sw_off,fsw,Qg,Idr)
        elif (ZVS==1 and ZCS==0):
            Pmos_sw=sw_losses(Vds_sw_off,Ids_sw_off,fsw,Qg,Idr)
        else:
            Pmos_sw=sw_losses(Vds_sw_on,Ids_sw_on,fsw,Qg,Idr)
        
        Tj=To
        Tj_end=Tj+0.2
        while (Tj_end<Tj-0.1 or Tj_end>Tj+0.1):
            Tj=Tj+0.1
            Rdson=Rdson_temp(MOSFETs['Rdson (ohms @25oC)'][i],MOSFETs['Rdson (ohms @125oC)'][i],Tj)
            Pmos_cond=cond_losses(Rdson,Irms)
            Pmos_tot=Pmos_cond+Pmos_sw+Pmos_gate+Pmos_coss
            Tj_end=To+Pmos_tot*MOSFETs['Rth jc (oC/W)'][i]
            if (Tj_end>MOSFETs['Tj_max (oC)'][i]):
                break
        
        if (Tj_end>MOSFETs['Tj_max (oC)'][i]):
            continue
        
        if ((Id_max_T(MOSFETs['Id (A @25oC)'][i],MOSFETs['Id (A @100oC)'][i],Tj_end)<Id)):
            continue
        
        if Pmos_tot<Pmos_opt:
            Pmos_opt=Pmos_tot
            Pmos_sw_opt=Pmos_sw
            Pmos_coss_opt=Pmos_coss
            Pmos_gate_opt=Pmos_gate
            Pmos_cond_opt=Pmos_cond
            Tj_opt=Tj_end
            sel=i
    if sel==666:
        return "NMV"
    return (MOSFETs['Part'][sel],Pmos_cond_opt,Pmos_sw_opt,Pmos_coss_opt,Pmos_gate_opt,Pmos_opt,Tj_opt)
            
        


# In[57]:


opt_MOSFET(100,10,12,0,0,100,100,30,30,20,1,100000,40)

