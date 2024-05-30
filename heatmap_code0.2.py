


# %% [markdown]
# HEATMAP

# %%
D1_r=np.linspace(D1_min_sh,D1_max_sh,size_D1)
D2_r=np.linspace(D2_min_sh,D2_max_sh,size_D2)

CSV = pd.read_csv("Measurements\Bbheatmap3.csv", delimiter=';')
D1_r = np.array(CSV['d1'])
D2_r=np.array(CSV['d2'])
from itertools import product
from typing import ValuesView
#values = np.array(list(product(D1_r,D2_r))) #to get all possible combinations. this if using linspace. if loading from csv use the one bleow
values=np.column_stack((D1_r,D2_r))#I already have all combination, so i only stack both columns

df_hm=pd.DataFrame(values)
df_hm.columns=['D1','D2']
D1_r=df_hm.loc[:,'D1']
D2_r=df_hm.loc[:,'D2']
phi=0
phi_bob=phi*np.ones(df_hm.shape[0])
fsw=180e3
T=1/fsw
P=500
Vbus=100
L_ind=83.5e-6
# df_hm['D1']=D1_r
# df_hm['D2']=D2_r
df_hm['Gain']=D1_r/(1-D2_r)
df_hm['phi']=phi_bob*np.ones(df_hm.shape[0])
df_hm['T']=T*np.ones(df_hm.shape[0])
df_hm['fsw']=fsw*np.ones(df_hm.shape[0])
df_hm['L']=L_ind*np.ones(df_hm.shape[0])
df_hm['Vin']=Vbus*np.ones(df_hm.shape[0])/df_hm.loc[:,'Gain']
df_hm['Vbus']=Vbus*np.ones(df_hm.shape[0])
df_hm['P']=P*np.ones(df_hm.shape[0])
df_hm['Iin(A)']=P/(Vbus/df_hm.loc[:,'Gain'])



c=currents_L(df_hm['L'],df_hm['Vin'],Vbus,df_hm['Iin(A)'],D1_r,D2_r,phi_bob,T)
interv=intervals(D1_r,D2_r,phi_bob,T)
df_hm_los=df_hm.copy()

Irms_M1,Irms_M2,Irms_M3,Irms_M4=Irms_Mi(D1_r,D2_r,phi_bob,fsw,interv[2],interv[5],interv[8],interv[11],interv[1],c[0],c[1],c[2],c[3])
aux_M1,Temp_M1=Losses_M1(D1_r,Vin_f(Vbus,D1_r,D2_r),c[3],c[1],fsw,M1_index,Vdrvr,Idrvr,Irms_M1,aging_n,Rjc_M1,Rca_M1,Tini,Tjmax_M1,df_hm_los)
aux_M2,Temp_M2=Losses_M2(D1_r,Vin_f(Vbus,D1_r,D2_r),c[1],c[3],fsw,M2_index,Vdrvr,Idrvr,Irms_M2,aging_n,Rjc_M2,Rca_M2,Tini,Tjmax_M2,df_hm_los)
aux_M3,Temp_M3=Losses_M3(D2_r,Vbus,c[2],c[0],fsw,M3_index,Vdrvr,Idrvr,Irms_M3,aging_n,Rjc_M3,Rca_M3,Tini,Tjmax_M3,df_hm_los)
aux_M4,Temp_M4=Losses_M4(D2_r,Vbus,c[0],c[2],fsw,M4_index,Vdrvr,Idrvr,Irms_M4,aging_n,Rjc_M4,Rca_M4,Tini,Tjmax_M4,df_hm_los)
aux_protM3,Temp_protM3=Losses_prot(M3prot_index,Irms_M3,aging_n,Rjc_M3prot,Rca_M3prot,Tini,Tjmax_M3prot,df_hm_los)
aux_protM1,Temp_protM1=Losses_prot(M1prot_index,Irms_M1,aging_n,Rjc_M1prot,Rca_M1prot,Tini,Tjmax_M1prot,df_hm_los)

Iac_ind=c[-1].max(axis=1)-c[-1].min(axis=1) #Validated with other implementation. Works
df_hm_los['Iac_ind']=Iac_ind

Iavg_ind=Iind_avg(c[3],c[0],c[1],c[2],interv[2]-interv[1],interv[5]-interv[4],interv[8]-interv[7],interv[11]-interv[10])

Irms1=Irms_sec(1/T,interv[2],interv[1],c[0],c[3])
Irms2=Irms_sec(1/T,interv[5],interv[4],c[1],c[0])
Irms3=Irms_sec(1/T,interv[8],interv[7],c[2],c[1])
Irms4=Irms_sec(1/T,interv[11],interv[10],c[3],c[2])
Irms_ind=np.sqrt(Irms1**2+Irms2**2+Irms3**2+Irms4**2)

df_hm_los['Irms_ind']=Irms_ind
df_hm_los['Ipk_ind']=c[-1].max(axis=1)
df_hm_los['Iavg']=(c[-1].max(axis=1)+c[-1].min(axis=1))/2
aux_L=Losses_ind(c_ind,x_ind,y_ind,Vol_ind,fsw,L_ind,Iac_ind,N_ind,Ae_ind,Irms_ind,Awire_ind,Lwire_ind,rho)

aux=aux_M1+aux_M2+aux_M3+aux_M4+aux_L+aux_protM1+aux_protM3
losses=aux

df_hm_los['losses M1']=aux_M1
df_hm_los['losses M2']=aux_M2
df_hm_los['losses M3']=aux_M3
df_hm_los['losses M4']=aux_M4
df_hm_los['losses L']=aux_L
df_hm_los['losses M1 prot']=aux_protM1
df_hm_los['Losses M3 prot']=aux_protM3

df_hm_los['Losses']=losses
df_hm_los['Temp_M1']=Temp_M1
df_hm_los['Temp_M2']=Temp_M2
df_hm_los['Temp_M3']=Temp_M3
df_hm_los['Temp_M4']=Temp_M4
df_hm_los['Temp_M1prot']=Temp_protM1
df_hm_los['Temp_M3prot']=Temp_protM3

df_th_meas = pd.DataFrame({
    'D1': D1_r,
    'D2': D2_r,
    'Theoretical losses': losses,

})

# %%
df_hm_los

# %%
heat=df_hm_los.loc[:,['D1','D2','Losses']]
heat=heat.set_index('D2')
heat=heat.pivot(columns='D1')
heat_val=heat.values
D1_x=heat.columns.get_level_values(1)
D2_y=heat.index.values

fig, ax = plt.subplots()
#Plot the surface.
gr=ax.contourf(D1_x, D2_y, heat_val, 50, cmap='hot')#,vmin=6.9,vmax=30)  #'RdYlBu'
cbar=fig.colorbar(gr)
# ticklabs = cbar.ax.get_yticklabels()
# cbar.ax.set_yticklabels(ticklabs, fontsize=20)
# fig.canvas.toolbar_visible = True
# fig.canvas.header_visible = True
# fig.canvas.resizable = True

# Black lines -> max and min gain
#ax.plot(D1_r,D2_f(Vsa_max,Vbus,D1_r),'k')
#ax.plot(D1_r,D2_f(Vsa_min,Vbus,D1_r),'k')

# White lines-> limit the duty cyles due to dead times

#ax.plot(D1_r,D2_f(Vbus,Lim_sup_D1*Vbus,D1_r),'w')
#ax.plot(D1_f(Vbus*(1-Lim_inf_D2),Vbus,D2_r),D2_r,'w')


#ax.plot(D1_r,D2_r,'g')
plt.xlim(0.64,1)
plt.ylim(0,0.4)
# plt.xticks(fontsize = 20)
# plt.yticks(fontsize = 20)
# f = plt.figure()
# f.set_figwidth(20)
# f.set_figheight(10)


plt.show()




# %%
heat_val
 

# %% [markdown]
# MEASUREMETNS

# %%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read CSV file
CSVM = pd.read_csv("Measurements\Bbheatmap3.csv", delimiter=';')
vout = CSV['vout']
iout = CSV['iout']
vin = CSV['vin']
iin = CSV['iin']
d1 = CSV['d1']
d2 = CSV['d2']

# Calculate additional parameters
pin = vin * iin
pout = iout * vout
losses = pin - pout
CSV['pin'] = pin
CSV['pout'] = pout
CSV['losses'] = losses
CSV = CSV[['vin', 'iin', 'iout', 'vout', 'd1', 'd2', 'pin', 'pout', 'losses']]
CSV.to_excel("output_data2.xlsx", index=False)
# Reshape losses array
losses_rearr = np.reshape(losses.to_numpy(), (len(np.unique(d1)), len(np.unique(d2))))
losses_rearr=losses_rearr.T
# Create contour plot
plt.figure()
contour = plt.contourf(np.unique(d1), np.unique(d2),losses_rearr,70,cmap='hot', linewidth=0)
plt.colorbar()

# Plot additional lines
d1b = np.linspace(np.min(np.unique(d1)), np.max(np.unique(d1)))
d2b = 1 - (d1b * 80 / 100)
plt.plot(d1b, d2b, 'k', linewidth=2)
d2b = 1 - (d1b * 150 / 100)
plt.plot(d1b, d2b, 'k', linewidth=2)

# Set axis limits
plt.axis([0.64,1,0, 0.4])

# Set colorbar range
#plt.clim(6.9, 30)

# Show the plot
plt.show()

df_th_meas['Meas. losses']=losses


# %%
#df_th_meas=df_th_meas.drop(['Loss diff_norm','Meas. losses_norm','Loss diff_norm','Meas. losses2','Loss diff2','Loss diff_norm2'],axis=1)
df_th_meas['Loss diff']=df_th_meas['Meas. losses']-df_th_meas['Theoretical losses']
df_th_meas['Loss diff_norm']=df_th_meas['Meas. losses']/df_th_meas['Meas. losses'].max()
df_th_meas['Theoretical losses_norm']=df_th_meas['Theoretical losses']/df_th_meas['Meas. losses'].max()
df_th_meas['Meas. losses_norm']=df_th_meas['Meas. losses']/df_th_meas['Meas. losses'].max()
df_th_meas['Loss diff_norm']=df_th_meas['Loss diff']/df_th_meas['Meas. losses'].max()
df_th_meas['Meas. losses2']=(df_th_meas['Meas. losses']-df_th_meas['Loss diff'].min())
df_th_meas['Meas. losses2_norm']=(df_th_meas['Meas. losses']-df_th_meas['Loss diff'].min())/df_th_meas['Meas. losses2'].max()
df_th_meas['Loss diff2']=(df_th_meas['Loss diff']-df_th_meas['Loss diff'].min())
df_th_meas['Theoretical losses_norm2']=df_th_meas['Theoretical losses']/df_th_meas['Meas. losses2'].max()
df_th_meas['Loss diff_norm2']=df_th_meas['Loss diff2']/df_th_meas['Meas. losses2'].max()
df_th_meas['Loss diff_norm3']=df_th_meas['Loss diff2']/df_th_meas['Meas. losses2'].max()
df_th_meas['Theoretical losses_norm3']=df_th_meas['Theoretical losses']/df_th_meas['Theoretical losses'].max()
df_th_meas['Meas. losses_norm3']=(df_th_meas['Meas. losses']-df_th_meas['Loss diff'].min())/df_th_meas['Meas. losses2'].max()
print(df_th_meas['Loss diff'].min())
df_th_meas

# %%
# losses_rearro = np.reshape(df_th_meas['Theoretical losses'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
# losses_rearro=losses_rearro.T
# contour = plt.contourf(np.unique(df_th_meas['D1']),np.unique(df_th_meas['D2']),losses_rearro,70,cmap='hot', linewidth=0)
# plt.colorbar()

# # Set axis limits
# plt.axis([0.64,1,0, 0.4])

# # Set colorbar range
# plt.clim(6.9, 30)

losses_rearro2 = np.reshape(df_th_meas['Loss diff'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro2=losses_rearro2.T
contour = plt.contourf(np.unique(df_th_meas['D1']),np.unique(df_th_meas['D2']),losses_rearro2,70,cmap='hot', linewidth=0)
plt.axis([0.64,1,0, 0.4])
plt.colorbar()

# Set colorbar range
#plt.clim(6.9, 30)

# %%
#reploteamos todo pero con la escala de colorbar del measured
#Teorico
losses_rearro = np.reshape(df_th_meas['Theoretical losses'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro=losses_rearro.T
contour = plt.contourf(np.unique(df_th_meas['D1']),np.unique(df_th_meas['D2']),losses_rearro,70,cmap='hot', linewidth=0)


# %%
import matplotlib.pyplot as plt
import numpy as np

# Assuming you have df_th_meas DataFrame and necessary data

# Create a figure with two subplots side by side
fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

# Plot the first contourf in the first subplot
losses_rearro = np.reshape(df_th_meas['Theoretical losses'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro = losses_rearro.T
contour1 = axs[0].contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro, 70, cmap='hot', linewidth=0)
axs[0].set_title('Theoretical Losses')
axs[0].set_xlabel('D1')
axs[0].set_ylabel('D2')

# Plot the second contourf in the second subplot
losses_rearro2 = np.reshape(df_th_meas['Meas. losses'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro2 = losses_rearro2.T
contour2 = axs[1].contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro2, 70, cmap='hot', linewidth=0)
axs[1].set_title('Measured Losses')
axs[1].set_xlabel('D1')
axs[1].set_ylabel('D2')

# Determine the common color limits based on the maximum values
max_value = max(losses_rearro.max(), losses_rearro2.max())

# Set color limits for both subplots
contour1.set_clim(0, max_value)
contour2.set_clim(0, max_value)

# Create a single colorbar on the right side
cbar_ax = fig.add_axes([1, 0.15, 0.02, 0.7])  # [x, y, width, height]
cbar = fig.colorbar(contour2, cax=cbar_ax, orientation='vertical')

# Add a title to the colorbar
cbar.set_label('Losses (W)', rotation=270, labelpad=15)

# Adjust layout for better spacing
plt.tight_layout()

# Show the plots
plt.show()



# %%
import matplotlib.pyplot as plt
import numpy as np

# Assuming you have df_th_meas DataFrame and necessary data

# Create a figure with two subplots side by side
fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

# Plot the first contourf in the first subplot
losses_rearro = np.reshape(df_th_meas['Theoretical losses'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro = losses_rearro.T
contour1 = axs[0].contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro, 70, cmap='hot', linewidth=0)
axs[0].set_title('Theoretical Losses')
axs[0].set_xlabel('D1')
axs[0].set_ylabel('D2')

# Plot the second contourf in the second subplot
losses_rearro2 = np.reshape(df_th_meas['Meas. losses2'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro2 = losses_rearro2.T
contour2 = axs[1].contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro2, 70, cmap='hot', linewidth=0)
axs[1].set_title('Measured Losses')
axs[1].set_xlabel('D1')
axs[1].set_ylabel('D2')

# Determine the common color limits based on the maximum values
max_value = max(losses_rearro.max(), losses_rearro2.max())

# Set color limits for both subplots
contour1.set_clim(0, max_value)
contour2.set_clim(0, max_value)

# Create a single colorbar on the right side
cbar_ax = fig.add_axes([1, 0.15, 0.02, 0.7])  # [x, y, width, height]
cbar = fig.colorbar(contour2, cax=cbar_ax, orientation='vertical')

# Add a title to the colorbar
cbar.set_label('Losses (W)', rotation=270, labelpad=15)

# Adjust layout for better spacing
plt.tight_layout()

# Show the plots
plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np

# Assuming you have df_th_meas DataFrame and necessary data

# Create a figure with two subplots side by side
fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

# Plot the first contourf in the first subplot
losses_rearro = np.reshape(df_th_meas['Theoretical losses_norm'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro = losses_rearro.T
contour1 = axs[0].contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro, 70, cmap='hot', linewidth=0)
axs[0].set_title('Theoretical Losses')
axs[0].set_xlabel('D1')
axs[0].set_ylabel('D2')

# Plot the second contourf in the second subplot
losses_rearro2 = np.reshape(df_th_meas['Meas. losses2_norm'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro2 = losses_rearro2.T
contour2 = axs[1].contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro2, 70, cmap='hot', linewidth=0)
axs[1].set_title('Measured Losses')
axs[1].set_xlabel('D1')
axs[1].set_ylabel('D2')

# Determine the common color limits based on the maximum values
max_value = max(losses_rearro.max(), losses_rearro2.max())

# Set color limits for both subplots
contour1.set_clim(0, max_value)
contour2.set_clim(0, max_value)

# Create a single colorbar on the right side
cbar_ax = fig.add_axes([1, 0.15, 0.02, 0.7])  # [x, y, width, height]
cbar = fig.colorbar(contour2, cax=cbar_ax, orientation='vertical')

# Add a title to the colorbar
cbar.set_label('Losses (W)', rotation=270, labelpad=15)

# Adjust layout for better spacing
plt.tight_layout()

# Show the plots
plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np

# Assuming you have df_th_meas DataFrame and necessary data

# Create a figure with two subplots side by side
fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)

# Plot the first contourf in the first subplot
losses_rearro = np.reshape(df_th_meas['Theoretical losses_norm3'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro = losses_rearro.T
contour1 = axs[0].contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro, 100, cmap='hot', linewidth=0)
axs[0].set_title('Theoretical Losses')
axs[0].set_xlabel('D1')
axs[0].set_ylabel('D2')

# Plot the second contourf in the second subplot
losses_rearro2 = np.reshape(df_th_meas['Meas. losses_norm'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro2 = losses_rearro2.T
contour2 = axs[1].contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro2, 100, cmap='hot', linewidth=0)
axs[1].set_title('Measured Losses')
axs[1].set_xlabel('D1')
axs[1].set_ylabel('D2')

# Determine the common color limits based on the maximum values
max_value = max(losses_rearro.max(), losses_rearro2.max())

# Set color limits for both subplots
contour1.set_clim(0, losses_rearro.max())
contour2.set_clim(0, losses_rearro2.max())

# Create a single colorbar on the right side
cbar_ax = fig.add_axes([1, 0.15, 0.02, 0.7])  # [x, y, width, height]
cbar = fig.colorbar(contour2, cax=cbar_ax, orientation='vertical')

# Add a title to the colorbar
cbar.set_label('Losses (W)', rotation=270, labelpad=15)

# Adjust layout for better spacing
plt.tight_layout()

# Show the plots
plt.show()

# %%
import matplotlib.pyplot as plt
import numpy as np

# Assuming you have df_th_meas DataFrame and necessary data

# Plot the first contourf
losses_rearro1 = np.reshape(df_th_meas['Loss diff'].to_numpy(), (len(np.unique(df_th_meas['D1'])), len(np.unique(df_th_meas['D2']))))
losses_rearro1 = losses_rearro1.T
contour1 = plt.contourf(np.unique(df_th_meas['D1']), np.unique(df_th_meas['D2']), losses_rearro1, 100, cmap='hot', linewidth=0)

plt.xlabel('D1')
plt.ylabel('D2')

# Create a colorbar for the first plot
cbar1 = plt.colorbar(contour1, orientation='vertical')
cbar1.set_label('Losses difference (W)', rotation=270, labelpad=15)

# Show the first plot
plt.show()




# %% [markdown]
# ## HEATMAP WITH MASK

# %%
fig, ax = plt.subplots()
#añadir protection mosfets OJO
phi_bob=phi*np.ones(df_hm.shape[0])
c=currents_L(df_hm['L'],df_hm['Vin'],Vbus,df_hm['Iin(A)'],D1_r,D2_r,phi_bob,T)
interv=intervals(D1_r,D2_r,phi_bob,T)
intervdf=pd.DataFrame(interv[-1])
df_hm_los1=df_hm.copy()
Irms_M1,Irms_M2,Irms_M3,Irms_M4=Irms_Mi(D1_r,D2_r,phi_bob,fsw,interv[2],interv[5],interv[8],interv[11],interv[1],c[0],c[1],c[2],c[3])
aux_M1_1,Temp_M1_1=Losses_M1(D1_r,Vin_f(Vbus,D1_r,D2_r),c[3],c[1],fsw,M1_index,Vdrvr,Idrvr,Irms_M1,max_age,Rjc_M1,Rca_M1,Tini,Tjmax_M1,df_hm_los1)
aux_M2_1,Temp_M2_1=Losses_M2(D1_r,Vin_f(Vbus,D1_r,D2_r),c[1],c[3],fsw,M2_index,Vdrvr,Idrvr,Irms_M2,max_age,Rjc_M2,Rca_M2,Tini,Tjmax_M2,df_hm_los1)
aux_M3_1,Temp_M3_1=Losses_M3(D2_r,Vbus,c[2],c[0],fsw,M3_index,Vdrvr,Idrvr,Irms_M3,max_age,Rjc_M3,Rca_M3,Tini,Tjmax_M3,df_hm_los1)
aux_M4_1,Temp_M4_1=Losses_M4(D2_r,Vbus,c[0],c[2],fsw,M4_index,Vdrvr,Idrvr,Irms_M4,max_age,Rjc_M4,Rca_M4,Tini,Tjmax_M4,df_hm_los1)

Iac_ind=c[-1].max(axis=1)-c[-1].min(axis=1) #Validated with other implementation. Works

df_hm_los1['Iac_ind']=Iac_ind
Iavg_ind=Iind_avg(c[3],c[0],c[1],c[2],interv[2]-interv[1],interv[5]-interv[4],interv[8]-interv[7],interv[11]-interv[10])
Irms1=Irms_sec(1/T,interv[2],interv[1],c[0],c[3])
Irms2=Irms_sec(1/T,interv[5],interv[4],c[1],c[0])
Irms3=Irms_sec(1/T,interv[8],interv[7],c[2],c[1])
Irms4=Irms_sec(1/T,interv[11],interv[10],c[3],c[2])
Irms_ind=np.sqrt(Irms1**2+Irms2**2+Irms3**2+Irms4**2)

df_hm_los1['Irms_ind']=Irms_ind
df_hm_los1['Ipk_ind']=c[-1].max(axis=1)
df_hm_los1['Iavg']=(c[-1].max(axis=1)+c[-1].min(axis=1))/2
aux_L_1=Losses_ind(c_ind,x_ind,y_ind,Vol_ind,fsw,L_ind,Iac_ind,N_ind,Ae_ind,Irms_ind,Awire_ind,Lwire_ind,rho)
aux1=aux_M1_1+aux_M2_1+aux_M3_1+aux_M4_1+aux_L_1-aux
losses_1=aux_M1_1+aux_M2_1+aux_M3_1+aux_M4_1+aux_L_1
loss_difference=aux1
df_hm_los1['Losses']=df_hm_los['Losses']
df_hm_los1['Losses_aged']=losses_1
df_hm_los1['Loss_diff']=loss_difference
df_hm_los1['aux_M1']=aux_M1_1
df_hm_los1['aux_M2']=aux_M2_1
df_hm_los1['aux_M3']=aux_M3_1
df_hm_los1['aux_M4']=aux_M4_1
df_hm_los1['aux_L']=aux_L_1
df_hm_los1['Temp_M1']=Temp_M1_1
df_hm_los1['Temp_M2']=Temp_M2_1
df_hm_los1['Temp_M3']=Temp_M3_1
df_hm_los1['Temp_M4']=Temp_M4_1
df_hm_los1['Small losses']=df_hm_los1['Loss_diff'].mask((df_hm_los1['Loss_diff']<0.5)|(df_hm_los1['Temp_M1']>70)|(df_hm_los1['Temp_M2']>70)|(df_hm_los1['Temp_M3']>70)|(df_hm_los1['Temp_M4']>70) )

heat=df_hm_los1.loc[:,['D1','D2','Loss_diff']]
heat=heat.set_index('D2')
heat=heat.pivot(columns='D1')
heat_val=heat.values
D1_x=heat.columns.get_level_values(1)
D2_y=heat.index.values

impossible_values=df_hm_los1.loc[:,['D1','D2','Small losses']]
impossible_values=impossible_values.set_index('D2')
impossible_values=impossible_values.pivot(columns='D1')
impossible_values_val=impossible_values.values
D1_x=impossible_values.columns.get_level_values(1)
D2_y=impossible_values.index.values

#Plot the surface.
# gr=ax.contourf(D1_x, D2_y, heat_val, 50, cmap='coolwarm')#,vmin=6.9,vmax=30)  #'RdYlBu'
# fig.colorbar(gr)

gr_imp=ax.contourf(D1_x, D2_y, impossible_values_val, 50, cmap='hot',corner_mask=True)#,vmin=6.9,vmax=30)  #'RdYlBu'
# gr_imp=ax.contour(D1_x, D2_y, impossible_values_val, 50, colors='k')
fig.colorbar(gr_imp)
# cbar=fig.colorbar(gr_imp)
# ticklabs = cbar.ax.get_yticklabels()
# cbar.ax.set_yticklabels(ticklabs, fontsize=20)

# Las líneas negras representan las parejas de puntos que dan la mínima y máxima ganancia necesarias.
ax.plot(D1_r,D2_f(Vsa_max,Vbus,D1_r),'k')
ax.plot(D1_r,D2_f(Vsa_min,Vbus,D1_r),'k')

# Las líneas blancas representan la ganancia máxima que podría alcanzarse si se imponen límites del 0.95 a D1 y D2.
# Es decir, se puede fijar un D de 1 o menor de 0.95,pero nunca entre ambos valores por provlemas con retrassos, tiempos
#muertos, etc.

ax.plot(D1_r,D2_f(Vbus,Lim_sup_D1*Vbus,D1_r),'w')
ax.plot(D1_f(Vbus*(1-Lim_inf_D2),Vbus,D2_r),D2_r,'w')


#ax.plot(D1_r,D2_r,'g')
plt.xlim(D1_min_sh,D1_max_sh)
plt.ylim(D2_min_sh,D2_max_sh)
# plt.xticks(fontsize = 20)
# plt.yticks(fontsize = 20)

# f = plt.figure()
# f.set_figwidth(10)
# f.set_figheight(5)

plt.show()



# %%


# %%




# %%
heat

# %%
impossible_values




# %%


# %%
gain_min=Vbus/Vsa_max
gain_max=Vbus/Vsa_min
gain_vert=D1_max/(1-0) #Máxima ganancia que se puede conseguir operando con D2=0 (modo Buck)
gain_horz=1/(1-D2_min) #Mínima ganancia que se puede conseguir operando con D1=1 (modo Boost)


# %% [markdown]
# NEW DEFINITION OF VARIABLES

# %%
elements=30

# %%
L_ind=83.5e-6
i=0
losses=[0 for i in range(elements)]
gainr=np.linspace(gain_min,gain_max,elements)
phi_bob=0.4*0
phi_bob=phi_bob*np.ones(elements)
D1x=np.where(gainr<gain_vert,gainr,np.where(gainr<gain_horz,D1_max,1))
D2x=np.where(gainr<gain_vert,0,1-D1x/gainr)
Vbus=Vbus
Vin=Vbus/gainr
fsw=180000
T=1/fsw
#We create a dataframe and add the different values
df=pd.DataFrame(columns=['Gain','D1','D2'])
df['Gain']=gainr
df['D1']=D1x
df['D2']=D2x
df['phi']=phi_bob*np.ones(df.shape[0])
df['T']=T*np.ones(df.shape[0])
df['fsw']=fsw*np.ones(df.shape[0])
df['L']=L_ind*np.ones(df.shape[0])
df['Vin']=Vbus*np.ones(df.shape[0])/gainr
df['Vbus']=Vbus*np.ones(df.shape[0])
df['P']=P*np.ones(df.shape[0])
df['Iin(A)']=P/(Vbus/gainr)
# df_Gain=df.set_index('Gain')
Iin=df['Iin(A)'].to_numpy()


#change to multiindex with D1 and D2 because this way it would be easier to develop the heatmap
indexmap=pd.MultiIndex.from_arrays([df['D1'],df['D2']],names=('D1','D2'))
df1mi=df.copy()#If not copied the original can be modified
df1mi=df1mi.drop(columns=['D1','D2'])
df1mi.index=indexmap
pd.Series(indexmap.get_level_values(0)) #If I use d1 and d2 as multiindex, it is to get those indexes to use them in t1

c=currents_L(df['L'],df['Vin'],Vbus,df['Iin(A)'],D1x,D2x,phi_bob,T)

interv=intervals(D1x,D2x,phi_bob,T)
intervdf=pd.DataFrame(interv[-1])
Irms_M1,Irms_M2,Irms_M3,Irms_M4=Irms_Mi(D1x,D2x,phi_bob,fsw,interv[2],interv[5],interv[8],interv[11],interv[1],c[0],c[1],c[2],c[3])
aux_M1=Losses_M1(D1x,Vbus/gainr,c[3],c[1],fsw,M1_index,Vdrvr,Idrvr,Irms_M1,aging,Rjc_M1,Rca_M1,Tini,Tjmax_M1,df)
aux_M2=Losses_M2(D1x,Vbus/gainr,c[1],c[3],fsw,M2_index,Vdrvr,Idrvr,Irms_M2,aging,Rjc_M2,Rca_M2,Tini,Tjmax_M2,df)
aux_M3=Losses_M3(D2x,Vbus,c[2],c[0],fsw,M3_index,Vdrvr,Idrvr,Irms_M3,aging,Rjc_M3,Rca_M3,Tini,Tjmax_M3,df)
aux_M4=Losses_M4(D2x,Vbus,c[0],c[2],fsw,M4_index,Vdrvr,Idrvr,Irms_M4,aging,Rjc_M4,Rca_M4,Tini,Tjmax_M4,df)
Iac_ind=c[-1].max(axis=1)-c[-1].min(axis=1) #Validated with other implementation. Works
df_los=df.copy()
df_los['Iac_ind']=Iac_ind
Iavg_ind=Iind_avg(c[3],c[0],c[1],c[2],interv[2]-interv[1],interv[5]-interv[4],interv[8]-interv[7],interv[11]-interv[10])
Irms1=Irms_sec(1/T,interv[2],interv[1],c[0],c[3])
Irms2=Irms_sec(1/T,interv[5],interv[4],c[1],c[0])
Irms3=Irms_sec(1/T,interv[8],interv[7],c[2],c[1])
Irms4=Irms_sec(1/T,interv[11],interv[10],c[3],c[2])
Irms_ind=np.sqrt(Irms1**2+Irms2**2+Irms3**2+Irms4**2)
df_los['Irms_ind']=Irms_ind
df_los['Ipk_ind']=c[-1].max(axis=1)
df_los['Iavg']=(c[-1].max(axis=1)+c[-1].min(axis=1))/2
aux_L=Losses_ind(c_ind,x_ind,y_ind,Vol_ind,fsw,L_ind,Iac_ind,N_ind,Ae_ind,Irms_ind,Awire_ind,Lwire_ind,rho)
aux=aux_M1+aux_M2+aux_M3+aux_M4+aux_L
losses=aux
df_los['Losses']=losses
df_los['aux_M1']=aux_M1
df_los['aux_M2']=aux_M2
df_los['aux_M3']=aux_M3
df_los['aux_M4']=aux_M4
df_los['aux_L']=aux_L
# interv[-1]
# c[-1]
# interv
plt.plot(interv[-1].loc[0,'t1_0':'t1_49'], Ip1(interv[-1].iloc[0,0:50],L_ind,df_los.loc[0,'Vin'],Vbus,c[-1].loc[0,'Ip4_end'],phi_bob[0],D1x[0],D2x[0]),'blue')
plt.plot(interv[-1].iloc[0,52:102],Ip2(interv[-1].iloc[0,52:102]-interv[-1].iloc[0,51],L_ind,df_los.loc[0,'Vin'],Vbus,c[-1].loc[0,'Ip1_end'],phi_bob[0],D1x[0],D2x[0]),'red')
plt.plot(interv[-1].iloc[0,104:154],Ip3(interv[-1].iloc[0,104:154]-interv[-1].iloc[0,103],L_ind,df_los.loc[0,'Vin'],Vbus,c[-1].loc[0,'Ip2_end'],phi_bob[0],D1x[0],D2x[0]),'orange')
plt.plot(interv[-1].iloc[0,156:206],Ip4(interv[-1].iloc[0,156:206]-interv[-1].iloc[0,155],L_ind,df_los.loc[0,'Vin'],Vbus,c[-1].loc[0,'Ip3_end'],phi_bob[0],D1x[0],D2x[0]),'green')
plt.show()
df_los.head() 




# %%


# %%


# %%


# %%


# %%


# %%



# %% [markdown]
# ANNEX: RMS VALUE VALIDATION

# %%
# c[-1].min(axis=1)
voboost=100
ioutboost=5
vsboost=80
dutyboost=(voboost-vsboost)/voboost
periodoboost=1/180000
inducboost=83.5e-6

Iavgboost=voboost*ioutboost/vsboost
Iripboost=(vsboost*dutyboost)*periodoboost/inducboost
iminboost=Iavgboost-Iripboost/2
Irmsboost=(np.sqrt(iminboost*iminboost+iminboost*Iripboost+Iripboost*Iripboost/3))
Irmsboost2=np.sqrt(Iavgboost*Iavgboost+Iripboost*Iripboost/12)
print('Iavg',Iavgboost)
print('Irms',Irmsboost)
print('irip',Iripboost)
print('Ipk', iminboost+Iripboost)
print('Irmsboost',Irmsboost2)

# %%
vobuck=100
ioutbuck=5
vsbuck=150
dutybuc=(vobuck/vsbuck)
periodobuck=1/180000
inducbuck=83.5e-6

ioutbuck=500/vobuck
Iavgbuck=ioutbuck
Iripbuck=(vsbuck-vobuck)*periodobuck*dutybuc/inducbuck
iminbuck=Iavgbuck-Iripbuck/2
Irmsbuck=(np.sqrt(iminbuck*iminbuck+iminbuck*Iripbuck+Iripbuck*Iripbuck/3))
Irmsbuck2=np.sqrt(Iavgbuck*Iavgbuck+Iripbuck*Iripbuck/12)
print('Iavg',Iavgbuck)
print('Irms',Irmsbuck)
print('irip',Iripbuck)
print('Ipk', iminbuck+Iripbuck)
print('Irmsbuck',Irmsbuck2)

# %% [markdown]
# MATRIZ-MODEL

# %%


# %%
from itertools import product
values = list(product(range(aging_steps), range(aging_steps), range(aging_steps),range(aging_steps),gainr)) #to get all possible combinations
dfrds=pd.DataFrame(values)
dfrds.columns=['iRds1','iRds2','iRds3','iRds4','Gain']
df_aged=pd.merge(dfrds,df)
c=currents_L(df_aged['L'],df_aged['Vin'],df_aged['Vbus'],df_aged['Iin(A)'],df_aged['D1'],df_aged['D2'],df_aged['phi'],df_aged['T'])
interv=intervals(df_aged['D1'],df_aged['D2'],df_aged['phi'],df_aged['T'])
Irms_M1,Irms_M2,Irms_M3,Irms_M4=Irms_Mi(df_aged['D1'],df_aged['D2'],df_aged['phi'],1/df_aged['T'],interv[2],interv[5],interv[8],interv[11],interv[1],c[0],c[1],c[2],c[3])
g=df_aged['iRds1']
aux_M1=Losses_M1(df_aged['D1'],df_aged['Vin'],c[3],c[1],fsw,M1_index,Vdrvr,Idrvr,Irms_M1,1+g/aging_steps*aging_max)
h=df_aged['iRds2']
aux_M2=Losses_M2(df_aged['D1'],df_aged['Vin'],c[1],c[3],fsw,M2_index,Vdrvr,Idrvr,Irms_M2,1+h/aging_steps*aging_max)
i=df_aged['iRds3']
aux_M3=Losses_M3(df_aged['D2'],df_aged['Vbus'],c[2],c[0],1/df_aged['T'],M3_index,Vdrvr,Idrvr,Irms_M3,1+i/aging_steps*aging_max)
j=df_aged['iRds4']
aux_M4=Losses_M4(df_aged['D2'],df_aged['Vbus'],c[0],c[2],1/df_aged['T'],M4_index,Vdrvr,Idrvr,Irms_M4,1+j/aging_steps*aging_max)
Iac_ind=c[-1].max(axis=1)-c[-1].min(axis=1) 
dfrds['Iac_ind']=Iac_ind
Irms_ind=Iind_avg(c[3],c[0],c[1],c[2],interv[2]-interv[1],interv[5]-interv[4],interv[8]-interv[7],interv[11]-interv[10])
df_aged['Irms_ind']=Irms_ind
df_aged['Ipk_ind']=c[-1].max(axis=1)

aux_L=Losses_ind(c_ind,x_ind,y_ind,Vol_ind,fsw,L_ind,Iac_ind,N_ind,Ae_ind,Irms_ind,Awire_ind,Lwire_ind,rho)
losses_aged=aux_M1+aux_M2+aux_M3+aux_M4+aux_L
df_aged['Losses aged']=losses_aged #Este es el aging data del codgio de manu

fig, ax = plt.subplots()

ax.plot(df_aged.loc[(df_aged['iRds1']==0)&(df_aged['iRds2']==0)&(df_aged['iRds3']==0)&(df_aged['iRds4']==0),'Gain'],df_aged.loc[(df_aged['iRds1']==0)&(df_aged['iRds2']==0)&(df_aged['iRds3']==0)&(df_aged['iRds4']==0),'Losses aged'],'y')
ax.plot(df_aged.loc[(df_aged['iRds1']==1)&(df_aged['iRds2']==0)&(df_aged['iRds3']==0)&(df_aged['iRds4']==0),'Gain'],df_aged.loc[(df_aged['iRds1']==1)&(df_aged['iRds2']==1)&(df_aged['iRds3']==1)&(df_aged['iRds4']==1),'Losses aged'],'k')
ax.plot(df_aged.loc[(df_aged['iRds1']==2)&(df_aged['iRds2']==0)&(df_aged['iRds3']==0)&(df_aged['iRds4']==0),'Gain'],df_aged.loc[(df_aged['iRds1']==2)&(df_aged['iRds2']==2)&(df_aged['iRds3']==2)&(df_aged['iRds4']==2),'Losses aged'],'g')
ax.plot(df_aged.loc[(df_aged['iRds1']==3)&(df_aged['iRds2']==3)&(df_aged['iRds3']==0)&(df_aged['iRds4']==0),'Gain'],df_aged.loc[(df_aged['iRds1']==3)&(df_aged['iRds2']==3)&(df_aged['iRds3']==3)&(df_aged['iRds4']==3),'Losses aged'],'r')
# Ya tengo una base de datos con todos los posibles valores de pérdidas en función del envejecimiento de cada MOSFET.
# aux1[age_MOSFET1][age_MOSFET2][age_MOSFET3][Age_MOSFET4][gain]
# Se tienen n edades de envejecimiento para cada MOSFET

#repeat what was done with the other 

def aging_estimation(coef_1_ini,coef_2_ini,coef_3_ini,coef_4_ini,error,df_losses):
    muestra_r=[0 for sample in range(df_losses.shape[0])]
    muestra=[0 for sample in range(df_losses.shape[0])]
    losses_not_damaged=[]
    df_not_damaged=pd.DataFrame()
    # Se toma una muestra de la diferencia entre un convertidor sano y uno envejecido según los coeficientes de arriba (muestra)
    # y la misma muestra con errores aleatorios (muestra_r) debido a medida con un máximo de (error%)
    for gain in df_losses['Gain'].unique():#variación de Gain
        #print("gain is:", gain)
        losses_not_damaged=np.append(losses_not_damaged,df_losses.loc[(df_losses['iRds1']==0)&(df_losses['iRds2']==0)&(df_losses['iRds3']==0)&(df_losses['iRds4']==0)& (df_losses['Gain']==gain),'Losses aged'])
    df_not_damaged['Gain']=df_losses['Gain'].unique()
    df_not_damaged['Losses not damaged']=losses_not_damaged
    df_losses=pd.merge(df_losses,df_not_damaged)
    df_losses['Diff_losses']=df_losses.loc[:,'Losses aged']-df_losses['Losses not damaged']
    df_losses['Diff_losses noise']=df_losses['Diff_losses']*(1+random.randint(-error, error)/1000)
    
    error_ind=((df_losses['Diff_losses noise']-df_losses['Diff_losses'])**2)
    df_losses['MSE_single']=error_ind 

    # Se busca por mínimos cuadrados la combinación de coefficentes que más aproximan el modelo a la muestra con errores de medición a
    # Es importante resaltar  Mque no se ha hecho con el modelo matemático, sino explorando la matrix de reusltados que se han sacado
    # con el modelo. Básicamente es lo mismo, pero sería bueno explorar si a partir de las ecuaciones se puede ganar algo.

    values_dummy=pd.DataFrame()

    r_valores = np.array(list(product(df_losses['iRds1'].unique(),df_losses['iRds2'].unique(), df_losses['iRds3'].unique(),df_losses['iRds4'].unique())))
    MSE=np.zeros(r_valores.shape[0]) 
    #Me quedo con todas las ocurrencias 
    for i in range(r_valores.shape[0]):
        sel_aging=df_aged.loc[(df_losses['iRds1']==r_valores[i][0])&(df_losses['iRds2']==r_valores[i][1])&(df_losses['iRds3']==r_valores[i][2])&(df_losses['iRds4']==r_valores[i][3])]

        suma_MSE=sel_aging['MSE_single'].sum()
        MSE[i]=(suma_MSE)#/(r_valores.shape[0])
    values_dummy['iRds1','iRds2','iRds3','iRds4','MSE']=r_valores[:,0],r_valores[:,1],r_valores[:,2],r_valores[:,3],MSE

    df_losses=pd.merge(df_losses,values_dummy)
    coef_1_out,coef_2_out,coef_3_out,coef_4_out=df_losses[df_losses['MSE']==df_losses['MSE'].min(),['iRds1','iRds2','iRds3','iRds4']] 
    return coef_1_out,coef_2_out,coef_3_out,coef_4_out,muestra,muestra_r


# %%
# MSE=np.zeros(r_valores.shape[0])
# r_valores=np.array(r_valores)
# for i in range(r_valores.shape[0]):
#     sel_aging=df_aged.loc[(df_aged['iRds1']==r_valores[i][0])&(df_aged['iRds2']==r_valores[i][1])&(df_aged['iRds3']==r_valores[i][2])&(df_aged['iRds4']==r_valores[i][3])]
#     MSE[i]=sel_aging['Losses aged'].sum() 
# MSE

# %%
# sel_aging

# %%
coef_1_ini=5
coef_2_ini=2
coef_3_ini=2
coef_4_ini=5

# muestra_r=[0 for gainr in range(elements)]
# muestra=[0 for gainr in range(elements)]
df_losses=df_aged.copy()
error=error_sense
# df_aged2=aging_estimation(coef_1_ini,coef_2_ini,coef_3_ini,coef_4_ini,error_sense,df_aged)
# df_aged2
# coef_1,coef_2,coef_3,coef_4,muestra,muestra_r=aging_estimation2(df_aged2)


muestra_r=[0 for sample in range(df_losses.shape[0])]
muestra=[0 for sample in range(df_losses.shape[0])]
losses_not_damaged=[]
df_not_damaged=pd.DataFrame()
    # Se toma una muestra de la diferencia entre un convertidor sano y uno envejecido según los coeficientes de arriba (muestra)
    # y la misma muestra con errores aleatorios (muestra_r) debido a medida con un máximo de (error%)
for gain in df_losses['Gain'].unique():#variación de Gain
       #print("gain is:", gain)
    losses_not_damaged=np.append(losses_not_damaged,df_losses.loc[(df_losses['iRds1']==0)&(df_losses['iRds2']==0)&(df_losses['iRds3']==0)&(df_losses['iRds4']==0)& (df_losses['Gain']==gain),'Losses aged'])
df_not_damaged['Gain']=df_losses['Gain'].unique()
df_not_damaged['Losses not damaged']=losses_not_damaged
df_losses=pd.merge(df_losses,df_not_damaged)
df_losses['Diff_losses']=df_losses.loc[:,'Losses aged']-df_losses['Losses not damaged']
df_losses['Diff_losses noise']=df_losses['Diff_losses']*(1+random.randint(-error, error)/1000)
    
error_ind=((df_losses['Diff_losses noise']-df_losses['Diff_losses'])**2)
df_losses['MSE_single']=error_ind 

    # Se busca por mínimos cuadrados la combinación de coefficentes que más aproximan el modelo a la muestra con errores de medición a
    # Es importante resaltar  Mque no se ha hecho con el modelo matemático, sino explorando la matrix de reusltados que se han sacado
    # con el modelo. Básicamente es lo mismo, pero sería bueno explorar si a partir de las ecuaciones se puede ganar algo.

values_dummy=pd.DataFrame()

r_valores = np.array(list(product(df_losses['iRds1'].unique(),df_losses['iRds2'].unique(), df_losses['iRds3'].unique(),df_losses['iRds4'].unique())))
MSE=np.zeros(r_valores.shape[0]) 
for i in range(r_valores.shape[0]):
    sel_aging=df_aged.loc[(df_losses['iRds1']==r_valores[i][0])&(df_losses['iRds2']==r_valores[i][1])&(df_losses['iRds3']==r_valores[i][2])&(df_losses['iRds4']==r_valores[i][3])]

    suma_MSE=sel_aging['MSE_single'].sum()
    MSE[i]=(suma_MSE)#/(r_valores.shape[0])
values_dummy['iRds1']=r_valores[:,0]
values_dummy['iRds2']=r_valores[:,1]
values_dummy['iRds3']=r_valores[:,2]
values_dummy['iRds4']=r_valores[:,3]
values_dummy['MSE']=MSE

df_losses=pd.merge(df_losses,values_dummy)
df_losses=df_losses.drop(df_losses[df_losses['MSE']==0].index)
coef_1_out,coef_2_out,coef_3_out,coef_4_out=df_losses.loc[df_losses['MSE']==df_losses['MSE'].min(),['iRds1','iRds2','iRds3','iRds4']] 
df_losses.loc[df_losses['MSE']==df_losses['MSE'].min()]

# %%
df_losses.loc[(df_losses['MSE']==df_losses['MSE'].min())&(df_losses['MSE']!=0),['Gain','iRds1','iRds2','iRds3','iRds4']] 



# %%


# %%




# %% [markdown]
# Generacion base de datos para un perfil de ganancia (ya hecho)

# %% [markdown]
# 

# %% [markdown]
# 

# %% [markdown]
# 

# %% [markdown]
# 

# %%



# %%


# %%



