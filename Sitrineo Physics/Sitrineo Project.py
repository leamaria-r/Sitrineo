#%% Projet Sitrineo

# =============================================================================
# The objective of this code is firstly to calculate the activity of our radioactive source and display its activity,
# secondly, to simulate the evolution of the electrons produced (Energy, Scattering) as well as obtain the physical limitations of the experiment such as the impact zone on the different planes, the cuts, which could be useful
# for the reconstruction of traces in the 'TAF' software
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
#%% Constants Source
N_Avo=6.02214076e23
#Strontium 90
Activ_Sr=5.11e12 #Bq/g
Activ_echantil_a_t_0=37e6 #Bq
m_Sr=Activ_echantil_a_t_0/Activ_Sr #g
M_mol_Sr=89.90793 #g/mol
tau_Sr=907921440/np.log(2) #en s
N_Sr0=m_Sr/M_mol_Sr*N_Avo #Nbr Atome
#Yttrium 90
Activ_Y=10.22e15 #Bq/g
M_mol_Y=89.907142*1.660538921e-24*N_Avo #g/mol
tau_Y=233280/np.log(2) #en s

#%% Evolution Function

def N_Sr(x): #Strontium 90 population evolution through time
    return N_Sr0*np.exp(-x/tau_Sr)

def dN_Sr(x): #Strontium 90 variation evolution through time
    return -N_Sr0/tau_Sr*np.exp(-x/tau_Sr)

def N_Y(x): #Yttrium 90 population evolution through time
    return N_Sr0*(np.exp(-x/tau_Sr)-np.exp(-x/tau_Y))/(tau_Sr/tau_Y-1)

def dN_Y(x): #Yttrium 90 variation evolution through time
    return -N_Sr0*(np.exp(-x/tau_Sr)/tau_Sr-np.exp(-x/tau_Y)/tau_Y)/(tau_Sr/tau_Y-1)

def N_Zr(x): #Zirconium 90 population evolution through time
    return N_Sr0*(-tau_Sr*np.exp(-x/tau_Sr)+tau_Y*np.exp(-x/tau_Y))/(tau_Sr-tau_Y)+N_Sr0
    

t=np.linspace(0, 1e10,num=1000000)
tab_Sr=N_Sr(t)
tab_Y=N_Y(t)
tab_Zr=N_Zr(t)

#%% Time
#Source Actuelle : t0=05/09/2017
start_date=dt.date(2017,9,5)
end_date=dt.date.today()
now=(end_date-start_date).days #days
print("Day difference between source reception and today: {}".format(now))

t_now=now*24*3600 #seconds

#Points to draw a line on the time plot to see our current position in time
barre=np.ones(1000)*t_now
k=np.linspace(0,N_Sr0,1000)

#%% Plot Evolution Source

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = False
fig, ax1 = plt.subplots()
l1,l2,l3,l4=ax1.plot(t,tab_Sr,'r',t,tab_Y,'g',t,tab_Zr,'b',barre,k,'y')
ax1.set_xscale('log')
plt.title("Source's composition evolution")
plt.xlabel('Time (s)')
plt.ylabel('Abondancy (Atoms)')
l, b, h, w = .18, .25, .3, .2
ax2 = fig.add_axes([l, b, w, h])
ax2.plot(t,tab_Y,'g')
ax2.plot(barre,k,'y')
ax2.set_ylim((0,0.0005*N_Sr0))
ax2.set_xlim((1e6,1e11))
ax2.set_xscale('log')
ax1.legend((l1,l2,l3,l4),('Strontium 90','Yttrium 90','Zirconium 90','Currently'),loc='upper left')

#%% Plot Activité Source

l=np.linspace(0,7.5e7,1000)
plt.figure()
plt.plot(t,-dN_Sr(t),'b',label='Activity Strontium')
plt.plot(t,-dN_Sr(t)-dN_Y(t),'r',label='Activity Yttrium')
plt.plot(t,-2*dN_Sr(t)-dN_Y(t),'k',label='Activity Source')
plt.plot(barre,l,'y',label='Currently')
plt.xscale('log')
plt.xlabel('Time (s)')
plt.ylabel('Activity (Bq)')
plt.title('Evolution of the activity of the source as a function of time')
plt.legend(loc='lower left')
#%% Activité Source at time y

#Our Source currently:
y=t_now
#print(N_Y(t_now))
Activ_Sr_now=np.round(dN_Sr(y),0) #Bq
Activ_Y_now=np.round(dN_Y(y)+dN_Sr(y),0) #Bq
Activ_Source=np.abs(Activ_Sr_now+Activ_Y_now)
ratio_Y=np.abs(Activ_Y_now/Activ_Source)
ratio_Sr=np.abs(Activ_Sr_now/Activ_Source)

print('Currently there are {} decays per second, {}% from Yttrium and {} from Strontium'.format(Activ_Source,ratio_Y,ratio_Sr))
#%% Class Sitrineo

# =============================================================================
# The Sitrineo class was thought for two mains aspects:
#     1-Calculate the number of particule going through a collimator with a certain radius at a certain distance from the source
#     2-Calculate, for a certain kinetic energy T of our electron, the amount of energy lost in the tracking device and the deviation induced
# =============================================================================

class Sitrineo:
    
    def __init__(self,radius_collimateur,distance):
        self.radius_collimateur = radius_collimateur #Radius of Collimator 
        self.distance = distance #Distance between source and exit point of the collimator
        
    def find_nearest(self,array, values): #Function used in stopping power function to get the stopping power value depending on the T value 
        indices = np.abs(np.subtract.outer(array, values)).argmin(0)
        return indices
        
    def increase_point(self,x,n=1000): #Function used in stopping power function to increase the number of points from the original file #A little archaic but if we have an x ​​and a y of the same length then we will have the same increase in the number of points
            y=np.zeros(2*len(x)) #New array of zeros two times longer than the original array
            for i in range (len(x)-1):
                y[2*i]=x[i] #The even indices of our new table with the old values
                y[2*i+1]=(x[i]+x[i+1])/2 #Les indices impaires avec la moyenne entre deux anciennes valeurs
            y[-2]=x[-1] #Our last two values ​​are 0, correction to have a well-defined table
            y=y[:-1]
            if len(y)<n:
                tab=self.increase_point(y,n)
                return tab
            else:
                return y
    
    def Area_circle(self,r):
        return np.pi*r**2
    
    def Radius_of_Sphere(self,round=10): #Function that calculate the radius of the sphere that represents 4pi rad
        return np.round(np.sqrt(self.radius_collimateur**2+self.distance**2),round)
    
    def Steradian(self): #Function that calculate the solid angle due to collimator
        unit=self.Radius_of_Sphere()**2 #1 steradian = unit surface
        return self.Area_circle(self.radius_collimateur)/unit
    
    def Activite(self,source): #This function calculates the number of particule going through the collimator each second
        if source<0:
           source=abs(source)
        res=np.round(source*self.Steradian()/(4*np.pi),0)
        return res
    
    def Momentum(self,T,m):
        E=T+m
        p=np.sqrt(E**2-m**2)
        return p,E
    
    def Mult_Scat(self,T,matter,x,z=-1,m=511e3): #This function calculates the angle of deviation in rad of due to particule going through material depending on its T
        if matter=='air':
            l0=36.62 #g/cm²
            rho=1.204e-3 #g/cm3
            X0=l0/rho
        elif matter=='silicium':
            l0=21.82 #g/cm²
            rho=2.33 #g/cm3
            X0=l0/rho
        elif matter=='aluminium':
            l0=24.01 #g/cm²
            rho=2.698 #g/cm3
            X0=l0/rho
        
        p,E=self.Momentum(T,m)
        beta=p/E
        theta=13.6e6*abs(z)*np.sqrt(x/X0)*(1+0.038*np.log(x/X0))/(beta*p) #in °
        theta=theta*np.pi/180 #in rad
        return theta
    
    def Stopping_power(self,T,matter): #This function is used for the Attenuation function, returns the characteristics of the material encountered by the particule in function of its T
        if matter=='air':
            rho_air=1.204e-3 #g/cm3
            tab=np.loadtxt('Stopping Power Air.dat',delimiter='|',unpack=True)*1e6 #eV cm²/g per eV
            E_e=self.increase_point(tab[0])
            tab_stop_air=self.increase_point(tab[1])
            idx=self.find_nearest(E_e,T)
            stop_air=tab_stop_air[idx]
            return stop_air,rho_air
        
        elif matter=='silicium':
            rho_Si=2.33 #g/cm3
            tab=np.loadtxt('Stopping Power Silicium.dat',delimiter='|',unpack=True)*1e6 #eV cm²/g per eV
            E_e=self.increase_point(tab[0])
            tab_stop_Si=self.increase_point(tab[1])
            idx=self.find_nearest(E_e,T)
            stop_Si=tab_stop_Si[idx]
            return stop_Si,rho_Si
        
        elif matter=='aluminium':
            rho_Al= 2.698 #g/cm3
            tab=np.loadtxt('Stopping Power Aluminium.dat',delimiter='|',unpack=True)*1e6 #eV cm²/g per eV
            E_e=self.increase_point(tab[0])
            tab_stop_Al=self.increase_point(tab[1])
            idx=self.find_nearest(E_e,T)
            stop_Al=tab_stop_Al[idx]
            return stop_Al,rho_Al
        
        else:
            return "Not coded already"
        
    def Attenuation(self,T,matter,length): #This function calculates the energy deposited/lost by the particule due to the interactions with material
        stop_power,rho_matter=self.Stopping_power(T, matter)
        dE_dX=stop_power*rho_matter
        E_lost=dE_dX*length
        return E_lost
    
    def Physic(self,T,dist_plane): #This function aims to reproduce what is happening to our original electron with kinetic energy T at the source and gives us : -The angle of scattering due to each material, -The energy lost due to each material, -The kinetic energy remaining after the propagation through Sitrineo,-the gap due to scattering (Not good because it's more difficult than the sum of all deviations)
        #Could be coded more properly with loops but lack of time to make it work
        dist_source_colli=4.1500 #4.15cm = distance between collimator source and output
        dist_colli_plane=0.6 #0.6cm = distance between collimator and plane 1
        dist=dist_source_colli+dist_colli_plane
        z_plane=dist_plane+dist #cm
        
        x_Al=20e-4 #cm, thickness of our Aluminium layer, covering the detector against the photons
        x_Si=50e-4 #cm, thickness of our Silicium detector, MIMOSA 28
        T_lost_air=[]
        T_lost_Al=[]
        T_lost_Si=[]
        theta_air=[]
        theta_Al=[]
        theta_Si=[]
        
        #Between Source and Plane 1
        
        attenuation_air_0=self.Attenuation(T, 'air',z_plane[0])
        theta_air_0=self.Mult_Scat(T, 'air', z_plane[0])
        T=T-attenuation_air_0
        
        T_lost_air.append(attenuation_air_0)
        theta_air.append(theta_air_0)
        decalage=(z_plane[3])*np.tan(theta_air_0)
        
        #Plane 1
        #Layer Aluminium
        attenuation_Al_1=self.Attenuation(T, 'aluminium',x_Al)
        theta_Al_1=self.Mult_Scat(T, 'aluminium', x_Al)
        T=T-attenuation_Al_1
        
        T_lost_Al.append(attenuation_Al_1)
        theta_Al.append(theta_Al_1)
        
        #Layer Silicium
        attenuation_Si_1=self.Attenuation(T,'silicium', x_Si)
        theta_Si_1=self.Mult_Scat(T, 'silicium', x_Si)
        T=T-attenuation_Si_1
        
        T_lost_Si.append(attenuation_Si_1)
        theta_Si.append(theta_Si_1)
        
        #Layer Air until next plane
        attenuation_air_1=self.Attenuation(T, 'air',z_plane[1]-z_plane[0])
        theta_air_1=self.Mult_Scat(T, 'air', z_plane[1]-z_plane[0])
        T=T-attenuation_air_1 
        
        T_lost_air.append(attenuation_air_1)
        theta_air.append(theta_air_1)
        decalage+=(z_plane[3]-z_plane[0])*(np.tan(theta_air_1)+np.tan(theta_Si_1)+np.tan(theta_Al_1))
        
        #Plane 2
        #Layer Silicium
        attenuation_Si_2=self.Attenuation(T,'silicium', x_Si)
        theta_Si_2=self.Mult_Scat(T, 'silicium', x_Si)
        T=T-attenuation_Si_2
        
        T_lost_Si.append(attenuation_Si_2)
        theta_Si.append(theta_Si_2)
        
        #Layer Aluminium
        attenuation_Al_2=self.Attenuation(T, 'aluminium',x_Al)
        theta_Al_2=self.Mult_Scat(T, 'aluminium', x_Al)
        T=T-attenuation_Al_2
          
        T_lost_Al.append(attenuation_Al_2)
        theta_Al.append(theta_Al_2)
        
        #Layer Air until next plane
        attenuation_air_2=self.Attenuation(T, 'air',z_plane[2]-z_plane[1])
        theta_air_2=self.Mult_Scat(T, 'air', z_plane[2]-z_plane[1])
        T=T-attenuation_air_2 
          
        T_lost_air.append(attenuation_air_2)
        theta_air.append(theta_air_2)
        decalage+=(z_plane[3]-z_plane[1])*(np.tan(theta_air_2)+np.tan(theta_Si_2)+np.tan(theta_Al_2))
        
        #Plane 3
        #Layer Aluminium
        attenuation_Al_3=self.Attenuation(T, 'aluminium',x_Al)
        theta_Al_3=self.Mult_Scat(T, 'aluminium', x_Al)
        T=T-attenuation_Al_3
         
        T_lost_Al.append(attenuation_Al_3)
        theta_Al.append(theta_Al_3)
       
        #Layer Silicium
        attenuation_Si_3=self.Attenuation(T,'silicium', x_Si)
        theta_Si_3=self.Mult_Scat(T, 'silicium', x_Si)
        T=T-attenuation_Si_3
          
        T_lost_Si.append(attenuation_Si_3)
        theta_Si.append(theta_Si_3)
        
        #Layer Air until next plane
        attenuation_air_3=self.Attenuation(T, 'air',z_plane[3]-z_plane[2])
        theta_air_3=self.Mult_Scat(T, 'air', z_plane[3]-z_plane[2])
        T=T-attenuation_air_3 
           
        T_lost_air.append(attenuation_air_3)
        theta_air.append(theta_air_3)
        decalage+=(z_plane[3]-z_plane[2])*(np.tan(theta_air_3)+np.tan(theta_Si_3)+np.tan(theta_Al_3))
        
        #Plane 4
        #Layer Silicium
        attenuation_Si_4=self.Attenuation(T,'silicium', x_Si)
        theta_Si_4=self.Mult_Scat(T, 'silicium', x_Si)
        T=T-attenuation_Si_4
        
        T_lost_Si.append(attenuation_Si_4)
        theta_Si.append(theta_Si_4)
        
        #Layer Aluminium
        attenuation_Al_4=self.Attenuation(T, 'aluminium',x_Al)
        theta_Al_4=self.Mult_Scat(T, 'aluminium', x_Al)
        T=T-attenuation_Al_4
        
        T_lost_Al.append(attenuation_Al_4)
        theta_Al.append(theta_Al_4)
        
        theta=np.array([theta_air,theta_Si,theta_Al])
        T_lost=np.array([T_lost_air,T_lost_Si,T_lost_Al])
        
        return theta,T_lost,T,decalage

#Our source:
Q_Y=2.28e6 #eV
Q_Sr=546e3 #eV
T_min_to_plane4=235e3 #eV #Below this energy, particules cannot reach the last plane due to the energy lost in the differents parts (air,aluminium,silicium) 
T_Sr=np.linspace(T_min_to_plane4,546e3,10000) 
T_Y=np.linspace(T_min_to_plane4,2.28e6,10000)
plane=np.array([0,0.1,2.7,2.8]) #Distance of planes according to z, with z=0 in the foreground
dist_source_colli=4.15 #cm
radius_colli=0.1 #cm

a=Sitrineo(radius_colli, dist_source_colli)

angle=np.round(a.Steradian(),6)
passage=a.Activite(Activ_Source)
Simu=a.Physic(1e6,plane)
print(f"The solid angle associated to our {radius_colli}cm opening window collimator at {dist_source_colli}cm is: {angle}sr \n This collimator reduce our number of electron going through to {passage} electrons")
print("An electron of 1 MeV kinetic energy at source should keep {} MeV when arriving at the last plane of the Sitrineo Tracker".format(np.round(Simu[2]/1e6,3)))
#%% Calcul GeoCut Plan 1 et 4
cut_plane1=0.11447 #cm #We can calculate that with this radius of circle, on the first plane (0.6 cm after collimator), we keep the same solid angle as previously, it means that our particule should, neglecting scattering in the air, be in this circle
b=Sitrineo(cut_plane1,4.75) 
cercle=np.round(b.Steradian(),6)
print("Our particule should remain in a circle of {} cm in the first plane".format(cut_plane1))
theta1=b.Physic(235e3,plane)[0]
decalage=b.Physic(235e3, plane)[3]
cut_plane4=np.round(decalage+cut_plane1,6)
print("On the 4th plane, our particule should be within a circle of radius {} cm".format(cut_plane4))
#%% Relation between Q and T, Q~=T when T>>m_e

u_kg= 1.6605402e-27 #kg/u
c=299792458 #m/s
eV_J=6.242e18 #eV/J
M_Y=89.9071417 #u
M_Sr=89.9077279 #u
M_e=511e3 #eV

q_max=np.round((M_Sr-M_Y)*u_kg*c**2*eV_J,1) #eV
t_max=np.round(q_max*(1-(q_max+2*M_e)/(2*M_Sr*u_kg*c**2*eV_J)),1) #eV
print(f"Q of reaction is : {q_max}eV, the associated kinetic energy is: {t_max}eV ,which represents a difference of {q_max-t_max}eV")
#%% Beta Spectrum

def beta_spectrum(T,Q,Z): #This function is using the fermi distribution to calculate the associated probability to each T value to Q_max of the decay
    #Constants
    c=1  #299792458 # m/s
    hbar=1   #1.05457182e-34 # m²kg/s
    J_eV=1   #1.6022e-19 # J/eV
    m_e=511e3 #eV/c²
    alpha=1/137
    
    E=T+m_e*c**2 #eV
    p=np.sqrt((E/c)**2-(m_e*c)**2) #eV/c
    v=p/E*c #m/s
    S=np.sqrt(1-(alpha*Z)**2)-1 
    eta=4*np.pi*alpha/(hbar*v)
    
    proba=(Q-T)**2*p**2*(2*np.pi*eta)*(1/(1-np.exp(-2*np.pi*eta)))*(T**2 +4*(T*Z*alpha)**2 -1)**(S)
    
    if (type(T)==float) or (type(T)==int):
        return proba
    else:
        norm=proba/np.sum(proba)
        return norm

T_Sr_spectre=np.linspace(1e1,546e3,10000) 
T_Y_spectre=np.linspace(1e1,2.28e6,10000) 

spectre_Sr=beta_spectrum(T_Sr_spectre, 546e3, 39)
spectre_Y=beta_spectrum(T_Y_spectre, 2.28e6,40)


E_combined=np.linspace(1e1, 546e3,12000)
spectre_Sr_combined=beta_spectrum(E_combined, 546e3, 39)
spectre_Y_combined=beta_spectrum(E_combined, 2.28e6, 40)

spectre_combined1=(spectre_Sr_combined+spectre_Y_combined)

E_2=np.linspace(546e3,2.28e6,9500)
spectre_combined2=beta_spectrum(E_2, 2.28e6, 40)

plt.figure()
plt.plot(T_Sr_spectre,spectre_Sr,'b',label='Spectre Strontium')
plt.plot(T_Y_spectre,spectre_Y,'r',label='Spectre Yttrium')
plt.plot(E_combined,spectre_combined1,'k',label='Spectre Cumulé')
plt.plot(E_2,spectre_combined2,'k')
plt.xlabel('Energie (eV)')
plt.ylabel('Probabilité')
plt.title("Probability associated with the kinetic energy of beta decay emission")
plt.legend()
