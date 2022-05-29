# Berechnung von Stoffdaten für Systeme mit feuchter Luft

import numpy as np
import CoolProp.CoolProp as CP
from scipy import optimize
import matplotlib.pyplot as plt
#  from matplotlib.lines import Line2D
#  import matplotlib.patches as patches
import itertools

#from CoolProp.HumidAirProp import HAPropsSI  #can be used to check values
''' see http://www.coolprop.org/dev/fluid_properties/HumidAir.html#table-of-inputs-outputs-to-hapropssi '''

hfg0 = 2500  # kJ/kg Verdampfungsenthalpie, 0 °C
CpL = 1.004  # kJ/(kg K) Cp air
CpW = 1.86  # kJ/(kg K) Cp water vapor
CW = 4.19  # kJ/(kg K) Cp liquid water

class Moist_air:
    '''object for moist air state'''
    def __init__(self, theta, phi, p, name): # t in °C, rel. humidity, p in bar
        self.theta = theta
        self.p = p
        self.phi = phi
        self.name = name
        self.h1x = h1x_phi(theta, p, phi)
        self.x = x(theta, p, phi)
        self.dew_T = dew_T(theta, phi)
        self.wet_bulb_T = wet_bulb_T(theta, phi, p)

    def Mollier_Diagramm(self):
        p = self.p
        marker = itertools.cycle(('o', 's', 'v', 'p', '*')) 
        def isenthalp(theta,x):
            # Erzeugung von Isenthalp-Kurven
            return 1.004*theta - 2500*x

        def isotherm(theta,x):
            # Erzeugung von Isotherm-Kurven
            return 1.004*theta + x*(2500+1.86*theta)-2500*x

        def isotherm_sat(theta,x,p):
            #Erzeugung von Isotherm-Kurven-Sat
            Psat = Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")/10**5
            x_sat = 0.622*Psat/(p-Psat)
            return 1.004*theta + x_sat*(2500+1.86*theta)+ (x - x_sat)*4.19*theta-2500*x

        def isofeucht(phi,x,p):
            # Erzeugung von Isotherm-Kurven
            Psat = x*p/(phi*(0.622 +x))
            theta = CP.PropsSI("T","P",Psat*10**5,"Q",1,"water")-273.15 #in °C
            return 1.004*theta + x*(2500+1.86*theta)-2500*x
        x = np.linspace(0.001,0.1,200)
        phi = np.linspace(0.02,1,12)
        theta_isenthalp = np.linspace(0.,350.,36)
        theta = np.linspace(0.,200.,41)
        # Nehmen den minimale Wert für die Isotherme als fcn von x
        isotherm_plot = np.minimum([isotherm(theta,x) for x in x],[isotherm_sat(theta,x,p) for x in x])
        #fig = plt.figure()
        plt.plot(x,[isenthalp(theta_isenthalp,x) for x in x],'g--')
        plt.plot(x,[isofeucht(phi,x,p) for x in x],'b--')
        plt.plot(x, isotherm_plot,'r')
        plt.ylim(0,100)
        plt.xlim(0,0.1)
        plt.xlabel('absolute Feuchtigkeit, kg Wasser/kg tr. Luft')
        plt.ylabel('h$_1+x$')
        plt.plot(self.x, self.h1x - 2500 * self.x, marker = next(marker), label = self.name)
        plt.legend()
        plt.ion()
        plt.show()

    def Connect_States(self, endpoint):
        # to be used when a first point is already drawn in a Mollier Diagram
        def drawArrow(A, B):
            plt.arrow(A[0], A[1], B[0] - A[0], B[1] - A[1],
                      width=0.001, length_includes_head=True)
        
        marker = itertools.cycle(('o', 's', 'v', 'p', '*'))
        plt.plot(endpoint.x, endpoint.h1x - 2500 * endpoint.x, marker = next(marker), label = endpoint.name)
        plt.legend()
        
        pt1 = np.array([self.x, self.h1x - 2500 * self.x])
        pt2 = np.array([endpoint.x, endpoint.h1x - 2500 * endpoint.x])
        #print(pt1, pt2)
        drawArrow(pt1, pt2)
        plt.ion()
        plt.show()

class Fog_air:
    '''object for moist air state'''
    def __init__(self, theta, x, p, name): # t in °C, abs. humidity, p in bar
        Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")/10**5
        x_sat = 0.622*Psat/(p-Psat)
        self.x = x
        if self.x < x_sat:
            print('Increase x to at least', x_sat, 'State is not in fog-area')
        self.h1x = 1.004*theta + x_sat*(2500+1.86*theta)+ (self.x - x_sat)*4.19*theta
        self.p = p
        self.name = name

    def Mollier_Diagramm(self):
        p = self.p
        def isenthalp(theta,x):
            # Erzeugung von Isenthalp-Kurven
            return 1.004*theta - 2500*x

        def isotherm(theta,x):
            # Erzeugung von Isotherm-Kurven
            return 1.004*theta + x*(2500+1.86*theta)-2500*x

        def isotherm_sat(theta,x,p):
            #Erzeugung von Isotherm-Kurven-Sat
            Psat = Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")/10**5
            x_sat = 0.622*Psat/(p-Psat)
            return 1.004*theta + x_sat*(2500+1.86*theta)+ (x - x_sat)*4.19*theta-2500*x

        def isofeucht(phi,x,p):
            # Erzeugung von Isotherm-Kurven
            Psat = x*p/(phi*(0.622 +x))
            theta = CP.PropsSI("T","P",Psat*10**5,"Q",1,"water")-273.15 #in °C
            return 1.004*theta + x*(2500+1.86*theta)-2500*x
        x = np.linspace(0.001,0.1,200)
        phi = np.linspace(0.02,1,12)
        theta_isenthalp = np.linspace(0.,350.,36)
        theta = np.linspace(0.,200.,41)
        # Nehmen den minimale Wert für die Isotherme als fcn von x
        isotherm_plot = np.minimum([isotherm(theta,x) for x in x],[isotherm_sat(theta,x,p) for x in x])
        #fig = plt.figure()
        plt.plot(x,[isenthalp(theta_isenthalp,x) for x in x],'g--')
        plt.plot(x,[isofeucht(phi,x,p) for x in x],'b--')
        plt.plot(x, isotherm_plot,'r')
        plt.ylim(0,100)
        plt.xlim(0,0.1)
        plt.xlabel('absolute Feuchtigkeit, kg Wasser/kg tr. Luft')
        plt.ylabel('h$_1+x$')
        plt.plot(self.x, self.h1x - 2500 * self.x, 'ro', label = self.name)
        plt.legend()
        plt.ion()
        plt.show()

    def Connect_States(self, endpoint):
        def drawArrow(A, B):
            plt.arrow(A[0], A[1], B[0] - A[0], B[1] - A[1],
                      width=0.001, length_includes_head=True)
        
        plt.plot(endpoint.x, endpoint.h1x - 2500 * endpoint.x, 'ro', label = endpoint.name)
        plt.legend()
        
        pt1 = np.array([self.x, self.h1x - 2500 * self.x])
        pt2 = np.array([endpoint.x, endpoint.h1x - 2500 * endpoint.x])
        #print(pt1, pt2)
        drawArrow(pt1, pt2)
        plt.ion()
        plt.show()

class Water:
    '''object for pure water states, liquid and vapor'''
    def __init__(self, theta, p): # t in °C, p in bar
        self.t_sat = CP.PropsSI("T","P",p * 1e5,"Q",1,"water") - 273.15 # Tsat in °C
        if theta > self.t_sat:
            self.h_w = hfg0 + CpW * theta
        if theta < self.t_sat:
            self.h_w = CW * theta        
        
def x(theta,p,phi):
    Psat = Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")/10**5
    return 0.622 * phi*Psat/(p - phi*Psat)# kg Wasser/kg tr. Luft

def CP_x(theta,p,phi):
    return CP.HAPropsSI('HumRat','Tdb',theta + 273.15,'RelHum',phi,'P',p*1e5)

def h1x_x(theta,p,x):
    Psat = Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")/10**5
    x_sat = 0.622*Psat/(p-Psat)
    if x < x_sat:
        return CpL*theta + x*(hfg0+CpW*theta) # kJ/kg tr. Luft
    else:
        return CpL*theta + x_sat*(hfg0+CpW*theta)+ (x - x_sat)*4.19*theta

def CP_h1x_x(theta,p,x):
    return CP.HAPropsSI('H','Tdb',theta + 273.15,'HumRat',x,'P',p*1e5)/1000 # kJ/kg tr.
    
def v1x(theta,p,phi):
    Psat = Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")/10**5
    T = theta + 273.15
    return 287 * T/(p*1e5) * (1 + (28.97/18.02) * x(theta,p,phi))

def CP_v1x(theta,p,phi):
    return CP.HAPropsSI('Vda','Tdb',theta + 273.15,'RelHum',phi,'P',p*1e5)

def h1x_phi(theta,p,phi):
    Psat = Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")/10**5
    x_sat = 0.622*Psat/(p-Psat)
    x = 0.622 * phi*Psat/(p - phi*Psat)
    return CpL*theta + x*(hfg0+CpW*theta) # kJ/kg tr. Luft

def CP_h1x_phi(theta,p,phi):
    return CP.HAPropsSI('H','Tdb',theta + 273.15,'RelHum',phi,'P',p*1e5)/1000 # kJ/kg tr.

def psychro_x(theta_dry,theta_wet,p):
    h_w = CW * theta_wet
    Psat = CP.PropsSI("P","T",theta_wet+273.15,"Q",1,"water")/10**5
    x_sat = 0.622*Psat/(p-Psat)
    return (CpL*(theta_wet-theta_dry)+x_sat*((hfg0+1.86*theta_wet)-h_w))/(hfg0+CpW*theta_dry-h_w)

def CP_psychro_x(theta_dry,theta_wet,p):
    return CP.HAPropsSI('HumRat','Tdb',theta_dry + 273.15,'Twb',theta_wet + 273.15,'P',p*1e5) # kg/kg tr.

def psychro_phi(theta_dry,theta_wet,p):
    h_w = CW * theta_wet
    Psat = CP.PropsSI("P","T",theta_wet+273.15,"Q",1,"water")/10**5
    x_sat = 0.622*Psat/(p-Psat)
    x = (CpL*(theta_wet-theta_dry)+x_sat*((hfg0+1.86*theta_wet)-h_w))/(hfg0+CpW*theta_dry-h_w)
    Psat = CP.PropsSI("P","T",theta_dry + 273.15,"Q",1,"water")/10**5
    return x*p/((0.622 + x)*Psat)

def CP_psychro_phi(theta_dry,theta_wet,p):
    return CP.HAPropsSI('RelHum','Tdb',theta_dry + 273.15,'Twb',theta_wet + 273.15,'P',p*1e5)

def psychro_h1x(theta_dry,theta_wet,p):
    h_w = CW * theta_wet
    Psat = CP.PropsSI("P","T",theta_wet+273.15,"Q",1,"water")/10**5
    x_sat = 0.622*Psat/(p-Psat)
    x = (CpL*(theta_wet-theta_dry)+x_sat*((hfg0+CpW*theta_wet)-h_w))/(hfg0+CpW*theta_dry-h_w)
    return CpL*theta_dry + x*(hfg0+CpW*theta_dry) # kJ/kg tr. Luft

def CP_psychro_h1x(theta_dry,theta_wet,p):
    return CP.HAPropsSI('H','Tdb',theta_dry + 273.15,'Twb',theta_wet+273.15,'P',p*1e5)/1000 # kJ/kg tr. Luft

def phi_xT(theta,x,p):
    Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")/10**5
    return x*p/((0.622 + x)*Psat)

def CP_phi_xT(theta,x,p):
    return CP.HAPropsSI('RelHum','Tdb',theta + 273.15,'HumRat',x,'P',p*1e5)

def dew_T(theta, phi):
    Psat = CP.PropsSI("P","T",theta+273.15,"Q",1,"water")
    Pw = phi * Psat
    if Pw < 611.655:
        return CP.HAProps("Tdp","T", theta+273.15, "P", 1e5, "R", phi) - 273.15 # in °C
    return CP.PropsSI("T","P",Pw,"Q",1,"water") - 273.15

def CP_dew_T(theta, phi):
    return CP.HAPropsSI('Tdp','Tdb',theta + 273.15,'RelHum',phi,'P',1e5) - 273.15 # in °C

def wet_bulb_T(theta_dry, phi, p):
    Psat = CP.PropsSI("P","T",theta_dry+273.15,"Q",1,"water")/10**5
    x = 0.622 * (phi * Psat)/(p - Psat)
    if phi * Psat < 611.655/1e5:
        return #CP.HAProps("Twb","T", theta_dry+273.15, "P", p * 1e5, "R", phi) - 273.15 # in °C
    def Balance(theta_wet):
        return h1x_phi(theta_dry,p,phi) + (x - psychro_x(theta_dry,theta_wet,p))*CW*theta_wet -\
               psychro_h1x(theta_wet, theta_wet,p)
    theta_wet = theta_dry -5    # intial guess for wet bulb temperature
    return optimize.fsolve(Balance, theta_wet)

def CP_wet_bulb_T(theta_dry, phi, p):
    return CP.HAPropsSI('Twb','Tdb',theta_dry + 273.15,'RelHum',phi,'P',p*1e5) - 273.15 # in °C

