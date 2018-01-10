
import os
import sys
import numpy as np


def _getXmax(primarytype, energy, zen2):
    # type of primary (electron or pion, energy in EeV, zenith (GRAND) in rad
    # factor approximated from https://pos.sissa.it/301/1100/pdf
    if primarytype=='electron': # aprroximated by gamma shower
        a=82.5 # g/cm2
        c=342.5 #g/cm2
    if primarytype=='pion' or primarytype=='proton': # pion, kaon .... aprroximated by proton
        a=62.5 # g/cm2
        c=357.5 #g/cm2
    if primarytype=='iron': # aprroximated by proton
        a=60 # g/cm2 # just approximated 
        c=177.5 #g/cm2
        
    Xmax= a*np.log10(energy*10**6.)+c # E/EeV* 10**6. to be in TeV

    return Xmax#/abs(np.cos(np.pi-zen2)) # TODO: how to correct for slanted shower

def _dist_decay_Xmax(zen2, injh2, Xmax_primary): #zen2: zenith of target shower
    #% Using isothermal Model
    rho_0 = 1.225*0.001#; % kg/m3 to 0.001g/cm3: 1g/cm3=1000kg/m3, since X given in g/cm2
    M = 0.028966#;  %kg/mol - 1000g/mol
    g = 9.81#; %ms-2
    T = 288.#; % K
    R = 8.32#; J/K/mol , J=kg m2/s2

    hD=injh2
    step=10 #m
    if hD>10000:
        step=100 #m
    Xmax_primary= Xmax_primary#* 10. # g/cm2 to kg/m2: 1g/cm2 = 10kg/m2
    gamma=np.pi-zen2 # counterpart of where it goes to
    Re= 6370949 # m, Earth radius
    X=0.
    i=0.
    h=hD
    ai=0
    while X< Xmax_primary:
        i=i+1
        ai=i*step #100. #m
        hi= -Re+np.sqrt(Re**2. + ai**2. + hD**2. + 2.*Re*hD - 2*ai*np.cos(gamma) *(Re+hD))## cos(gamma)= + to - at 90dg
        deltah= abs(h-hi) #(h_i-1 - hi)= delta h
        h=hi # new height
        X=X+ rho_0*np.exp(-g*M*hi/(R*T)) * ai*100. #(deltah*100) *abs(1./np.cos(np.pi-zen2)) # Xmax in g/cm2, slanted = Xmax, vertical/ cos(theta); density in g/cm3, h: m->100cm, np.pi-zen2 since it is defined as where the showers comes from, abs(cosine) so correct for minus values

    return h, ai # Xmax_height in m, Xmax_distance in m