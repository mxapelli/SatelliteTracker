
import math
from coord import *
from calc import *

#Function that computes the visibility area of the user
def visibilidadObs(a, latObs, longObs, altObs, name):
    # Visibilidad observador
    pi = math.pi
    REarth = 6378.e3;                    #meters
    h = a-REarth #Orbital Height Satellite

    if ("ISS" in name or "IRIDIUM" in name or "STARLINK" in name) and (abs(longObs)<150) and (abs(latObs)<70):             
        Emin = 10                                            #min. elevation angle (degrees)
        Emin = Emin * pi / 180
        gamma = math.asin(math.cos(Emin) * REarth / (REarth + h))
        beta = pi / 2 - gamma - Emin
        r = math.sin(beta) * (REarth + h) / math.cos(Emin)
        area = []
        areapos=[]
        areaneg=[]
        userECEF = LLA2ECEF(latObs, longObs, altObs)
        for i in range(1,362):
            a = (i - 1) * pi / 180
            NEDvis = ELEV_AZ2NED(r, Emin, a)
            ECEFvis = NED2ECEF(NEDvis, latObs*pi/180, longObs*pi/180)         #rotate from NED to ECEF
            ECEFvis = sum(userECEF,ECEFvis)  #add ECEF coords. of observer
            lla= ECEF2LLA(ECEFvis[0], ECEFvis[1], ECEFvis[2])
            LAT=lla[0]
            LON=lla[1]
            area.append([LAT,LON])
    else:
        Emin = 10                                            #min. elevation angle (degrees)
        Emin = Emin * pi / 180
        gamma = math.asin(math.cos(Emin) * REarth / (REarth + h))
        beta = pi / 2 - gamma - Emin
        r = math.sin(beta) * (REarth + h) / math.cos(Emin)
        area = []
        areapos=[]
        areaneg=[]
        userECEF = LLA2ECEF(latObs, longObs, altObs)
        for i in range(1,362):
            a = (i - 1) * pi / 180
            NEDvis = ELEV_AZ2NED(r, Emin, a)
            ECEFvis = NED2ECEF(NEDvis, latObs*pi/180, longObs*pi/180)         #rotate from NED to ECEF
            ECEFvis = sum(userECEF,ECEFvis)  #add ECEF coords. of observer
            lla= ECEF2LLA(ECEFvis[0], ECEFvis[1], ECEFvis[2])
            LAT=lla[0]
            LON=lla[1]
            if LON>0:
                if len(areapos)>0 and abs(areapos[-1][1]-LON)>100 and latObs>=0:
                    areapos.append([90,areapos[-1][1]])
                    areapos.append([90,0])
                    areapos.append([LAT,0])
                    areapos.append([LAT, LON])
                elif len(areapos)>0 and abs(areapos[-1][1]-LON)>100 and latObs<0:
                    areapos.append([areapos[-1][0],0])
                    areapos.append([-90,0])
                    areapos.append([-90,LON])
                    areapos.append([LAT, LON])
                elif len(areapos)>0 and abs(areapos[-1][0]-LAT)>10 and abs(LON)<4:
                    areapos.append([areapos[-1][0],0])
                    areapos.append([LAT,0])
                    areapos.append([LAT,LON])
                elif len(areapos)>0 and abs(areapos[-1][0]-LAT)>10 and abs(LON)>170:
                    areapos.append([areapos[-1][0],180])
                    areapos.append([LAT,180])
                    areapos.append([LAT,LON])
                else:
                    areapos.append([LAT, LON])
            else:
                if len(areaneg)>0 and abs(areaneg[-1][1]-LON)>100 and latObs>=0:
                    areaneg.append([areaneg[-1][0],0])
                    areaneg.append([90,0])
                    areaneg.append([90,-180])
                    areaneg.append([LAT, -180])
                elif len(areaneg)>0 and abs(areaneg[-1][1]-LON)>100 and latObs<0:
                    areaneg.append([-90,areaneg[-1][1]])
                    areaneg.append([-90,0])
                    areaneg.append([LAT,0])
                    areaneg.append([LAT, LON])
                elif len(areaneg)>0 and abs(areaneg[-1][0]-LAT)>10 and abs(LON)<4:
                    areaneg.append([areaneg[-1][0],0])
                    areaneg.append([LAT,0])
                    areaneg.append([LAT,LON])
                elif len(areaneg)>0 and abs(areaneg[-1][0]-LAT)>10 and abs(LON)>173:
                    areaneg.append([areaneg[-1][0],-180])
                    areaneg.append([LAT,-180])
                    areaneg.append([LAT,LON])
                else:
                    areaneg.append([LAT, LON])
         
        if len(areaneg)>1 and abs(areaneg[-1][1]-areaneg[0][1])>100 and latObs>=0:
            areaneg.append([areaneg[-1][0],0])
            areaneg.append([90,0])
            areaneg.append([90,areaneg[0][1]])
        if len(areaneg)>1 and abs(areaneg[-1][1]-areaneg[0][1])>100 and latObs<0:
            areaneg.append([-90,areaneg[-1][1]])
            areaneg.append([-90,0])
            areaneg.append([areaneg[0][0],0])
        if len(areaneg)>1 and abs(areaneg[-1][0]-areaneg[0][0])>5 and abs(areaneg[-1][1])<5:
            areaneg.append([areaneg[-1][0],0])
            areaneg.append([areaneg[0][0],0])
        if len(areaneg)>1 and abs(areaneg[-1][0]-areaneg[0][0])>5 and abs(areaneg[-1][1])>175:
            areaneg.append([areaneg[-1][0],-180])
            areaneg.append([areaneg[0][0],-180])
        
        

        if len(areapos)>1 and abs(areapos[-1][1]-areapos[0][1])>100 and latObs>=0:
            areapos.append([areapos[-1][0],180])
            areapos.append([90,180])
            areapos.append([90,0])
            areapos.append([areapos[0][0],0])
        if len(areapos)>1 and abs(areapos[-1][1]-areapos[0][1])>100 and latObs<0:
            areapos.append([areapos[-1][0],0])
            areapos.append([-90,0])
            areapos.append([-90,areapos[0][1]])
        if len(areapos)>1 and abs(areapos[-1][0]-areapos[0][0])>5 and abs(areapos[-1][1])<5:
            areapos.append([areapos[-1][0],0])
            areapos.append([areapos[0][0],0])
        if len(areapos)>1 and abs(areapos[-1][0]-areapos[0][0])>5 and abs(areapos[-1][1])>175:
            areapos.append([areapos[-1][0],180])
            areapos.append([areapos[0][0],180])
        

        if len(areapos)>0:
            area.append(areapos)   
        if len(areaneg)>0:
            area.append(areaneg) 
           #visibility zone
    
    return area

