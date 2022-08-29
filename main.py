import os
from time import time
from flask import Flask, render_template, request, session, redirect, url_for
from dotenv import load_dotenv
import json
import os
import pymongo
import matplotlib
import matplotlib.patches as mpatches
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO
import base64
from datetime import datetime,timedelta
from apscheduler.schedulers.background import BackgroundScheduler

import time
import calendar
import aiohttp
import asyncio
import urllib.request
import math
import numpy
import requests


load_dotenv() # use dotenv to hide sensitive credential as environment variables
DATABASE_URL=f'mongodb+srv://MXapelli:{os.environ.get("passwordDB")}'\
	      '@mongo-heroku-cluster.aiiqhhv.mongodb.net/satellites?retryWrites=true&w=majority'# get connection url from environment

client=pymongo.MongoClient(DATABASE_URL) # establish connection with database
mongo_db=client.db

app = Flask(__name__)

app.secret_key = os.environ.get("passwordSess")

@app.route('/',methods = ['POST', 'GET'])
def index():
    if request.method == 'POST':
      coordsData = request.form['nm']
      print("Coordenadas ",coordsData)
      if isinstance(coordsData, str):
        coords=coordsData.split("/")
        urlAlt="https://api.open-elevation.com/api/v1/lookup?locations="+coords[0]+","+coords[1]
        result = requests.get(urlAlt)
        alt=result.json()
        altM=alt["results"][0]['elevation']
        session['latUser'] = float(coords[0])
        session['longUser'] = float(coords[1])
        session['altUser'] = float(altM)
        return redirect(url_for('index'))

    jsonSat=[]
    for sat in mongo_db.satellites.find():
        jsonSat.append({"OBJECT_NAME": str(sat["name"]), "NORAD_CAT_ID": str(sat["noradID"])})

    return render_template('index.html',datos=jsonSat)

@app.route('/<int:catnr>')
def satellite(catnr):
    xmap=[] 
    ymap=[]
    zmap=[]
    pi=math.pi
    incTime=30
    Ge = 6.67384*10**(-11)					#Gravitational constant
    Me = 5.972*10**24
    if 'latUser' in session:	
        latUser=session['latUser']
        longUser=session['longUser']
        altUser=session['altUser']
    else:
        latUser = 41.2757116961354
        longUser= 1.9872286269285848
        altUser= 4

    sat=mongo_db.satellites.find_one({"noradID": catnr})
    if sat is None:
        print("Error satellite does not exist")
        return render_template('error.html')

    name=sat['name']
    print("The user has selected satellite: ",name)
    epoch=sat['epoch']
    epochT=str(epoch).split("T")
    epochT[1]=epochT[1].split(".")
    hour=epochT[1][0]
    incl=sat['incl']
    omega=sat['omega']
    ecc=sat['ecc']
    w=sat['w']
    M=sat['M']
    n=sat['n']

    #TIME TO TOA
    date=epochT[0]
    x=time.strptime(date+" "+ hour, "%Y-%m-%d %H:%M:%S") 
    s=x.tm_yday
    s=s+x.tm_hour/24+x.tm_min/(24*60)+x.tm_sec/(24*60*60)
    ToA=s
    result=esecGAST(ToA)
    esec=result[0]
    GASTdegToA=result[1]
    
    #Process satellite kep. elements
    T = 86400 / n                    #orbit period (secs.)
    n = 2 * pi / T                          #mean movement [rad/s]
    a = (Ge * Me / n**2)**(1/3)                   #semi-major axis [m]
    io =  incl * pi / 180                   #orbit inclination [rad]
    ecc =  ecc                     #orbit eccentricity
    #GAST at a given epoch is the RA of the Greenwich meridian at that epoch
    Oo = (omega - GASTdegToA) * pi / 180      #Longitude of AN at ToA [rad] Aqui hay que poner GASTdegToA
    dO = 0                                  #Rate of Right ascension
    w =  w * pi / 180                   #argument of perigee [rad]
    Mo =  M* pi / 180                   #Mean anomaly at ToA [rad]

    #Visibility Area of User
    visCoord=visibilidadObs(a,latUser,longUser,altUser);
    Xvis=visCoord[0]
    Yvis=visCoord[1]

    #Groundtrack for NT periods
    coordinates=computeCoordinates(esec,T,n,a,io,ecc,Oo,dO,w,Mo,incTime)
    xmap=coordinates[0]
    ymap=coordinates[1]
    zmap=coordinates[2]
    latmap=coordinates[3]
    longmap=coordinates[4]
    altmap=coordinates[5]

    #Check Visibility of satellite
    elev=[]
    for i in range(len(xmap)):
        satECEF=[xmap[i], ymap[i], zmap[i]]
        userECEF=LLA2ECEF(latUser,longUser,altUser)
        pointV=sub(satECEF,userECEF)
        NED=ECEF2NED(pointV,latUser*pi/180,longUser*pi/180)
        d=math.sqrt(NED[0]**2+NED[1]**2+NED[2]**2)
        alpha= math.asin((-NED[2])/d)*180/pi
        if(alpha>=10):
            elev.append(alpha)
    if len(elev)>0:
        vis=1
    else:
        vis=0
    
    atime=time.localtime()
    st=time.strftime("%a, %d %b %Y %H:%M:%S ",atime)
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    session['time']=dt_string
    print ("session info",dt_string)
    text=[("You have selected the "+name+" satellite with Catalog Number: "+str(catnr)),("The local time is: " +st)]
    return render_template('satellite.html',entries=text,longs=longmap,lats=latmap,worldLa=worldLat,worldLo=worldLong,number=str(catnr),userlat=latUser,userlong=longUser,xvis=Xvis,yvis=Yvis,vis=vis)

@app.route('/<int:catnr>/doppler')
def doppler(catnr):
    pi=math.pi
    Ge = 6.67384*10**(-11)					#Gravitational constant
    Me = 5.972*10**24
    if 'latUser' in session:	
        latUser=session['latUser']
        longUser=session['longUser']
        altUser=session['altUser']
    else:
        latUser = 41.2757116961354
        longUser= 1.9872286269285848
        altUser= 4

    sat=mongo_db.satellites.find_one({"noradID": catnr})
    if sat is None:
        return render_template('error.html')
    name=sat["name"]
    freq=sat["freq"]
    if freq==1575.42*10**6:
        incTime=30
    else:
        incTime=10

    datosObtenidos=requests.get('https://celestrak.com/NORAD/elements/gp.php?CATNR='+str(catnr)+'&FORMAT=JSON')
    check=jsonCheck(datosObtenidos)
    if check=="error":
        print("Error satellite does not exist")
        return render_template('error.html')
    datosSat= datosObtenidos.json()
    print("The user has selected satellite: ",datosSat[0]['OBJECT_NAME'])
    epoch=datosSat[0]['EPOCH']
    epochT=str(epoch).split("T")
    epochT[1]=epochT[1].split(".")
    hour=epochT[1][0]
    incl=datosSat[0]['INCLINATION']
    omega=datosSat[0]['RA_OF_ASC_NODE']
    ecc=datosSat[0]['ECCENTRICITY']
    w=datosSat[0]['ARG_OF_PERICENTER']
    M=datosSat[0]['MEAN_ANOMALY']
    n=datosSat[0]['MEAN_MOTION']

    if sat is not None:
        name=sat["name"]
        if sat["epoch"]!=epoch:
            id = { "noradID": catnr }
            newvalues = { "$set": { "epoch": epoch,"incl": incl,"ecc": ecc,"omega": omega,"w": w,"M":M,"n":n} }
            mongo_db.satellites.update_one(id,newvalues)
            print(name," has been updated")


    #TIME TO TOA
    date=epochT[0]
    x=time.strptime(date+" "+ hour, "%Y-%m-%d %H:%M:%S") 
    s=x.tm_yday
    s=s+x.tm_hour/24+x.tm_min/(24*60)+x.tm_sec/(24*60*60)
    ToA=s
    result=esecGAST(ToA)
    esec=result[0]
    GASTdegToA=result[1]
    
    #Process satellite kep. elements
    T = 86400 / n                    #orbit period (secs.)
    n = 2 * pi / T                          #mean movement [rad/s]
    a = (Ge * Me / n**2)**(1/3)                   #semi-major axis [m]
    io =  incl * pi / 180                   #orbit inclination [rad]
    ecc =  ecc                     #orbit eccentricity
    #GAST at a given epoch is the RA of the Greenwich meridian at that epoch
    Oo = (omega - GASTdegToA) * pi / 180      #Longitude of AN at ToA [rad] Aqui hay que poner GASTdegToA
    dO = 0                                  #Rate of Right ascension
    w =  w * pi / 180                   #argument of perigee [rad]
    Mo =  M* pi / 180                   #Mean anomaly at ToA [rad]

    #Groundtrack for NT periods
    coordinates=computeCoordinates(esec,T,n,a,io,ecc,Oo,dO,w,Mo,incTime)
    xmap=coordinates[0]
    ymap=coordinates[1]
    zmap=coordinates[2]


    print("Frequency of satellite ",freq)
    r=(1-ecc)*a
    v=math.sqrt(Ge*Me*((2/r)-(1/a)))
    Lm=r*v
    GAST1=GAST(0)
    t4=round(len(xmap)/4)
    GAST2=GAST(incTime*t4)
    ECEF1=ECEF2ECI(xmap[0],ymap[0],zmap[0],GAST1*math.pi/180)
    ECEF2=ECEF2ECI(xmap[t4],ymap[t4],zmap[t4],GAST2*math.pi/180)

    # Orbital Plane Vector
    C= cross(ECEF1,ECEF2)
    V=norm(C)
    LmV=prod(V,Lm)

    #Doppler
    vSatV=[]
    ECI=[]
    GastV=[]
    for n in range(len(xmap)):
        t=n*incTime
        GAST1=GAST(t)
        GastV.append(GAST1)
        rV=ECEF2ECI(xmap[n],ymap[n],zmap[n],GAST1*math.pi/180)
        ECI.append(rV)
        cros = cross(LmV,rV)
        absrV=math.sqrt(rV[0]**2+rV[1]**2+rV[2]**2)**2
        vSat=div(cros,absrV)
        vSatV.append(vSat)
        vSatAbs=math.sqrt(vSat[0]**2+vSat[1]**2+vSat[2]**2)


    #Vector Velocidad Observador
    Rearth=6378*1000;
    Tearth=(23*60+56)*60+4;#segons
    vearth=2*math.pi*(Rearth/Tearth)
    latObs=latUser
    longObs=longUser
    altObs=altUser

    dSatObs=[]
    dObs=[]
    vinst=[]
    vinstReal=[]
    fD=[]
    fDReal=[]
    for n in range(len(ECI)):

        # Velocidad Observador en Funcion de la Latitud
        Gast=GastV[n]
        GaLo=(Gast+longObs)*math.pi/180
        vObs=[vearth*-math.sin(GaLo),vearth*math.cos(GaLo),0]
        vObsLat=prod(vObs,math.cos(latObs*math.pi/180))

        # Velocidad Satelite respecto Observador
        vResta=sub(vSatV[n],vObsLat)

        # Posicion Satelite respecto Observador
        userECEF=LLA2ECEF(latObs,longObs,altObs)
        #Checkpoint Todo Correcto print(userECEF)
        userECI=ECEF2ECI(userECEF[0],userECEF[1],userECEF[2],Gast*math.pi/180)
        satV=sub(userECI,ECI[n])
        #print(satV)
        dSatObs.append(modulo(satV))
        dObs.append(modulo(userECEF))
        satVnorm=norm(satV)
        vinst.append(dot(vSatV[n],satVnorm))
        vinstReal.append(dot(vResta,satVnorm))
        fD.append(freq*vinst[n]/(3*10**8))
        fDReal.append((freq*vinstReal[n])/(3*10**8))

    print ("Max Freq. Doppler: ", max(fDReal))
    print ("Max Speed ", max(vinstReal))

    id = { "noradID": catnr }
    newvalues = { "$set": { "vDoppler": fDReal} }
    mongo_db.satellites.update_one(id,newvalues)


    #Visibility of satellite
    elev=[]
    az=[]
    ind=[]
    for i in range(len(xmap)):
        satECEF=[xmap[i], ymap[i], zmap[i]]
        userECEF=LLA2ECEF(latUser,longUser,altUser)
        pointV=sub(satECEF,userECEF)
        NED=ECEF2NED(pointV,latUser*pi/180,longUser*pi/180)
        d=math.sqrt(NED[0]**2+NED[1]**2+NED[2]**2)
        alpha= math.asin((-NED[2])/d)*180/pi
        beta= math.atan2(NED[1],NED[0])*180/pi;
        if(alpha>=10):
            elev.append(alpha)
            ind.append(i)
            if beta<0:
                az.append((beta+360)*pi/(180))
            else:
                az.append(beta*pi/180)

    dopplerVis=[]
    timeP=[]
    cat=[]
    c=0
    catN=0
    change=[]
    pointChange=0
    for i in range(len(ind)):
        dopplerVis.append(fDReal[ind[i]])
        timeP.append(ind[i]*incTime)
        if (c!=ind[i]-1 and c!=0):
            catN=1
            #To detect visibility changes
            change.append(ind[i-1]*incTime)
            change.append(ind[i]*incTime)
            pointChange=i
        cat.append(catN)
        c=ind[i]

    now = datetime.now()
    actual_time = now.strftime("%d/%m/%Y %H:%M:%S")
    if len(timeP)<1:
        return  render_template('errorDoppler.html')
    init=timeP[0]
    fin = timeP[-1]
    duration=(fin-init)/7

    #It convert the strings to datetime obj
    actual_time_obj = datetime.strptime(actual_time, '%d/%m/%Y %H:%M:%S')
    times=[]
    days=[]
    for i in range(8):
        t=init+duration*(i)
        time_obj=actual_time_obj + timedelta(seconds=t)
        timeText=str(time_obj).split()
        timev=timeText[1].split('.')
        times.append(timev[0])
        if (i==0):
            days.append(time_obj.strftime("%d %b"))
    
    
    visStep=[]
    visStep.append(times[0])
    if len(change)>0:
        fin1 = actual_time_obj + timedelta(seconds=change[0])
        days.append(fin1.strftime("%d %b"))
        in2=actual_time_obj+timedelta(seconds=change[1])
        days.append(in2.strftime("%d %b"))
        timeText=str(fin1).split()
        timev=timeText[1].split('.')
        visStep.append(timev[0])
        timeText=str(in2).split()
        timev=timeText[1].split('.')
        visStep.append(timev[0])
    visStep.append(times[-1])
    fin2 = actual_time_obj + timedelta(seconds=fin)
    days.append(fin2.strftime("%d %b"))
    visText="The satellite will be visible from "+visStep[0]+" ("+days[0]+") "+" to "+ visStep[1]+" ("+days[1]+") "
    if len(visStep)>2:
        visText=visText+ " and from "+visStep[2]+" ("+days[2]+") "+"to "+ visStep[3]+" ("+days[3]+")"+"."
    


    #New Polar Chart Visibility
    img = BytesIO()
    colormap = numpy.array(['yellow','red'])

    fig = plt.figure(facecolor='#383A3F')
    ax = fig.add_subplot(projection='polar')
    ax.grid(c='white')
    d = ax.scatter(az, elev, c=colormap[cat], s=30, alpha=1)
    if len(change)>0:
        azfirst=[]
        elevfirst=[]
        azfirst.append(az[0])
        azfirst.append(az[pointChange])
        elevfirst.append(elev[0])
        elevfirst.append(elev[pointChange])
        e = ax.scatter(azfirst, elevfirst, c='black', s=40, alpha=1)
    else:
        e = ax.scatter(az[0], elev[0], c='black', s=40, alpha=1)

    ax.set_title("Visibility of Satellite "+name+" with ID "+str(catnr), va='bottom',c='white')
    ax.set_theta_zero_location('N')
    ax.set_rlabel_position(-90)
    ax.set_theta_direction(-1)  # theta increasing clockwise
    ax.set_xticklabels(['N', '45','E', '135', 'S', '225', 'W', '315'],c='white',weight='bold')
    ax.set_ylim(90, 0)
    ax.set_yticklabels(['0','20','40','60','80'],color='white',weight='bold')
    ax.set_yticks([0,20,40,60,80])
    ax.set_facecolor((0, 0, 1))
    ax.set_axisbelow(True)
    pop_a = mpatches.Patch(color='yellow', label='First Stretch')
    pop_b = mpatches.Patch(color='red', label='Second Stretch')
    pop_c = mpatches.Patch(color='black', label='First point of Visibility')
    if len(change)>0:
        fig.legend(handles=[pop_a,pop_b,pop_c],loc="lower center",ncol=3)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.12)
    else:
        fig.legend(handles=[pop_a,pop_c],loc="lower center",ncol=2)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.12)
    
    plt.savefig(img, format='png',bbox_inches = "tight")
    plt.close()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf8')

    #Doppler new graph
    img2 = BytesIO()
        
    fig,ax = plt.subplots(facecolor='#383A3F')
    ax.grid(c='white')
    k=ax.scatter(timeP, dopplerVis, c=colormap[cat], s=40, alpha=1)
    ax.set_title("Doppler frequency of Satellite "+name+" with ID "+str(catnr), va='bottom',c='white')
    ax.set_xlabel('Time (UTC)',c='white')  # Add an x-label to the axes.
    xticks=numpy.arange(init, fin+1, step=duration)
    ax.set_xticks(xticks,c='white')
    ax.set_xticklabels(times,c='white',rotation = 90)

    # setting ticks for y-axis
    yticks=ax.get_yticks()
    yticklabel=[]
    for i in range(len(yticks)):
        yticklabel.append(yticks[i]/1000)

    ax.set_ylabel('Doppler Frequency (kHz)',c='white')  # Add a y-label to the axes
    ax.set_yticklabels(yticklabel,color='white')
    ax.set_facecolor((0, 0, 1))
    ax.set_axisbelow(True)

    plt.savefig(img2, format='png',bbox_inches = "tight")
    plt.close()
    img2.seek(0)
    plotDoppler_url = base64.b64encode(img2.getvalue()).decode('utf8')

    text=[("The frequency for this satellite is: "+str(freq/10**6)+" MHz"),("Max Doppler frequency: "+str(round(max(dopplerVis)/1000,2))+" kHz"),("Min Doppler frequency: "+str(round(min(dopplerVis)/1000,2))+" kHz"),("Max radial speed: " +str(round(max(vinstReal),2))+ " m/s"),(visText)]
    return render_template('doppler.html',entries=text,doppler=fDReal, plot_url=plot_url,plotDoppler_url=plotDoppler_url)

def test_job():
    satsDB=[]
    satsCelestrak=[]
    start_time=time.time()
    sats=mongo_db.satellites.find()
    for i in sats:
        satsDB.append(i)
    
    constellations = ['GROUP=starlink','GROUP=iridium-next', 'GROUP=gps-ops', 'CATNR=25544','GROUP=iridium']
    for i in range(len(constellations)):
        result_data=requests.get('https://celestrak.com/NORAD/elements/gp.php?'+constellations[i]+'&FORMAT=JSON')
        check=jsonCheck(result_data)
        if check!="error":
            satsCelestrak.append(result_data.json())

    #Updating existing satellites
    for sat in satsDB:
        catnr=sat['noradID']
        name=sat['name']
        found=0
        while found==0:
            for i in range(len(satsCelestrak)):
                for j in range(len(satsCelestrak[i])):
                    noradID=int(satsCelestrak[i][j]['NORAD_CAT_ID'])
                    epoch=satsCelestrak[i][j]['EPOCH']
                    if (catnr==noradID):
                        if "IRIDIUM" in name:
                            name=name.strip("[-]")
                        if sat['epoch']!=epoch:
                            incl=satsCelestrak[i][j]['INCLINATION']
                            omega=satsCelestrak[i][j]['RA_OF_ASC_NODE']
                            ecc=satsCelestrak[i][j]['ECCENTRICITY']
                            w=satsCelestrak[i][j]['ARG_OF_PERICENTER']
                            M=satsCelestrak[i][j]['MEAN_ANOMALY']
                            n=satsCelestrak[i][j]['MEAN_MOTION']
                            id = { "noradID": catnr }
                            newvalues = { "$set": { "name": name,"epoch": epoch,"incl": incl,"ecc": ecc,"omega": omega,"w": w,"M":M,"n":n} }
                            mongo_db.satellites.update_one(id,newvalues)
                            print(name," has been updated")
                        del satsCelestrak[i][j]
                        break
            found=1
    print("Updated Completed")

    #Cleaning list of sats
    i=0
    while i<(len(satsCelestrak)-1):
        if len(satsCelestrak[i])==0:
            del satsCelestrak[i]
            i=i-1
        i=i+1

    #Adding new satellites
    for i in range(len(satsCelestrak)):
        for j in range(len(satsCelestrak[i])):
            name=satsCelestrak[i][j]['OBJECT_NAME']
            catnr=int(satsCelestrak[i][j]['NORAD_CAT_ID'])
            epoch=satsCelestrak[i][j]['EPOCH']
            incl=satsCelestrak[i][j]['INCLINATION']
            omega=satsCelestrak[i][j]['RA_OF_ASC_NODE']
            ecc=satsCelestrak[i][j]['ECCENTRICITY']
            w=satsCelestrak[i][j]['ARG_OF_PERICENTER']
            M=satsCelestrak[i][j]['MEAN_ANOMALY']
            n=satsCelestrak[i][j]['MEAN_MOTION']

            fre=0
            if "ISS" in name:
                fre=145.8*10**6
            elif "STARLINK" in name:
                fre=10950*10**6
            elif "NAVSTAR" in name:
                fre=1575.42*10**6
            elif "IRIDIUM" in name:
                fre=1626.1042*10**6
            elif "GPS" in name:
                fre=1575.42*10**6

            if fre!=0:
                sat={"noradID": catnr, "name": name,"epoch": epoch,"ecc": ecc,"incl": incl,"omega": omega,"w": w,"M": M,
                "n": n,"freq": fre,"xECEF": [],"yECEF": [],"zECEF": [],"vDoppler": []}
                mongo_db.satellites.insert_one(sat)
                print(name," has been inserted in DB")
    print("Added Completed")
    now=datetime.now()
    t=time.time()-start_time
    print("Database updated at",now,"It took",t,"seconds to complete the task.")

scheduler = BackgroundScheduler()
job = scheduler.add_job(test_job, 'interval', minutes=10)
scheduler.start()

def GAST(esec):
    #Computes GAST [deg] at current time + esec [sec]
    #esec = ToA time in secs. from/to now (esec < 0 means ToA is in the past)
    #Find Julian Date now
    unixepoch = 2440587.5                #JD at unix epoch: 0h (UTC) 1/1/1970
    unixseconds=time.time()        
    JD = unixepoch + unixseconds / (86400); 
    JToA = JD + esec / 86400;                             #Julian Date at ToA

    #Greenwich Mean Sidereal Time (GMST) is the hour angle of the average position of the vernal equinox,
    #neglecting short term motions of the equinox due to nutation. GAST is GMST corrected for
    #the shift in the position of the vernal equinox due to nutation.
    #GAST at a given epoch is the RA of the Greenwich meridian at that epoch (usually in time units).

    #Find GAST in degrees at ToA
    J2000 = 2451545.0;                                   # epoch is 1/1/2000 at 12:00 UTC
    midnight = round(JToA) - 0.5;                   # midnight of JToA
    days_since_midnight = JToA - midnight;
    hours_since_midnight = days_since_midnight * 24;
    days_since_epoch = JToA - J2000;
    centuries_since_epoch = days_since_epoch / 36525;
    whole_days_since_epoch = midnight - J2000;
    GAST = 6.697374558 + 0.06570982441908 * whole_days_since_epoch + 1.00273790935 * hours_since_midnight + 0.000026 * centuries_since_epoch**2; #GAST in hours from ?
    GASTh = GAST % 24;                #GAST in hours at ToA
    GASTdeg = 15 * 1.0027855 * GASTh;      #GAST in degrees at ToA (approx. 361º/24h)
    return GASTdeg

def esecGAST(ToA):
    unixepoch = 2440587.5                #JD at unix epoch: 0h (UTC) 1/1/1970
    unixseconds=time.time()        
    JD = unixepoch + unixseconds / (86400);         #Julian Date now in UTC

    #Compute current time in secs. from ToA
    MJ2022 = 59580#MJD on 1/1/2022 at 00:00 UTC (see http://leapsecond.com/java/cal.htm)
    J2022 = MJ2022 + 2400000.5 #JD on 1/1/2022 at 00:00 UTC
    JToA = J2022 + ToA - 1#ToA in JD 
    esec = 86400 * (JD - JToA) #temps transcorregut des de'l ToA en segons

    #Find GAST in degrees at ToA
    J2000 = 2451545.0                                #epoch is 1/1/2000 at 12:00 UTC
    midnight = round(JToA) - 0.5                       #midnight of JToA
    days_since_midnight = JToA - midnight
    hours_since_midnight = days_since_midnight * 24
    days_since_epoch = JToA - J2000
    centuries_since_epoch = days_since_epoch / 36525
    whole_days_since_epoch = midnight - J2000
    GAST = 6.697374558 + 0.06570982441908 * whole_days_since_epoch + 1.00273790935 * hours_since_midnight + 0.000026 * centuries_since_epoch**2 #GAST in hours from ?
    GASTh = GAST%24       #GAST in hours at ToA
    GASTdegToA = 15 * 1.0027855 * GASTh     #GAST in degrees at ToA (approx. 361º/24h)
    esecGAST=[]
    esecGAST.append(esec)
    esecGAST.append(GASTdegToA)
    return esecGAST

def ECEF2ECI(X, Y, Z, B):
    ECI=[]
    Xe=math.cos(B)*X-math.sin(B)*Y
    Ye=math.sin(B)*X+math.cos(B)*Y
    Ze=Z
    ECI.append(Xe)
    ECI.append(Ye)
    ECI.append(Ze)
    return ECI

def LLA2ECEF(lat,lon,alt):
    a=6378137
    e2=0.00669437999014
    e=math.sqrt(e2)
    lat=lat*math.pi/180
    lon=lon*math.pi/180
    
    coco=math.sqrt(1-e2*math.sin(lat)**2)
    clat=math.cos(lat)
    clong=math.cos(lon)
    slat=math.sin(lat)
    slong=math.sin(lon)
    x=(a/coco+alt)*clat*clong;
    y=(a/coco+alt)*clat*slong;
    z=(a*(1-e2)/(coco)+alt)*slat;
    LLA=[x,y,z]
    return LLA

def ECEF2NED(ECEF,phi,lamda):

    #1) Rotation Matrix from NED to ECEF: ECEF = M*NED
    r=-math.sin(phi)*math.cos(lamda)
    M = numpy.array([[(-math.sin(phi)*math.cos(lamda)),(-math.sin(lamda)),(-math.cos(phi)*math.cos(lamda))],[(-math.sin(phi)*math.sin(lamda)),(math.cos(lamda)),(-math.cos(phi)*math.sin(lamda))],[(math.cos(phi)),0,(-math.sin(phi))]],dtype='f')
    #2) Compute NED coordinates
    #print(M[0][0])
    #print(r)
    ECEFt=numpy.array([ECEF]).T
    invM=numpy.linalg.inv(M)

    NED=numpy.matmul(invM,ECEFt)
    return NED

def computeCoordinates(esec,T,n,a,io,ecc,Oo,dO,w,Mo,incTime):
    t = esec
    NT = 1
    j=0
    xmap=[]
    ymap=[]
    zmap=[]
    latmap=[]
    longmap=[]
    altmap=[]
    R=[]
    while t < esec + T * NT:

        #ECEF coordinates
        #Kepler2ECEF ##[x, y, z] = Kepler2ECEF(a, io, e, Oo, dO, w, Mo, n, t)
        M = Mo + n * t #Mean anomaly now [rad]
        #Solve Keppler equation for eccentric anomally iteratively
        # Loop is exited when we reach the desired tolerance
        tol = 10**(-8) #Converging tolerance
        E = M
        while True:
            Eo = E
            E = M + ecc * math.sin(Eo)
            if abs(Eo - E) < tol: 
                break

        #True anomaly
        cosv = (math.cos(E) - ecc) / (1 - ecc * math.cos(E))
        sinv = math.sqrt(1 - ecc**2) * math.sin(E) / (1 - ecc * math.cos(E))
        v = math.atan2(sinv, cosv)

        #Cartesian coordinates within the orbital plane
        u = v + w
        r = a * (1 - ecc * math.cos(E))
        O = Oo - 7.2921151467e-5 * t + dO * t #corrected right ascension (7.2921151467e-5 = rate of Earth rotation [rad/s])
        xp = r * math.cos(u)
        yp = r * math.sin(u)
        
        #ECEF coordinates [m]
        x = xp * math.cos(O) - yp * math.cos(io) * math.sin(O)
        y = xp * math.sin(O) + yp * math.cos(io) * math.cos(O)
        z = yp * math.sin(io)

        xmap.append(x)
        ymap.append(y)
        zmap.append(z)
        R.append(math.sqrt(x**2+y**2+z**2))

        #LLA coordinates [lat, lon, alt] = ECEF2LLA(x, y, z)
        #Convert from cartesian (x,y,z) to LLA coordinates in datum WGS84
        A=6378137 
        E2=0.00669437999014  #datos elipsoide WGS84
        B=A*math.sqrt(1-E2) 
        EP=A*math.sqrt(E2)/B
        
        r=math.sqrt(x**2+y**2)
        EE2=A**2-B**2
        F=54*B**2*z**2
        G=r**2+(1-E2)*z**2-E2*EE2
        C=E2**2*(F*r**2)/G**3
        S=(1+C+math.sqrt(C**2+2*C))**(1/3)
        P=F/(3*G**2*(S+1/S+1)**2)
        Q=math.sqrt(1+2*E2**2*P)
        ro=-E2*P*r/(1+Q)+math.sqrt(0.5*A**2*(1+1/Q)-(1-E2)*P*(z**2)/Q/(1+Q)-0.5*P*r**2)
        U=math.sqrt(z**2+(r-ro*E2)**2)
        V=math.sqrt((1-E2)*z**2+(r-ro*E2)**2)
        Zo=(B**2*z)/(A*V)
        
        ALT=U*(1-B**2/(A*V))
        LAT=math.atan2((z+EP**2*Zo),r)
        LONG=math.atan2(y,x)
        latmap.append(LAT*180/pi)
        longmap.append(LONG*180/pi)
        altmap.append(ALT)
        
        #increase time
        t = t + incTime
        j=j+1
    coordinates=[]
    coordinates.append(xmap)
    coordinates.append(ymap)
    coordinates.append(zmap)
    coordinates.append(latmap)
    coordinates.append(longmap)
    coordinates.append(altmap)
    return coordinates


#define function to calculate cross product 
def cross(a, b):
    result = [a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]]

    return result

#Normalize Vector
def norm(a):
    module = math.sqrt(a[0]**2+a[1]**2+a[2]**2)
    V=[]
    V.append(a[0]/module)
    V.append(a[1]/module)
    V.append(a[2]/module)
    return V

#Product of vectors by value
def prod(a,val):
    products = []
    for num1 in a:
        products.append(num1*val)
    return products

#Division of vector by value
def div(a,val):
    products = []
    for num1 in a:
        products.append(num1/val)
    return products

#Sum of vectors
def sum(a,b):
    result = []
    for num1,num2 in zip(a,b):
        result.append(num1+num2)
    return result

#Substraction of vectors
def sub(a,b):
    result = []
    for num1,num2 in zip(a,b):
        result.append(num1-num2)
    return result

#Modulo of a vector
def modulo(a):
    result=math.sqrt(a[0]**2+a[1]**2+a[2]**2)
    return result

#Dot product
def dot(a,b):
    result=0.0
    for n in range(len(a)):
        result=result+a[n]*b[n]
    return result

def visibilidadObs(a,latObs,longObs,altObs):
    #Visibilidad observador
    pi=math.pi
    Rearth=6.378e+06
    r=a-Rearth
    
    latvis=[]
    longvis=[]
    k=1
    j=1
    for i in range(-720,720):
        #Radio latitud
        latV=i*pi/1440
        latV=latV*180/pi
        satECEF=LLA2ECEF(latV, longObs, r)
        userECEF=LLA2ECEF(latObs,longObs,altObs)
        pointV=sub(satECEF,userECEF)
        NED=ECEF2NED(pointV,latObs*pi/180,longObs*pi/180)
        d=math.sqrt(NED[0]**2+NED[1]**2+NED[2]**2)
        alpha= math.asin((-NED[2])/d)*180/pi
        if(alpha>=10):
            latvis.append(latV)
            k=k+1
        
        #Radio longitud
        longV=i*pi/720
        longV=longV*180/pi
        satECEF=LLA2ECEF(latObs, longV, r)
        userECEF=LLA2ECEF(latObs,longObs,altObs)
        pointV=sub(satECEF,userECEF)
        NED=ECEF2NED(pointV,latObs*pi/180,longObs*pi/180)
        d=math.sqrt(NED[0]**2+NED[1]**2+NED[2]**2)
        alpha= math.asin((-NED[2])/d)*180/pi

        if(alpha>=10):
            longvis.append(longV)
            j=j+1
    #Generació elipse visibilitat
    a=(abs(longvis[len(longvis)-1]-longvis[0])/2)*0.925
    b=abs(latvis[len(latvis)-1]-latvis[0])/2
    x0 = longObs
    y0 = latObs
    Xvis=[]
    Yvis=[]
    for n in range(-100,101):
        t=n*0.01*pi
        Xvis.append(x0 + a*math.cos(t))
        Yvis.append(y0 + b*math.sin(t))
    visCoord=[Xvis,Yvis]
    return visCoord

def jsonCheck(datosObtenidos):
    try:
        return datosObtenidos.json()
    except ValueError:
        return("error")

if __name__ == '__main__':
    app.run()
 