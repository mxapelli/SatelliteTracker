import os
from calc import *
from coord import *
from visibility import *
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
from datetime import datetime,timedelta, timezone
from apscheduler.schedulers.background import BackgroundScheduler
import numpy as np
from satellite_position import *


import time
import math
import numpy
import requests





load_dotenv()  # use dotenv to hide sensitive credential as environment variables


DATABASE_URL = f'mongodb+srv://MXapelli:{os.environ.get("passwordDB")}'\
    '@mongo-heroku-cluster.aiiqhhv.mongodb.net/satellites?retryWrites=true&w=majority'  # Get connection url from environment

# Establish connection with database
client = pymongo.MongoClient(DATABASE_URL)
mongo_db = client.db

app = Flask(__name__)

#
app.secret_key = os.environ.get("passwordSess")

#Global variables
pi = math.pi
Ge = 6.67384*10**(-11)  # Gravitational constant
Me = 5.972*10**24 # Mass of earth

@app.route('/', methods=['POST', 'GET'])
def index():
    # POST method to get coordinates from user
    if request.method == 'POST':
        coordsData = request.form['nm']
        print("Coordenadas ", coordsData)
        if isinstance(coordsData, str):
            if "/" in coordsData:
                coords = coordsData.split("/")
                urlAlt = "https://api.open-elevation.com/api/v1/lookup?locations=" + coords[0]+","+coords[1]
                result = requests.get(urlAlt)
                alt = result.json()
                altM = alt["results"][0]['elevation']
                session['latUser'] = float(coords[0])
                session['longUser'] = float(coords[1])
                session['altUser'] = float(altM)
            return redirect(url_for('index'))

    jsonSat = []
    for sat in mongo_db.satellites.find():
        jsonSat.append(
            {"OBJECT_NAME": str(sat["name"]), "NORAD_CAT_ID": str(sat["noradID"])})

    return render_template('index.html', datos=jsonSat)


@app.route('/<int:catnr>')
def satellite(catnr):
    start_time=time.time()
    xmap = []
    ymap = []
    zmap = []
    incTime = 5

    #Obtaining LLA coordinates from user session
    if 'latUser' in session:
        latUser = session['latUser']
        longUser = session['longUser']
        altUser = session['altUser']
    else:
        # Default coordinates of EETAC
        latUser = 41.2757116961354 
        longUser = 1.9872286269285848
        altUser = 4

    sat = mongo_db.satellites.find_one({"noradID": catnr})

    result= requests.get('https://celestrak.com/NORAD/elements/gp.php?CATNR='+str(catnr)+'&FORMAT=2le')
    print(result)

    tle_data = result.text.split('\n')
    print(tle_data)
    line1 = tle_data[0].strip()  # 1st TLE line
    line2 = tle_data[1].strip()
    print(line1)
    print(line2)
    
    lon_pos=[]
    lat_pos=[]

    positions=compute_satellite_positions(line1,line2)
    # Extracting latitude and longitude values
    lat_pos = [entry[1] for entry in positions]
    lon_pos = [entry[2] for entry in positions]

    if sat is None: # Handling error of non-existing satellite
        print("Error satellite does not exist")
        return render_template('error.html',error="Error 501",title="Error Satellite Not Found",text="We are sorry, there aren't any satellites with "+ str(catnr)+" as a Norad ID Number. We suggest that you try again with another ID number.")

    #Getting data of satellite
    name = sat['name']
    print("The user has selected satellite: ", name)
    epoch = sat['epoch']
    incl = sat['incl']
    omega = sat['omega']
    ecc = sat['ecc']
    w = sat['w']
    M = sat['M']
    n = sat['n']
    # Computing time to ToA
    ToA = time2toa(epoch)
    result = esecGAST(ToA)
    esec = result[0]
    GASTdegToA = result[1]

    # Process satellite kep. elements
    T = 86400 / n  # orbit period (secs.)
    n = 2 * pi / T  # mean movement [rad/s]
    a = (Ge * Me / n**2)**(1/3)  # semi-major axis [m]
    io = incl * pi / 180  # orbit inclination [rad]
    ecc = ecc  # orbit eccentricity
    # GAST at a given epoch is the RA of the Greenwich meridian at that epoch
    # Longitude of AN at ToA [rad] Aqui hay que poner GASTdegToA
    Oo = (omega - GASTdegToA) * pi / 180
    dO = 0  # Rate of Right ascension
    w = w * pi / 180  # argument of perigee [rad]
    Mo = M * pi / 180  # Mean anomaly at ToA [rad]

    # Groundtrack for period T in seconds
    coordinates = computeCoordinates(esec, T, n, a, io, ecc, Oo, dO, w, Mo, incTime)
    xmap = coordinates[0]
    ymap = coordinates[1]
    zmap = coordinates[2]
    latmap = coordinates[3]
    longmap = coordinates[4]
    altmap = coordinates[5]
    

    # Check Visibility of satellite
    elev = []
    vis = []
    visSat=[]
    for i in range(len(xmap)):
        satECEF = [xmap[i], ymap[i], zmap[i]]
        userECEF = LLA2ECEF(latUser, longUser, altUser)
        pointV = sub(satECEF, userECEF)
        NED = np.array(ECEF2NED(pointV, latUser*np.pi/180, longUser*np.pi/180))
        d = np.linalg.norm(NED)
        alpha = np.arcsin((-NED[2])/d) * 180/np.pi
        if (alpha >= 10):
            elev.append(alpha)
            visSat.append(1)
        else:
            visSat.append(0)
    if len(elev) > 0:
        vis = 1
    else:
        vis = 0

    # Visibility Area of User
    visPos = visibilidadObs(a, latUser, longUser, altUser, name)

    print("Time to process",name,time.time()-start_time)
    print(latmap[0],longmap[0])
    atime = time.localtime()
    st = time.strftime("%a, %d %b %Y %H:%M:%S ", atime)
    ##now = datetime.now()
    now = datetime.datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    session['time'] = dt_string
    print("session info", dt_string)
    text = [("You have selected the "+name+" satellite with Catalog Number: " + str(catnr))]
    return render_template('satellite.html', entries=text, longs=lon_pos, lats=lat_pos, number=str(catnr), userlat=latUser, userlong=longUser, vis=vis, visSat=visSat, satName=name,visPos=json.dumps(visPos))

@app.route('/<constellation_name>')
def constellation(constellation_name):
    start_time=time.time()

    constName = constellation_name.upper()

    print("Selected constellation:",constName)

    #Getting LLA Coordinates from user session
    if 'latUser' in session:
        latUser = session['latUser']
        longUser = session['longUser']
        altUser = session['altUser']
    else:
        latUser = 41.2757116961354
        longUser = 1.9872286269285848
        altUser = 4

    const="^"+constName
    listSat = []
    myquery = {"name": {"$regex": const}}
    for sat in mongo_db.satellites.find(myquery):
        listSat.append({"name": sat["name"], "catnr": sat["noradID"], "epoch": sat["epoch"], "ecc": sat["ecc"], "incl": sat["incl"], "omega": sat["omega"], "w": sat["w"], "M": sat["M"],
                        "n": sat["n"]})
    print("The constellation contains:",len(listSat),"satellites.")
    
    
    
    if len(listSat) == 0: #Handling error of non-existing constellation
        print("Error constellation does not exist")
        return render_template('error.html',error="Error 501",title="Error Constellation Not Found",text="We are sorry, the selected satellite constellation could not be processed as some information is missing. We suggest that you try a different one.")
    else:
        latmap = []
        longmap = []
        satname=[]
        satID=[]
        vis = []
        amax=0
        for i in range(len(listSat)):
            #Getting data of every satellite
            name = listSat[i]['name']
            satname.append(name)
            catnr = listSat[i]['catnr']
            satID.append(catnr)
            epoch = listSat[i]['epoch']
            incl = listSat[i]['incl']
            omega = listSat[i]['omega']
            ecc = listSat[i]['ecc']
            w = listSat[i]['w']
            M = listSat[i]['M']
            n = listSat[i]['n']

            #Getting time to ToA and GAST degrees at ToA
            ToA = time2toa(epoch)
            result = esecGAST(ToA)
            esec = result[0]
            GASTdegToA = result[1]
            
            # Process satellite kep. elements
            T = 86400 / n  # orbit period (secs.)
            n = 2 * pi / T  # mean movement [rad/s]
            a = (Ge * Me / n**2)**(1/3)  # semi-major axis [m]
            io = incl * pi / 180  # orbit inclination [rad]
            ecc = ecc  # orbit eccentricity
            # GAST at a given epoch is the RA of the Greenwich meridian at that epoch
            # Longitude of AN at ToA [rad] Aqui hay que poner GASTdegToA
            Oo = (omega - GASTdegToA) * pi / 180
            dO = 0  # Rate of Right ascension
            w = w * pi / 180  # argument of perigee [rad]
            Mo = M * pi / 180  # Mean anomaly at ToA [rad]

            t=esec
            if (constName=="GPS"):
                T=3600
                incTime=5
            if (constName=="STARLINK"):
                T=600
                incTime=10
            if (constName=="IRIDIUM"):
                T=600
                incTime=5
            if (constName=="INMARSAT"):
                T=7200
                incTime=5
            visSat=[]
            latSat=[]
            longSat=[]
            # Groundtrack for NT periods
            while t < esec + T:
                ecef= Kepler2ECEF(t, a, io, ecc, Oo, dO, w, Mo, n)
                # Check Visibility of satellite
                satECEF = [ecef[0], ecef[1], ecef[2]]
                userECEF = LLA2ECEF(latUser, longUser, altUser)
                pointV = sub(satECEF, userECEF)
                NED = ECEF2NED(pointV, latUser*pi/180, longUser*pi/180)
                d = math.sqrt(NED[0]**2+NED[1]**2+NED[2]**2)
                alpha = math.asin((-NED[2])/d)*180/pi
                if (alpha >= 10):
                    visSat.append(1)
                else:
                    visSat.append(0)
                lla=ECEF2LLA(ecef[0],ecef[1],ecef[2])
                lla[0] = round(lla[0], 10)
                lla[1] = round(lla[1], 10)
                latSat.append(lla[0])
                longSat.append(lla[1])
                t=t+incTime
            vis.append(visSat)
            longmap.append(longSat)
            latmap.append(latSat)
            
            if (a>amax):
                amax=a
            if((time.time()-start_time)>29):
                satmessage=", but just "+str(i)+" are represented.";
                break
            else:
                satmessage="."
        

        # Visibility Area of User
        visPos = visibilidadObs(amax, latUser, longUser, altUser, name)
   
        print("Time to process",constName,time.time()-start_time)
        atime = time.localtime()
        st = time.strftime("%a, %d %b %Y %H:%M:%S ", atime)
        text = [("You have selected the "+constName+" constellation, which contains "+str(len(listSat))+ " satellites"+satmessage)]
        return render_template('constellation.html', entries=text, longs=json.dumps(longmap), lats=json.dumps(latmap), satname=satname, userlat=latUser, userlong=longUser,constName=constName,vis=json.dumps(vis),satID=satID,visPos=json.dumps(visPos))


@app.route('/<int:catnr>/satdata')
def doppler(catnr):

    #Getting LLA coordinates from user session
    if 'latUser' in session:
        latUser = session['latUser']
        longUser = session['longUser']
        altUser = session['altUser']
    else:
        latUser = 41.2757116961354
        longUser = 1.9872286269285848
        altUser = 4
    
    #Checking if satellite exists
    sat = mongo_db.satellites.find_one({"noradID": catnr})
    

    if sat is None: # Handling non-existing satellite
        return render_template('error.html',error="Error 501",title="Error Satellite Not Found",text="We are sorry, there aren't any satellites with "+ str(catnr)+" as a Norad ID Number. We suggest that you try again with another ID number.")

    name = sat["name"]
    freq = sat["freq"]
    if freq == 1575.42*10**6:
        incTime = 30
    else:
        incTime = 10

    print("The user has selected visibility of satellite: ", name)
    epoch = sat['epoch']
    incl = sat['incl']
    omega = sat['omega']
    ecc = sat['ecc']
    w = sat['w']
    M = sat['M']
    n = sat['n']

    #Getting time to ToA and
    ToA=time2toa(epoch)
    result = esecGAST(ToA)
    esec = result[0]
    GASTdegToA = result[1]

    # Process satellite kep. elements
    T = 86400 / n  # orbit period (secs.)
    n = 2 * pi / T  # mean movement [rad/s]
    a = (Ge * Me / n**2)**(1/3)  # semi-major axis [m]
    io = incl * pi / 180  # orbit inclination [rad]
    ecc = ecc  # orbit eccentricity
    # GAST at a given epoch is the RA of the Greenwich meridian at that epoch
    # Longitude of AN at ToA [rad] Aqui hay que poner GASTdegToA
    Oo = (omega - GASTdegToA) * pi / 180
    dO = 0  # Rate of Right ascension
    w = w * pi / 180  # argument of perigee [rad]
    Mo = M * pi / 180  # Mean anomaly at ToA [rad]

    # Groundtrack for period T in seconds
    coordinates = computeCoordinates(esec, T, n, a, io, ecc, Oo, dO, w, Mo, incTime)
    xmap = coordinates[0]
    ymap = coordinates[1]
    zmap = coordinates[2]

    r = (1-ecc)*a
    v = math.sqrt(Ge*Me*((2/r)-(1/a)))
    Lm = r*v
    GAST1 = GAST(0)
    t4 = round(len(xmap)/4)
    GAST2 = GAST(incTime*t4)
    ECEF1 = ECEF2ECI(xmap[0], ymap[0], zmap[0], GAST1*math.pi/180)
    ECEF2 = ECEF2ECI(xmap[t4], ymap[t4], zmap[t4], GAST2*math.pi/180)

    # Orbital Plane Vector
    C = cross(ECEF1, ECEF2)
    V = norm(C)
    LmV = prod(V, Lm)

    # Doppler
    vSatV = []
    ECI = []
    GastV = []
    for n in range(len(xmap)):
        t = n*incTime
        GAST1 = GAST(t)
        GastV.append(GAST1)
        rV = ECEF2ECI(xmap[n], ymap[n], zmap[n], GAST1*math.pi/180)
        ECI.append(rV)
        cros = cross(LmV, rV)
        absrV = math.sqrt(rV[0]**2+rV[1]**2+rV[2]**2)**2
        vSat = div(cros, absrV)
        vSatV.append(vSat)
        vSatAbs = math.sqrt(vSat[0]**2+vSat[1]**2+vSat[2]**2)

    # Vector Velocidad Observador
    Rearth = 6378*1000
    Tearth = (23*60+56)*60+4  # segons
    vearth = 2*math.pi*(Rearth/Tearth)
    latObs = latUser
    longObs = longUser
    altObs = altUser

    dSatObs = []
    dObs = []
    vinst = []
    vinstReal = []
    fD = []
    fDReal = []
    for n in range(len(ECI)):

        # Velocidad Observador en Funcion de la Latitud
        Gast = GastV[n]
        GaLo = (Gast+longObs)*math.pi/180
        vObs = [vearth*-math.sin(GaLo), vearth*math.cos(GaLo), 0]
        vObsLat = prod(vObs, math.cos(latObs*math.pi/180))

        # Velocidad Satelite respecto Observador
        vResta = sub(vSatV[n], vObsLat)

        # Posicion Satelite respecto Observador
        userECEF = LLA2ECEF(latObs, longObs, altObs)
        userECI = ECEF2ECI(userECEF[0], userECEF[1],
                           userECEF[2], Gast*math.pi/180)
        satV = sub(userECI, ECI[n])
        # print(satV)
        dSatObs.append(modulo(satV))
        dObs.append(modulo(userECEF))
        satVnorm = norm(satV)
        vinst.append(dot(vSatV[n], satVnorm))
        vinstReal.append(dot(vResta, satVnorm))
        fD.append(freq*vinst[n]/(3*10**8))
        fDReal.append((freq*vinstReal[n])/(3*10**8))

    id = {"noradID": catnr}
    newvalues = {"$set": {"vDoppler": fDReal}}
    mongo_db.satellites.update_one(id, newvalues)

    # Visibility of satellite
    elev = []
    az = []
    ind = []
    for i in range(len(xmap)):
        satECEF = [xmap[i], ymap[i], zmap[i]]
        userECEF = LLA2ECEF(latUser, longUser, altUser)
        pointV = sub(satECEF, userECEF)
        NED = ECEF2NED(pointV, latUser*pi/180, longUser*pi/180)
        d = math.sqrt(NED[0]**2+NED[1]**2+NED[2]**2)
        alpha = math.asin((-NED[2])/d)*180/pi
        beta = math.atan2(NED[1], NED[0])*180/pi
        if (alpha >= 10):
            elev.append(alpha)
            ind.append(i)
            if beta < 0:
                az.append((beta+360)*pi/(180))
            else:
                az.append(beta*pi/180)

    dopplerVis = []
    timeP = []
    cat = []
    c = 0
    catN = 0
    change = []
    pointChange = 0
    for i in range(len(ind)):
        dopplerVis.append(fDReal[ind[i]])
        timeP.append(ind[i]*incTime)
        if (c != ind[i]-1 and c != 0):
            catN = 1
            # To detect visibility changes
            change.append(ind[i-1]*incTime)
            change.append(ind[i]*incTime)
            pointChange = i
        cat.append(catN)
        c = ind[i]

    now = datetime.now(timezone.utc)
    actual_time = now.strftime("%d/%m/%Y %H:%M:%S")
    if len(timeP) < 1:
        return render_template('error.html',error="Error Not Visible",title="Satellite Not Visible",text="We are sorry, the selected satellite "+name+" is not visible during his next orbit. We suggest that you try later or select a different satellite.")
    init = timeP[0]
    fin = timeP[-1]
    duration = (fin-init)/7
    print(actual_time)

    #Creating the texts for the ticks of Doppler graph
    actual_time_obj = datetime.strptime(actual_time, '%d/%m/%Y %H:%M:%S')+timedelta(hours=2)# It convert the strings to datetime obj
    times = []
    days = []
    for i in range(8):
        t = init+duration*(i)
        time_obj = actual_time_obj + timedelta(seconds=t)
        timeText = str(time_obj).split()
        timev = timeText[1].split('.')
        times.append(timev[0])
        if (i == 0):
            days.append(time_obj.strftime("%d %b"))

    visStep = []
    visStep.append(times[0])
    if len(change) > 0:
        fin1 = actual_time_obj + timedelta(seconds=change[0])
        days.append(fin1.strftime("%d %b"))
        in2 = actual_time_obj+timedelta(seconds=change[1])
        days.append(in2.strftime("%d %b"))
        timeText = str(fin1).split()
        timev = timeText[1].split('.')
        visStep.append(timev[0])
        timeText = str(in2).split()
        timev = timeText[1].split('.')
        visStep.append(timev[0])
    visStep.append(times[-1])
    fin2 = actual_time_obj + timedelta(seconds=fin)
    days.append(fin2.strftime("%d %b"))
    visText = "The satellite will be visible from " + \
        visStep[0]+" ("+days[0]+") "+" to " + visStep[1]+" ("+days[1]+") "
    if len(visStep) > 2:
        visText = visText + " and from " + \
            visStep[2]+" ("+days[2]+") "+"to " + \
            visStep[3]+" ("+days[3]+")"+"."

    # New Polar Chart Visibility
    img = BytesIO()
    colormap = numpy.array(['yellow', 'red'])

    fig = plt.figure(facecolor='#383A3F')
    ax = fig.add_subplot(projection='polar')
    ax.grid(c='white')
    d = ax.scatter(az, elev, c=colormap[cat], s=30, alpha=1)
    if len(change) > 0:
        azfirst = []
        elevfirst = []
        azfirst.append(az[0])
        azfirst.append(az[pointChange])
        elevfirst.append(elev[0])
        elevfirst.append(elev[pointChange])
        e = ax.scatter(azfirst, elevfirst, c='black', s=40, alpha=1)
    else:
        e = ax.scatter(az[0], elev[0], c='black', s=40, alpha=1)

    ax.set_title("Visibility of Satellite "+name+" with ID " +
                 str(catnr), va='bottom', c='white')
    ax.set_theta_zero_location('N')
    ax.set_rlabel_position(-90)
    ax.set_theta_direction(-1)  # theta increasing clockwise
    ax.set_xticklabels(['N', '45', 'E', '135', 'S', '225',
                       'W', '315'], c='white', weight='bold')
    ax.set_ylim(90, 0)
    ax.set_yticklabels(['0', '20', '40', '60', '80'],
                       color='white', weight='bold')
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_facecolor((0, 0, 1))
    ax.set_axisbelow(True)

    #Legend of plot
    pop_a = mpatches.Patch(color='yellow', label='First Stretch')
    pop_b = mpatches.Patch(color='red', label='Second Stretch')
    pop_c = mpatches.Patch(color='black', label='First point of Visibility')
    if len(change) > 0:
        fig.legend(handles=[pop_a, pop_b, pop_c], loc="lower center", ncol=3)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.12)
    else:
        fig.legend(handles=[pop_a, pop_c], loc="lower center", ncol=2)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.12)

    #Saving plot and generating url
    plt.savefig(img, format='png', bbox_inches="tight")
    plt.close()
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf8')

    # Doppler new graph
    img2 = BytesIO()

    fig, ax = plt.subplots(facecolor='#383A3F')
    ax.grid(c='white')
    k = ax.scatter(timeP, dopplerVis, c=colormap[cat], s=40, alpha=1)
    ax.set_title("Doppler frequency of Satellite "+name +
                 " with ID "+str(catnr), va='bottom', c='white')
    ax.set_xlabel('Time (UTC)', c='white')  # Add an x-label to the axes.
    xticks = numpy.arange(init, fin+1, step=duration)
    ax.set_xticks(xticks, c='white')
    ax.set_xticklabels(times, c='white', rotation=90)

    # Setting ticks for y-axis
    yticks = ax.get_yticks()
    yticklabel = []
    for i in range(len(yticks)):
        yticklabel.append(yticks[i]/1000)

    # Add a y-label to the axes
    ax.set_ylabel('Doppler Frequency (kHz)', c='white')
    ax.set_yticklabels(yticklabel, color='white')
    ax.set_facecolor((0, 0, 1))
    ax.set_axisbelow(True)

    #Saving plot and generating url
    plt.savefig(img2, format='png', bbox_inches="tight")
    plt.close()
    img2.seek(0)
    plotDoppler_url = base64.b64encode(img2.getvalue()).decode('utf8')

    init = 0
    fin = (len(vinstReal)-1)*incTime
    duration = (fin-init)/7
    timeT = list(range(init, fin+1,incTime))

    #Creating the texts for the ticks of Doppler graph
    # It convert the strings to datetime obj
    times2 = []
    for i in range(8):
        t = init+duration*(i)
        time_obj = actual_time_obj + timedelta(seconds=t)
        timeText = str(time_obj).split()
        timev = timeText[1].split('.')
        times2.append(timev[0])

    # New Graph Distance to Satellite
    # New plot3
    img3 = BytesIO()

    fig, ax = plt.subplots(facecolor='#383A3F')
    ax.grid(c='white')
    k = ax.scatter(timeT, dSatObs, c='magenta', s=40, alpha=1)
    ax.set_title("Distance to satellite "+name +
                 " with ID "+str(catnr), va='bottom', c='white')
    ax.set_xlabel('Time (UTC)', c='white')  # Add an x-label to the axes.
    xticks = numpy.arange(init, fin+1, step=duration)
    ax.set_xticks(xticks, c='white')
    ax.set_xticklabels(times2, c='white', rotation=90)

    # Setting ticks for y-axis
    yticks = ax.get_yticks()
    yticklabel = []
    for i in range(len(yticks)):
        yticklabel.append(yticks[i]/1000)

    # Add a y-label to the axes
    ax.set_ylabel('Distance (km)', c='white')
    ax.set_yticklabels(yticklabel, color='white')
    ax.set_facecolor((0, 0, 1))
    ax.set_axisbelow(True)

    #Saving plot and generating url
    plt.savefig(img3, format='png', bbox_inches="tight")
    plt.close()
    img3.seek(0)
    plotDistance_url = base64.b64encode(img3.getvalue()).decode('utf8')



    # New plot Radial speed of Satellite
    # New plot4
    img4 = BytesIO()

    fig, ax = plt.subplots(facecolor='#383A3F')
    ax.grid(c='white')
    k = ax.scatter(timeT, vinstReal, c='aqua', s=40, alpha=1)
    ax.set_title("Radial speed of Satellite "+name +
                 " with ID "+str(catnr), va='bottom', c='white')
    ax.set_xlabel('Time (UTC)', c='white')  # Add an x-label to the axes.
    xticks = numpy.arange(init, fin+1, step=duration)
    ax.set_xticks(xticks, c='white')
    ax.set_xticklabels(times2, c='white', rotation=90)

    # Setting ticks for y-axis
    yticks = ax.get_yticks()
    yticklabel = []
    for i in range(len(yticks)):
        yticklabel.append(yticks[i])

    # Add a y-label to the axes
    ax.set_ylabel('Radial speed (m/s)', c='white')
    ax.set_yticklabels(yticklabel, color='white')
    ax.set_facecolor((0, 0, 1))
    ax.set_axisbelow(True)

    #Saving plot and generating url
    plt.savefig(img4, format='png', bbox_inches="tight")
    plt.close()
    img4.seek(0)
    plotRadialSpeed_url = base64.b64encode(img4.getvalue()).decode('utf8')

    text1 = [("The frequency for this satellite is: "+str(freq/10**6)+" MHz"), (visText)]
    text2 =[("Maximum Doppler frequency: "+str(round(max(dopplerVis)/1000, 2))+" kHz"),
            ("Minimum Doppler frequency: "+str(round(min(dopplerVis)/1000, 2))+" kHz"),]
    text3 = [("Maximum radial speed: " + str(round(max(vinstReal), 2)) + " m/s"),("Minimum radial speed: " + str(round(min(vinstReal), 2)) + " m/s")]
    text4 = [("Maximum distance to satellite: " + str(round(max(dSatObs)/1000, 2)) + " km"),("Minimum distance to satellite: " + str(round(min(dSatObs)/1000, 2)) + " km")]
    return render_template('satAnalysis.html', entries1=text1, entries2=text2,speed=text3,distance=text4,doppler=fDReal, plot_url=plot_url, plotDoppler_url=plotDoppler_url, plotDistance_url=plotDistance_url,plotRadialSpeed_url=plotRadialSpeed_url,satname=name)

@app.route('/<int:catnr>/satdata/<texto>')
def error1(catnr,texto):
    text2=str(catnr)+"/satdata/"+texto
    return render_template('error.html',error="Error 404",title="URL Not Found",text="We are sorry, the requested URL "+ text2 +" was not found on the server.")

@app.route('/<int:catnr>/<texto>')
def error2(catnr,texto):
    text2=str(catnr)+"/"+texto
    return render_template('error.html',error="Error 404",title="URL Not Found",text="We are sorry, the requested URL "+ text2 +" was not found on the server.")

@app.route('/<texto1>/<texto2>')
def error3(texto1,texto2):
    text2=texto1+"/"+texto2
    return render_template('error.html',error="Error 404",title="URL Not Found",text="We are sorry, the requested URL "+ text2 +" was not found on the server.")

# Function to update database with new satellite data
def dbUpdate():

    satsDB = []
    satsCelestrak = []
    start_time = time.time()
    sats = mongo_db.satellites.find() #Find all satellites

    for i in sats: #Add DB satellites to array
        satsDB.append(i)

    constellations = ['GROUP=starlink', 'GROUP=iridium-next',
                      'GROUP=gps-ops', 'CATNR=25544', 'GROUP=iridium', 'NAME=inmarsat','CATNR=39215']
    
    #Download satellites from Celestrak
    for i in range(len(constellations)):
        result_data = requests.get(
            'https://celestrak.com/NORAD/elements/gp.php?'+constellations[i]+'&FORMAT=JSON')
        check = jsonCheck(result_data)
        if check != "error":
            satsCelestrak.append(result_data.json())

    # Updating existing satellites
    for sat in satsDB:
        catnr = sat['noradID']
        name = sat['name']
        found = 0
        while found == 0:
            for i in range(len(satsCelestrak)):
                for j in range(len(satsCelestrak[i])):
                    noradID = int(satsCelestrak[i][j]['NORAD_CAT_ID'])
                    epoch = satsCelestrak[i][j]['EPOCH']
                    if (catnr == noradID):
                        if "IRIDIUM" in name:  # Esto es momentaneo, luego se borrara
                            name = name.strip("[-]")
                        if "STARLINK" in name:  # Esto es momentaneo, luego se borrara
                            name = satsCelestrak[i][j]['OBJECT_NAME']
                        if sat['epoch'] != epoch:
                            incl = satsCelestrak[i][j]['INCLINATION']
                            omega = satsCelestrak[i][j]['RA_OF_ASC_NODE']
                            ecc = satsCelestrak[i][j]['ECCENTRICITY']
                            w = satsCelestrak[i][j]['ARG_OF_PERICENTER']
                            M = satsCelestrak[i][j]['MEAN_ANOMALY']
                            n = satsCelestrak[i][j]['MEAN_MOTION']
                            sats = mongo_db.satellites.find({"noradID": noradID})
                            count = len(list(sats))
                            print(count)
                            if count>1:
                                mongo_db.satellites.delete_one({"noradID": noradID})
                            id = {"noradID": catnr}
                            newvalues = {"$set": {"name": name, "epoch": epoch, "incl": incl,
                                                  "ecc": ecc, "omega": omega, "w": w, "M": M, "n": n}}
                            mongo_db.satellites.update_one(id, newvalues)
                            print(name, " has been updated")
                        del satsCelestrak[i][j]
                        break
            found = 1
    print("Updated Completed")

    # Cleaning list of sats
    i = 0
    while i < (len(satsCelestrak)-1):
        if len(satsCelestrak[i]) == 0:
            del satsCelestrak[i]
            i = i-1
        i = i+1

    # Adding new satellites
    for i in range(len(satsCelestrak)):
        for j in range(len(satsCelestrak[i])):
            name = satsCelestrak[i][j]['OBJECT_NAME']
            catnr = int(satsCelestrak[i][j]['NORAD_CAT_ID'])
            epoch = satsCelestrak[i][j]['EPOCH']
            incl = satsCelestrak[i][j]['INCLINATION']
            omega = satsCelestrak[i][j]['RA_OF_ASC_NODE']
            ecc = satsCelestrak[i][j]['ECCENTRICITY']
            w = satsCelestrak[i][j]['ARG_OF_PERICENTER']
            M = satsCelestrak[i][j]['MEAN_ANOMALY']
            n = satsCelestrak[i][j]['MEAN_MOTION']

            fre = 0
            if "ISS" in name:
                fre = 145.8*10**6
            elif "STARLINK" in name:
                fre = 10950*10**6
            elif "NAVSTAR" in name:
                fre = 1575.42*10**6
            elif "IRIDIUM" in name:
                fre = 1626.1042*10**6
            elif "GPS" in name:
                fre = 1575.42*10**6
            elif "INMARSAT" or "ALPHA" in name: ##Corregir freq
                fre = 1575.42*10**6

            if fre != 0:
                sat = {"noradID": catnr, "name": name, "epoch": epoch, "ecc": ecc, "incl": incl, "omega": omega, "w": w, "M": M,
                       "n": n, "freq": fre, "xECEF": [], "yECEF": [], "zECEF": [], "vDoppler": []}
                mongo_db.satellites.insert_one(sat)
                print(name, " has been inserted in DB")
    print("Added Completed")
    now = datetime.now()
    t = time.time()-start_time
    print("Database updated at", now, "It took",t, "seconds to complete the task.")

# Scheduler to program db update every 10 minutes
scheduler = BackgroundScheduler()
job = scheduler.add_job(dbUpdate, 'interval', minutes=360)
scheduler.start()


# Function that checks if the json object is correct
def jsonCheck(datosObtenidos):
    try:
        return datosObtenidos.json()
    except ValueError:
        return ("error")
    
def json_to_tle(json_obj):
    # Extract values from the JSON object
    norad_cat_id = json_obj['NORAD_CAT_ID']
    classification = json_obj['CLASSIFICATION_TYPE']
    int_designator_year = json_obj['OBJECT_ID'].split('-')[0][-2:] # Last two digits of year
    int_designator_launch_number = json_obj['OBJECT_ID'].split('-')[1]
    epoch_year = json_obj['EPOCH'][:4][-2:]
    epoch_day = str((int(json_obj['EPOCH'][5:7]) - 1) * 30 + int(json_obj['EPOCH'][8:10]) + float(json_obj['EPOCH'][11:23])/86400)[:11]
    mean_motion_dot = json_obj['MEAN_MOTION_DOT']
    bstar = json_obj['BSTAR']
    elset_num = json_obj['ELEMENT_SET_NO']
    inclination = json_obj['INCLINATION']
    ra_of_asc_node = json_obj['RA_OF_ASC_NODE']
    eccentricity = str(json_obj['ECCENTRICITY'])[2:]
    arg_of_pericenter = json_obj['ARG_OF_PERICENTER']
    mean_anomaly = json_obj['MEAN_ANOMALY']
    mean_motion = json_obj['MEAN_MOTION']
    rev_at_epoch = json_obj['REV_AT_EPOCH']

    # Format the TLE strings
    line1 = f"1 {norad_cat_id:5}U {int_designator_year}{int_designator_launch_number} {epoch_year}{epoch_day} {mean_motion_dot: .8f} 00000+0 {bstar: .8f} 0 {elset_num:4}"
    line2 = f"2 {norad_cat_id:5} {inclination:7.4f} {ra_of_asc_node:7.4f} {eccentricity:7} {arg_of_pericenter:7.4f} {mean_anomaly:7.4f} {mean_motion:11.8f} {rev_at_epoch:5}"

    return line1, line2


if __name__ == '__main__':
    app.run()
