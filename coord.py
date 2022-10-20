import math
import numpy
import time

pi=math.pi

#Function that transforms ECEF coordinates to ECI coordinates
def ECEF2ECI(X, Y, Z, B):
    ECI = []
    Xe = math.cos(B)*X-math.sin(B)*Y
    Ye = math.sin(B)*X+math.cos(B)*Y
    Ze = Z
    ECI.append(Xe)
    ECI.append(Ye)
    ECI.append(Ze)
    return ECI

#Function that transforms LLA coordinates to ECEF coordinates
def LLA2ECEF(lat, lon, alt):
    a = 6378137
    e2 = 0.00669437999014
    e = math.sqrt(e2)
    lat = lat*math.pi/180
    lon = lon*math.pi/180

    coco = math.sqrt(1-e2*math.sin(lat)**2)
    clat = math.cos(lat)
    clong = math.cos(lon)
    slat = math.sin(lat)
    slong = math.sin(lon)
    x = (a/coco+alt)*clat*clong
    y = (a/coco+alt)*clat*slong
    z = (a*(1-e2)/(coco)+alt)*slat
    ecef = [x, y, z]
    return ecef

#Function that transforms ECEF coordinates to LLA coordinates
def ECEF2LLA(x,y,z):
    lla=[]
    # LLA coordinates [lat, lon, alt] = ECEF2LLA(x, y, z)
    # Convert from cartesian (x,y,z) to LLA coordinates in datum WGS84
    A = 6378137
    E2 = 0.00669437999014  # datos elipsoide WGS84
    B = A*math.sqrt(1-E2)
    EP = A*math.sqrt(E2)/B

    r = math.sqrt(x**2+y**2)
    EE2 = A**2-B**2
    F = 54*B**2*z**2
    G = r**2+(1-E2)*z**2-E2*EE2
    C = E2**2*(F*r**2)/G**3
    S = (1+C+math.sqrt(C**2+2*C))**(1/3)
    P = F/(3*G**2*(S+1/S+1)**2)
    Q = math.sqrt(1+2*E2**2*P)
    ro = -E2*P*r/(1+Q)+math.sqrt(0.5*A**2*(1+1/Q) - (1-E2)*P*(z**2)/Q/(1+Q)-0.5*P*r**2)
    U = math.sqrt(z**2+(r-ro*E2)**2)
    V = math.sqrt((1-E2)*z**2+(r-ro*E2)**2)
    Zo = (B**2*z)/(A*V)

    ALT = U*(1-B**2/(A*V))
    LAT = math.atan2((z+EP**2*Zo), r)
    LONG = math.atan2(y, x)
    lla.append(LAT*180/pi)
    lla.append(LONG*180/pi)
    lla.append(ALT)
    return lla

#Function that omputes ECEF coordinates from  Kepler orbital elements
def Kepler2ECEF(t, a, io, ecc, Oo, dO, w, Mo, n):
    # ECEF coordinates
    ecef=[]
    M = Mo + n * t  # Mean anomaly now [rad]
    # Solve Keppler equation for eccentric anomally iteratively
    # Loop is exited when we reach the desired tolerance
    tol = 10**(-8)  # Converging tolerance
    E = M
    while True:
        Eo = E
        E = M + ecc * math.sin(Eo)
        if abs(Eo - E) < tol:
            break

    # True anomaly
    cosv = (math.cos(E) - ecc) / (1 - ecc * math.cos(E))
    sinv = math.sqrt(1 - ecc**2) * math.sin(E) / (1 - ecc * math.cos(E))
    v = math.atan2(sinv, cosv)

    # Cartesian coordinates within the orbital plane
    u = v + w
    r = a * (1 - ecc * math.cos(E))
    # Corrected right ascension (7.2921151467e-5 = rate of Earth rotation [rad/s])
    O = Oo - 7.2921151467e-5 * t + dO * t
    xp = r * math.cos(u)
    yp = r * math.sin(u)

    # ECEF coordinates [m]
    x = xp * math.cos(O) - yp * math.cos(io) * math.sin(O)
    y = xp * math.sin(O) + yp * math.cos(io) * math.cos(O)
    z = yp * math.sin(io)

    ecef.append(x)
    ecef.append(y)
    ecef.append(z)

    return ecef

#Function that transforms ECEF coordinates to NED
def ECEF2NED(ECEF, phi, lamda):

    # 1) Rotation Matrix from NED to ECEF: ECEF = M*NED
    r = -math.sin(phi)*math.cos(lamda)
    M = numpy.array([[(-math.sin(phi)*math.cos(lamda)), (-math.sin(lamda)), (-math.cos(phi)*math.cos(lamda))], [(-math.sin(phi)
                    * math.sin(lamda)), (math.cos(lamda)), (-math.cos(phi)*math.sin(lamda))], [(math.cos(phi)), 0, (-math.sin(phi))]], dtype='f')
    # 2) Compute NED coordinates
    ECEFt = numpy.array([ECEF]).T
    invM = numpy.linalg.inv(M)

    NED = numpy.matmul(invM, ECEFt)
    return NED

#Function that transforms ECEF coordinates to NED
def NED2ECEF(NED, phi, lamda):

    # 1) Rotation Matrix from NED to ECEF: ECEF = M*NED
    M = numpy.array([[(-math.sin(phi)*math.cos(lamda)), (-math.sin(lamda)), (-math.cos(phi)*math.cos(lamda))], [(-math.sin(phi)
                    * math.sin(lamda)), (math.cos(lamda)), (-math.cos(phi)*math.sin(lamda))], [(math.cos(phi)), 0, (-math.sin(phi))]], dtype='f')
    # 2) Compute NED coordinates
    NEDt = numpy.array([NED]).T

    ECEFT = numpy.matmul(M, NEDt)
    ECEF=[]
    ECEF.append(ECEFT[0][0])
    ECEF.append(ECEFT[1][0])
    ECEF.append(ECEFT[2][0])
    return ECEF

def ELEV_AZ2NED(r, E, a):
    theta = E
    phi = a
    x, y, z= sph2cart(phi, theta, r);
    NED=[]
    NED.append(x)   # North
    NED.append(-y)  # East
    NED.append(-z) # Down
    
    return NED

def sph2cart(azimuth,elevation,r):
    x = r * math.cos(elevation) * math.cos(azimuth)
    y = r * math.cos(elevation) * math.sin(azimuth)
    z = r * math.sin(elevation)
    return x, y, z

#Function that computes the coordinates of a satellite
def computeCoordinates(esec, T, n, a, io, ecc, Oo, dO, w, Mo, incTime):
    t = esec
    NT = 1
    j = 0
    xmap = []
    ymap = []
    zmap = []
    latmap = []
    longmap = []
    altmap = []
    R = []
    while t < esec + T * NT:
        #Getting all the ECEF coordinates
        ecef=Kepler2ECEF(t, a, io, ecc, Oo, dO, w, Mo, n)

        xmap.append(ecef[0])
        ymap.append(ecef[1])
        zmap.append(ecef[2])
        
        # increase time
        t = t + incTime
    
    #Getting all LLA coordinates
    for i in range(len(xmap)):
        lla=ECEF2LLA(xmap[i],ymap[i],zmap[i])
        latmap.append(lla[0])
        longmap.append(lla[1])
        altmap.append(lla[2])

    coordinates = []
    coordinates.append(xmap)
    coordinates.append(ymap)
    coordinates.append(zmap)
    coordinates.append(latmap)
    coordinates.append(longmap)
    coordinates.append(altmap)
    return coordinates


    #Function that computes GAST [deg] at current time + esec [sec]
def GAST(esec):
    # esec = ToA time in secs. from/to now (esec < 0 means ToA is in the past)
    # Find Julian Date now
    unixepoch = 2440587.5  # JD at unix epoch: 0h (UTC) 1/1/1970
    unixseconds = time.time()
    JD = unixepoch + unixseconds / (86400)
    JToA = JD + esec / 86400  # Julian Date at ToA

    # Greenwich Mean Sidereal Time (GMST) is the hour angle of the average position of the vernal equinox,
    # neglecting short term motions of the equinox due to nutation. GAST is GMST corrected for
    # the shift in the position of the vernal equinox due to nutation.
    # GAST at a given epoch is the RA of the Greenwich meridian at that epoch (usually in time units).

    # Find GAST in degrees at ToA
    # epoch is 1/1/2000 at 12:00 UTC
    J2000 = 2451545.0
    midnight = round(JToA) - 0.5                   # midnight of JToA
    days_since_midnight = JToA - midnight
    hours_since_midnight = days_since_midnight * 24
    days_since_epoch = JToA - J2000
    centuries_since_epoch = days_since_epoch / 36525
    whole_days_since_epoch = midnight - J2000
    GAST = 6.697374558 + 0.06570982441908 * whole_days_since_epoch + 1.00273790935 * \
        hours_since_midnight + 0.000026 * centuries_since_epoch**2  # GAST in hours from ?
    GASTh = GAST % 24  # GAST in hours at ToA
    # GAST in degrees at ToA (approx. 361º/24h)
    GASTdeg = 15 * 1.0027855 * GASTh
    return GASTdeg

#Function that computes esec and GAST from a satelite
def esecGAST(ToA):
    unixepoch = 2440587.5  # JD at unix epoch: 0h (UTC) 1/1/1970
    unixseconds = time.time()
    JD = unixepoch + unixseconds / (86400)  # Julian Date now in UTC

    # Compute current time in secs. from ToA
    # MJD on 1/1/2022 at 00:00 UTC (see http://leapsecond.com/java/cal.htm)
    MJ2022 = 59580
    J2022 = MJ2022 + 2400000.5  # JD on 1/1/2022 at 00:00 UTC
    JToA = J2022 + ToA - 1  # ToA in JD
    esec = 86400 * (JD - JToA)  # Time since ToA in seconds

    # Find GAST in degrees at ToA
    J2000 = 2451545.0  # epoch is 1/1/2000 at 12:00 UTC
    midnight = round(JToA) - 0.5  # midnight of JToA
    days_since_midnight = JToA - midnight
    hours_since_midnight = days_since_midnight * 24
    days_since_epoch = JToA - J2000
    centuries_since_epoch = days_since_epoch / 36525
    whole_days_since_epoch = midnight - J2000
    GAST = 6.697374558 + 0.06570982441908 * whole_days_since_epoch + 1.00273790935 * \
        hours_since_midnight + 0.000026 * centuries_since_epoch**2  # GAST in hours from ?
    GASTh = GAST % 24  # GAST in hours at ToA

    # GAST in degrees at ToA (approx. 361º/24h)
    GASTdegToA = 15 * 1.0027855 * GASTh
    esecGAST = []
    esecGAST.append(esec)
    esecGAST.append(GASTdegToA)
    return esecGAST


#Function that computes the time to ToA (Time of applicability)
def time2toa(epoch):
    #Getting date and time from epoch
    epochT = str(epoch).split("T")
    epochT[1] = epochT[1].split(".")
    hour = epochT[1][0]
    date = epochT[0]

    # Computing time to ToA
    x = time.strptime(date+" " + hour, "%Y-%m-%d %H:%M:%S")
    s = x.tm_yday
    s = s+x.tm_hour/24+x.tm_min/(24*60)+x.tm_sec/(24*60*60)
    ToA = s

    return ToA

