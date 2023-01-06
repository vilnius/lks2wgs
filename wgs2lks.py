import math


def deg2rad(degrees):
    return degrees * (math.pi / 180)


def rad2deg(radians):
    return radians * (180 / math.pi)


def roundoff(x, y):
    return round(x * math.pow(10, y)) / math.pow(10, y)


# WGS-84 to LKS-94
def geo2grid(lat, lon):
    distsize = 3
    latrad = deg2rad(lat)

    k = 0.9998
    a = 6378137
    f = 1 / 298.257223563
    b = a * (1 - f)
    e2 = (a * a - b * b) / (a * a)

    w = deg2rad(lon - 24)
    t = math.tan(latrad)
    rho = a * (1 - e2) / math.pow(1 - (e2 * math.sin(latrad) * math.sin(latrad)), (3 / 2))
    nu = a / math.sqrt(1 - (e2 * math.sin(latrad) * math.sin(latrad)))

    psi = nu / rho
    coslat = math.cos(latrad)
    sinlat = math.sin(latrad)

    A0 = 1 - (e2 / 4) - (3 * e2 * e2 / 64) - (5 * math.pow(e2, 3) / 256)
    A2 = (3 / 8) * (e2 + (e2 * e2 / 4) + (15 * math.pow(e2, 3) / 128))
    A4 = (15 / 256) * (e2 * e2 + (3 * math.pow(e2, 3) / 4))
    A6 = 35 * math.pow(e2, 3) / 3072
    m = a * ((A0 * latrad) - (A2 * math.sin(2 * latrad)) + (A4 * math.sin(4 * latrad)) - (A6 * math.sin(6 * latrad)))

    eterm1 = (w * w / 6) * coslat * coslat * (psi - t * t)
    eterm2 = (math.pow(w, 4) / 120) * math.pow(coslat, 4) * (4 * math.pow(psi, 3) * (1 - 6 * t * t) + psi * psi * (1 + 8 * t * t) - psi * 2 * t * t + math.pow(t, 4))
    eterm3 = (math.pow(w, 6) / 5040) * math.pow(coslat, 6) * (61 - 479 * t * t + 179 * math.pow(t, 4) - math.pow(t, 6))
    dE = k * nu * w * coslat * (1 + eterm1 + eterm2 + eterm3)
    east = roundoff(500000 + dE, distsize)

    nterm1 = (w * w / 2) * nu * sinlat * coslat
    nterm2 = (math.pow(w, 4) / 24) * nu * sinlat * math.pow(coslat, 3) * (4 * psi * psi + psi - t * t)
    nterm3 = (math.pow(w, 6) / 720) * nu * sinlat * math.pow(coslat, 5) * (8 * math.pow(psi, 4) * (11 - 24 * t * t) - 28 * math.pow(psi, 3) * (1 - 6 * t * t) + psi * psi * (1 - 32 * t * t) - psi * 2 * t * t + math.pow(t, 4))
    nterm4 = (math.pow(w, 8) / 40320) * nu * sinlat * math.pow(coslat, 7) * (1385 - 3111 * t * t + 543 * math.pow(t, 4) - math.pow(t, 6))
    dN = k * (m + nterm1 + nterm2 + nterm3 + nterm4)
    north = roundoff(0 + dN, distsize)
    return north, east
 

# LKS-94 to WGS-84
def grid2geo(x, y):
    k = 0.9998
    a = 6378137
    f = 1 / 298.257223563
    b = a * (1 - f)
    e2 = (a * a - b * b) / (a * a)
    n = (a - b) / (a + b)
    G = a * (1 - n) * (1 - n * n) * (1 + (9 / 4) * n * n + (255 / 64) * math.pow(n, 4)) * (math.pi / 180)
    north = (y - 0)
    east = (x - 500000)
    m = north / k
    sigma = (m * math.pi) / (180 * G)
    footlat = sigma + ((3 * n / 2) - (27 * math.pow(n, 3) / 32)) * math.sin(2 * sigma) + ((21 * n * n / 16) - (55 * math.pow(n, 4) / 32)) * math.sin(4 * sigma) + (151 * math.pow(n, 3) / 96) * math.sin(6 * sigma) + (1097 * math.pow(n, 4) / 512) * math.sin(8 * sigma)
    rho = a * (1 - e2) / math.pow(1 - (e2 * math.sin(footlat) * math.sin(footlat)), (3 / 2))
    nu = a / math.sqrt(1 - (e2 * math.sin(footlat) * math.sin(footlat)))
    psi = nu / rho
    t = math.tan(footlat)
    x = east / (k * nu)
    laterm1 = (t / (k * rho)) * (east * x / 2)
    laterm2 = (t / (k * rho)) * (east * math.pow(x, 3) / 24) * (-4 * psi * psi + 9 * psi * (1 - t * t) + 12 * t * t)
    laterm3 = (t / (k * rho)) * (east * math.pow(x, 5) / 720) * (8 * math.pow(psi, 4) * (11 - 24 * t * t) - 12 * math.pow(psi, 3) * (21 - 71 * t * t) + 15 * psi * psi * (15 - 98 * t * t + 15 * math.pow(t, 4)) + 180 * psi * (5 * t * t - 3 * math.pow(t, 4)) + 360 * math.pow(t, 4))
    laterm4 = (t / (k * rho)) * (east * math.pow(x, 7) / 40320) * (1385 + 3633 * t * t + 4095 * math.pow(t, 4) + 1575 * math.pow(t, 6))
    latrad = footlat - laterm1 + laterm2 - laterm3 + laterm4

    lat_deg = rad2deg(latrad)

    seclat = 1 / math.cos(footlat)
    loterm1 = x * seclat
    loterm2 = (math.pow(x, 3) / 6) * seclat * (psi + 2 * t * t)
    loterm3 = (math.pow(x, 5) / 120) * seclat * (-4 * math.pow(psi, 3) * (1 - 6 * t * t) + psi * psi * (9 - 68 * t * t) + 72 * psi * t * t + 24 * math.pow(t, 4))
    loterm4 = (math.pow(x, 7) / 5040) * seclat * (61 + 662 * t * t + 1320 * math.pow(t, 4) + 720 * math.pow(t, 6))
    w = loterm1 - loterm2 + loterm3 - loterm4
    longrad = deg2rad(24) + w

    lon_deg = rad2deg(longrad)
    return lat_deg, lon_deg
