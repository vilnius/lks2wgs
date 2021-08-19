const deg2rad = (degrees) => degrees * (Math.PI / 180)
const rad2deg = (radians) => radians * (180 / Math.PI)
const roundoff = (x, y) => Math.round(x * Math.pow(10, y)) / Math.pow(10, y)

// WGS-84 to LKS-94
function geo2grid(lat, lon) {
    const distsize = 3;
    const latrad = deg2rad(lat);

    const k = 0.9998;
    const a = 6378137;
    const f = 1 / 298.257223563;
    const b = a * (1 - f);
    const e2 = (a * a - b * b) / (a * a);

    const w = deg2rad(lon - 24);
    const t = Math.tan(latrad);
    const rho = a * (1 - e2) / Math.pow(1 - (e2 * Math.sin(latrad) * Math.sin(latrad)), (3 / 2));
    const nu = a / Math.sqrt(1 - (e2 * Math.sin(latrad) * Math.sin(latrad)));

    const psi = nu / rho;
    const coslat = Math.cos(latrad);
    const sinlat = Math.sin(latrad);

    const A0 = 1 - (e2 / 4) - (3 * e2 * e2 / 64) - (5 * Math.pow(e2, 3) / 256);
    const A2 = (3 / 8) * (e2 + (e2 * e2 / 4) + (15 * Math.pow(e2, 3) / 128));
    const A4 = (15 / 256) * (e2 * e2 + (3 * Math.pow(e2, 3) / 4));
    const A6 = 35 * Math.pow(e2, 3) / 3072;
    const m = a * ((A0 * latrad) - (A2 * Math.sin(2 * latrad)) + (A4 * Math.sin(4 * latrad)) - (A6 * Math.sin(6 * latrad)));

    const eterm1 = (w * w / 6) * coslat * coslat * (psi - t * t);
    const eterm2 = (Math.pow(w, 4) / 120) * Math.pow(coslat, 4) * (4 * Math.pow(psi, 3) * (1 - 6 * t * t) + psi * psi * (1 + 8 * t * t) - psi * 2 * t * t + Math.pow(t, 4));
    const eterm3 = (Math.pow(w, 6) / 5040) * Math.pow(coslat, 6) * (61 - 479 * t * t + 179 * Math.pow(t, 4) - Math.pow(t, 6));
    const dE = k * nu * w * coslat * (1 + eterm1 + eterm2 + eterm3);
    const east = roundoff(500000 + dE, distsize);

    const nterm1 = (w * w / 2) * nu * sinlat * coslat;
    const nterm2 = (Math.pow(w, 4) / 24) * nu * sinlat * Math.pow(coslat, 3) * (4 * psi * psi + psi - t * t);
    const nterm3 = (Math.pow(w, 6) / 720) * nu * sinlat * Math.pow(coslat, 5) * (8 * Math.pow(psi, 4) * (11 - 24 * t * t) - 28 * Math.pow(psi, 3) * (1 - 6 * t * t) + psi * psi * (1 - 32 * t * t) - psi * 2 * t * t + Math.pow(t, 4));
    const nterm4 = (Math.pow(w, 8) / 40320) * nu * sinlat * Math.pow(coslat, 7) * (1385 - 3111 * t * t + 543 * Math.pow(t, 4) - Math.pow(t, 6));
    const dN = k * (m + nterm1 + nterm2 + nterm3 + nterm4);
    const north = roundoff(0 + dN, distsize);

    return {north, east};
}

// LKS-94 to WGS-84
function grid2geo(x, y) {
    const k = 0.9998;
    const a = 6378137;
    const f = 1 / 298.257223563;
    const b = a * (1 - f);
    const e2 = (a * a - b * b) / (a * a);
    const n = (a - b) / (a + b);
    const G = a * (1 - n) * (1 - n * n) * (1 + (9 / 4) * n * n + (255 / 64) * Math.pow(n, 4)) * (Math.PI / 180);
    const north = (y - 0);
    const east = (x - 500000);
    const m = north / k;
    const sigma = (m * Math.PI) / (180 * G);
    const footlat = sigma + ((3 * n / 2) - (27 * Math.pow(n, 3) / 32)) * Math.sin(2 * sigma) + ((21 * n * n / 16) - (55 * Math.pow(n, 4) / 32)) * Math.sin(4 * sigma) + (151 * Math.pow(n, 3) / 96) * Math.sin(6 * sigma) + (1097 * Math.pow(n, 4) / 512) * Math.sin(8 * sigma);
    const rho = a * (1 - e2) / Math.pow(1 - (e2 * Math.sin(footlat) * Math.sin(footlat)), (3 / 2));
    const nu = a / Math.sqrt(1 - (e2 * Math.sin(footlat) * Math.sin(footlat)));
    const psi = nu / rho;
    const t = Math.tan(footlat);
    x = east / (k * nu);
    const laterm1 = (t / (k * rho)) * (east * x / 2);
    const laterm2 = (t / (k * rho)) * (east * Math.pow(x, 3) / 24) * (-4 * psi * psi + 9 * psi * (1 - t * t) + 12 * t * t);
    const laterm3 = (t / (k * rho)) * (east * Math.pow(x, 5) / 720) * (8 * Math.pow(psi, 4) * (11 - 24 * t * t) - 12 * Math.pow(psi, 3) * (21 - 71 * t * t) + 15 * psi * psi * (15 - 98 * t * t + 15 * Math.pow(t, 4)) + 180 * psi * (5 * t * t - 3 * Math.pow(t, 4)) + 360 * Math.pow(t, 4));
    const laterm4 = (t / (k * rho)) * (east * Math.pow(x, 7) / 40320) * (1385 + 3633 * t * t + 4095 * Math.pow(t, 4) + 1575 * Math.pow(t, 6));
    const latrad = footlat - laterm1 + laterm2 - laterm3 + laterm4;

    const lat_deg = rad2deg(latrad);

    const seclat = 1 / Math.cos(footlat);
    const loterm1 = x * seclat;
    const loterm2 = (Math.pow(x, 3) / 6) * seclat * (psi + 2 * t * t);
    const loterm3 = (Math.pow(x, 5) / 120) * seclat * (-4 * Math.pow(psi, 3) * (1 - 6 * t * t) + psi * psi * (9 - 68 * t * t) + 72 * psi * t * t + 24 * Math.pow(t, 4));
    const loterm4 = (Math.pow(x, 7) / 5040) * seclat * (61 + 662 * t * t + 1320 * Math.pow(t, 4) + 720 * Math.pow(t, 6));
    const w = loterm1 - loterm2 + loterm3 - loterm4;
    const longrad = deg2rad(24) + w;

    const lon_deg = rad2deg(longrad);
    return { lat_deg, lon_deg }
}