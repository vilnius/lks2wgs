<?php

	// WGS-84 to LKS-94
	function geo2grid($lat, $lon) {
		$pi = pi();
		$distsize = 3;

		$j = 0; 		
		$units = 1;

		$latddd = $lat;
		$latrad = deg2rad($latddd);
		
		$londdd = $lon;
		$lonrad = deg2rad($londdd);
		
				$k = 0.9998;
				$a = 6378137;
				$f = 1 / 298.257223563;
				$b = $a * (1 - $f);
				$e2 = ($a * $a - $b * $b) / ($a * $a);
				$e = sqrt($e2);
				$ei2 = ($a * $a - $b * $b) / ($b * $b);
				$ei = sqrt($ei2);
				$n = ($a - $b) / ($a + $b);
				$G = $a * (1 - $n) * (1 - $n * $n) * (1 + (9 / 4) * $n * $n + (255 / 64) * pow($n, 4)) * ($pi / 180);
		
				$w = $londdd - 24;
				$w = deg2rad($w);
				$t = tan($latrad);
				$rho = $a * (1 - $e2) / pow(1 - ($e2 * sin($latrad) * sin($latrad)), (3 / 2));
				$nu = $a / sqrt(1 - ($e2 * sin($latrad) * sin($latrad)));
				
		$psi = $nu / $rho;
		$coslat = cos($latrad);
		$sinlat = sin($latrad);

				$A0 = 1 - ($e2 / 4) - (3 * $e2 * $e2 / 64) - (5 * pow($e2, 3) / 256);
		$A2 = (3 / 8) * ($e2 + ($e2 * $e2 / 4) + (15 * pow($e2, 3) / 128));
		$A4 = (15 / 256) * ($e2 * $e2 + (3 * pow($e2, 3) / 4));
		$A6 = 35 * pow($e2, 3) / 3072;
		$m = $a * (($A0 * $latrad) - ($A2 * sin(2 * $latrad)) + ($A4 * sin(4 * $latrad)) - ($A6 * sin(6 * $latrad)));

				$eterm1 = ($w * $w / 6) * $coslat * $coslat * ($psi - $t * $t);
		$eterm2 = (pow($w, 4) / 120) * pow($coslat, 4) * (4 * pow($psi, 3) * (1 - 6 * $t * $t) + $psi * $psi * (1 + 8 * $t * $t) - $psi * 2 * $t * $t + pow($t, 4));
		$eterm3 = (pow($w, 6) / 5040) * pow($coslat, 6) * (61 - 479 * $t * $t + 179 * pow($t, 4) - pow($t, 6));
		$dE = $k * $nu * $w * $coslat * (1 + $eterm1 + $eterm2 + $eterm3);
		$east = roundoff(500000 + ($dE / $units), $distsize);

				$nterm1 = ($w * $w / 2) * $nu * $sinlat * $coslat;
		$nterm2 = (pow($w, 4) / 24) * $nu * $sinlat * pow($coslat, 3) * (4 * $psi * $psi + $psi - $t * $t);
		$nterm3 = (pow($w, 6) / 720) * $nu * $sinlat * pow($coslat, 5) * (8 * pow($psi, 4) * (11 - 24 * $t * $t) - 28 * pow($psi, 3) * (1 - 6 * $t * $t) + $psi * $psi * (1 - 32 * $t * $t) - $psi * 2 * $t * $t + pow($t, 4));
		$nterm4 = (pow($w, 8) / 40320) * $nu * $sinlat * pow($coslat, 7) * (1385 - 3111 * $t * $t + 543 * pow($t, 4) - pow($t, 6));
		$dN = $k * ($m + $nterm1 + $nterm2 + $nterm3 + $nterm4);
		$north = roundoff(0 + ($dN / $units), $distsize);

				return array($north, $east);
	}

	// LKS-94 to WGS-84
	function grid2geo($lat, $lon) {
		$pi = pi();
		$distsize = 3;

		$j = 0; 		$units = 1;

				$k = 0.9998;  		$a = 6378137;  		$f = 1 / 298.257223563;  		$b = $a * (1 - $f); 		$e2 = ($a * $a - $b * $b) / ($a * $a);  		$e = sqrt($e2);  		$ei2 = ($a * $a - $b * $b) / ($b * $b); 		$ei = sqrt($ei2); 		$n = ($a - $b) / ($a + $b);
		$G = $a * (1 - $n) * (1 - $n * $n) * (1 + (9 / 4) * $n * $n + (255 / 64) * pow($n, 4)) * ($pi / 180); 
		$north = ($lat - 0) * $units;  		$east = ($lon - 500000) * $units; 		$m = $north / $k;  		$sigma = ($m * $pi) / (180 * $G);

				$footlat = $sigma + ((3 * $n / 2) - (27 * pow($n, 3) / 32)) * sin(2 * $sigma) + ((21 * $n * $n / 16) - (55 * pow($n, 4) / 32)) * sin(4 * $sigma) + (151 * pow($n, 3) / 96) * sin(6 * $sigma) + (1097 * pow($n, 4) / 512) * sin(8 * $sigma);
		$rho = $a * (1 - $e2) / pow(1 - ($e2 * sin($footlat) * sin($footlat)), (3 / 2)); 		$nu = $a / sqrt(1 - ($e2 * sin($footlat) * sin($footlat))); 		$psi = $nu / $rho;
		$t = tan($footlat);  				$x = $east / ($k * $nu);
		$laterm1 = ($t / ($k * $rho )) * ( $east * $x / 2);
		$laterm2 = ($t / ($k * $rho )) * ( $east * pow($x, 3) / 24) * (-4 * $psi * $psi + 9 * $psi * (1 - $t * $t) + 12 * $t * $t );
		$laterm3 = ($t / ($k * $rho )) * ( $east * pow($x, 5) / 720) * (8 * pow($psi, 4) * (11 - 24 * $t * $t) - 12 * pow($psi, 3) * (21 - 71 * $t * $t) + 15 * $psi * $psi * (15 - 98 * $t * $t + 15 * pow($t, 4)) + 180 * $psi * (5 * $t * $t - 3 * pow($t, 4)) + 360 * pow($t, 4));
		$laterm4 = ($t / ($k * $rho )) * ( $east * pow($x, 7) / 40320) * (1385 + 3633 * $t * $t + 4095 * pow($t, 4) + 1575 * pow($t, 6));
		$latrad = $footlat - $laterm1 + $laterm2 - $laterm3 + $laterm4;

		$lat_deg = rad2deg($latrad);

				$seclat = 1 / cos($footlat);
		$loterm1 = $x * $seclat;
		$loterm2 = (pow($x, 3) / 6) * $seclat * ($psi + 2 * $t * $t);
		$loterm3 = (pow($x, 5) / 120) * $seclat * (-4 * pow($psi, 3) * (1 - 6 * $t * $t) + $psi * $psi * (9 - 68 * $t * $t) + 72 * $psi * $t * $t + 24 * pow($t, 4));
		$loterm4 = (pow($x, 7) / 5040) * $seclat * (61 + 662 * $t * $t + 1320 * pow($t, 4) + 720 * pow($t, 6));
		$w = $loterm1 - $loterm2 + $loterm3 - $loterm4;
		$longrad = deg2rad(24) + $w;

		$lon_deg = rad2deg($longrad);


		return array($lat_deg, $lon_deg);
	}

	function roundoff($x, $y) {
		$x = round($x * pow(10, $y)) / pow(10, $y);
		return $x;
	}

?>
