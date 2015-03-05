function gravity_mee(oemee)

# computes oblate earth gravitational acceleration 
# vector in the modified equinoctial element system

# input

#  oemee(1) = semiparameter (kilometers)
#  oemee(2) = f equinoctial element
#  oemee(3) = g equinoctial element
#  oemee(4) = h equinoctial element
#  oemee(5) = k equinoctial element
#  oemee(6) = true longitude (radians)

# output

#  gmee = gravity acceleration vector in modified
#         equinoctial element system (kilometers/second)^2

# global

#  req = earth equatorial radius (kilometers)
#  mu  = earth gravitational constant (kilometers^3/seconds^2)
#  j2, j3, j4 = zonal gravitation terms (non-dimensional)

#########################################

	xhat = zeros(1, 3);

	# compute eci state vector

	(reci, veci) = mee2eci(mu, oemee);

	geci = [0.; 0.; 0.];
	
	gravj = [0.; 0.; 0.; 0.];

	gravj[2] = j2;

	gravj[3] = j3;

	gravj[4] = j4;

	radius = norm(reci);

	# construct unit vectors for local horizontal frame

	zhat = -reci / radius;

	zdotr = zhat[3];

	for i = 1:3
		
		xhat[i] = -zdotr * zhat[i];
		
	end

	xhat[3] = 1.0 + xhat[3];

	xnorm = norm(xhat);

	xhat = xhat / xnorm;

	# compute sin and cosine of latitude

	sinphi = reci[3] / radius;

	cosphi = sqrt(1.0 - sinphi^2);

	glhz = 0.0;

	glhx = 0.0;

	aovrr = req / radius;

	aovrn = aovrr;

	pnm2 = 1.0;

	pnm1 = sinphi;

	ppnm2 = 0.0;

	ppnm1 = 1.0;

	for i = 2:4
		
		# evaluate legendre polynomial for this pass
		
		twonm1 = 2.0 * i - 1.0;
		
		onemn = 1.0 - i;
		
		denomi = i;
		
		pn = (sinphi * twonm1 * pnm1 + onemn * pnm2) / denomi;
		
		ppn = twonm1 * pnm1 + ppnm2;
		
		# construct local horizontal x direction sum
		
		aovrn = aovrn * aovrr;
		
		glhx = glhx + aovrn * ppn * gravj[i];
		
		# construct local horizontal z direction sum
		
		xnp1 = i + 1.0;
		
		glhz = glhz + xnp1 * aovrn * pn * gravj[i];
		
		# update legendre polynomials for next pass
		
		pnm2 = pnm1;
		
		pnm1 = pn;
		
		ppnm2 = ppnm1;
		
		ppnm1 = ppn;
		
	end

	cgovrr = mu / radius^2;

	glhx = -cosphi * cgovrr * glhx;

	glhz = -cgovrr * glhz;

	# complete the calculations by forming the gravitational
	# accelerations in the eci frame

	for i = 1:3
		
		geci[i] = glhx * xhat[i] + glhz * zhat[i];
		
	end

	# compute radial frame unit vectors

	(xrdl, yrdl, zrdl) = eci2rdl (reci, veci);

	# transform eci gravity vector to mee gravity components

	gmee = [0.; 0.; 0.];
	gmee[1] = dot(geci, xrdl);

	gmee[2] = dot(geci, yrdl);

	gmee[3] = dot(geci, zrdl);

	return gmee;

end;
