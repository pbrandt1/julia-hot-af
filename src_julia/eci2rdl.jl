function eci2rdl (reci, veci)

# compute radial frame unit vectors

# input

#  reci = eci position vector (kilometers)
#  veci = eci velocity vector (kilometers/second)

# output

#  xrdl, yrdl, zrdl = radial frame unit vectors

# Orbital Mechanics with MATLAB

###############################

	xrdl = reci / norm(reci);

	hvec = cross(reci, veci);

	zrdl = hvec / norm(hvec);

	yrdl = cross(zrdl, xrdl);

	return xrdl, yrdl, zrdl;
end;
