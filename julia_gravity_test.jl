# Wow such main

module Gravity

include("src_julia/eci2rdl.jl")
include("src_julia/mee2eci.jl")
include("src_julia/gravity_mee.jl")

scstate = [7088.42, -4.928, -1.579, -0.9797, 0.6117, 2.386]';
newstate = zeros(6,1);
#  [p, f, g, h, k, L]
#  [1, 2, 3, 4, 5, 6]
p = 1;
f = 2;
g = 3;
h = 4;
k = 5;
L = 6;

# global req mu j2 j3 j4 ;

# wgs-84 constants
mu = 398600.5; #// in km3 / s2
req = 6378.137; #// km
j2 = 0.00108262998905;
j3 = -0.00000253215306;
j4 = -0.00000161098761;

# Propagate the spacecraft for like a lot of steps or something


t = 0;
dt = .01;
t_max = 1000;

while (t < t_max)
	# calculate gravity
	gmee = gravity_mee(scstate);
	
	# propagate state
	s2 = 1 + scstate[h]*scstate[h] + scstate[k]*scstate[k];
	w = 1 + scstate[f] * cos(scstate[L]) + scstate[g] * sin(scstate[L]);
	r = scstate[p] / w;

	newstate[p] = scstate[p] + sqrt(scstate[p]/mu)*2*scstate[p]/w * gmee[2] * dt;
	newstate[f] = scstate[f] + sqrt(scstate[p]/mu)*( sin(scstate[L])*gmee[1] + ((w + 1)*cos(scstate[L]) + scstate[f])*gmee[2]/w - (scstate[h]*sin(scstate[L]) - k*cos(scstate[L]))*scstate[g] * gmee[3]/w ) * dt;
	newstate[g] = scstate[g] + sqrt(scstate[p]/mu)*( -cos(scstate[L])*gmee[1] + ((w + 1)*sin(scstate[L]) + scstate[g])*gmee[2]/w - (scstate[h]*sin(scstate[L]) - k*cos(scstate[L]))*scstate[f] * gmee[3]/w ) * dt;
	newstate[h] = scstate[h] + sqrt(scstate[p]/mu)*s2/(2*w)*cos(scstate[L])*gmee[3] * dt;
	newstate[k] = scstate[k] + sqrt(scstate[p]/mu)*s2/(2*w)*cos(scstate[L])*gmee[3] * dt;
	newstate[L] = scstate[L] + (sqrt(scstate[p]*mu)*w*w/scstate[p]/scstate[p] + 1/w*sqrt(scstate[p]/mu)*(scstate[h]*sin(scstate[L]) - scstate[k]*cos(scstate[L]))*gmee[3]) * dt;

	scstate = newstate;
	t += dt;
end

end
