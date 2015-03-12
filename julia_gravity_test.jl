# Wow such main

module Gravity

include("src_julia/eci2rdl.jl")
include("src_julia/mee2eci.jl")
include("src_julia/gravity_mee.jl")

#  [p, f, g, h, k, L]
scstate = [7088.42, -4.928, -1.579, -0.9797, 0.6117, 2.386]';
newstate = zeros(6,1);



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
	p = scstate[1];
	f = scstate[2];
	g = scstate[3];
	h = scstate[4];
	k = scstate[5];
	L = scstate[6];
	
	# propagate state
	s2 = 1 + h * h + k * k;
	w = 1 + f * cos(L) + g * sin(L);
	r = p / w;

	newstate[1] = p + sqrt(p/mu)*2*p/w * gmee[2] * dt;
	newstate[2] = f + sqrt(p/mu)*( sin(L)*gmee[1] + ((w + 1)*cos(L) + f)*gmee[2]/w - (h*sin(L) - k*cos(L))*g * gmee[3]/w ) * dt;
	newstate[3] = g + sqrt(p/mu)*( -cos(L)*gmee[1] + ((w + 1)*sin(L) + g)*gmee[2]/w - (h*sin(L) - k*cos(L))*f * gmee[3]/w ) * dt;
	newstate[4] = h + sqrt(p/mu)*s2/(2*w)*cos(L)*gmee[3] * dt;
	newstate[5] = k + sqrt(p/mu)*s2/(2*w)*cos(L)*gmee[3] * dt;
	newstate[6] = L + (sqrt(p*mu)*w*w/p/p + 1/w*sqrt(p/mu)*(h*sin(L) - k*cos(L))*gmee[3]) * dt;

	scstate = newstate;
	t += dt;
end

end
