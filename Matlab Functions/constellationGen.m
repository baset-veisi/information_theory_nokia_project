function x_m = constellationGen(N)
n = sqrt(N)/2; % n*n points in each quadrant
amp = 1+2*(n-1);
x_m = -amp:2:amp;

