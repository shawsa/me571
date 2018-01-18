N = 10^6;
v = 1:N;

tic
w = v+v;
toc

tic
w = v.*v;
toc

tic
w = sin(v);
toc