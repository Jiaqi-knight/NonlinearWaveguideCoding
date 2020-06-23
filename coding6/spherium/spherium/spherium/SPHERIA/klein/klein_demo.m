function klein_demo

% Granularity down the pipe (default 100)
pipe_N = 100;

% Rotational granularity (default 50)
rot_N = 50;

% Radius of the top bend relative to the small pipe (default 3)
rtb = 3;

% Radius of the base relative to the small pipe (default 5)
rb2sp = 7;

% Height of the conical base relative to the small pipe (default 10)
hcb2sp = 8;

figure('name','klein_demo','color',[1 1 1])
surface_handle = make_klein_bottle_plot( pipe_N, rot_N, rtb, rb2sp, hcb2sp );
camlight;
lighting phong;

%End of code