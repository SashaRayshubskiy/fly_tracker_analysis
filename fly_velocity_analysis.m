%%

basepath = '/home/sasha/fly_trackball_data/fly14/';
cd(basepath);

data = load([basepath '2014_0923_140441_raw_cummulative_xy.mat']);

t = data.t_all(1:end-1);
dx = data.dx_all(1:end-1);
dy = data.dy_all(1:end-1);

%%
integrateN = 18;
[t_d, dx_d, dy_d] = integrate_t_dx_dy(t,dx,dy,integrateN);

t_diff = diff(t_d);
t_zero = t_d-t_d(1);
v_f = dy_d(2:end) ./ t_diff;
v_l = dx_d(2:end) ./ t_diff;

hist(v_f, 100000);
xlim([0 500]);
    