%%

cfiles = dir([basepath '*raw_cummulative_xy*.mat']);

f = figure;
prev_x_start = 0;
prev_y_start = 0;
first_time = '';
last_time = '';
for i = 1:length(cfiles)

    filename = cfiles(i).name;
    disp(filename);
    fs = strsplit(filename, '_');
    
    hour = fs{3}(1:2);
    min  = fs{3}(3:4);
    sec  = fs{3}(5:6);
    date_time = [ fs{1} '-' fs{2} ' ' hour ':' min ':' sec ];
    
    if(i==1)
        first_time = date_time;
    end
    
    if(i==length(cfiles))
        last_time = date_time;
    end
    
    d = load([basepath filename]);
    
    t_all = d.t_all;
    dx_all = d.dx_all;
    dy_all = d.dy_all;
    
    [traj_x, traj_y] = calc_trial_trajectory(dx_all, dy_all, prev_x_start, prev_y_start);
    prev_x_start = traj_x(end);
    prev_y_start = traj_y(end);
    hold on;
    plot(traj_x, traj_y);
end

xlabel('Distance (au)', 'FontSize', 14);
ylabel('Distance (au)', 'FontSize', 14);
tt = title(['Cumulative run ( ' first_time ' to ' last_time ' )'], 'FontSize', 16);
set(tt,'interpreter','none')

figname = 'cumulative_run';
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps']);

%% Plot cummulative

[ traj_x, traj_y ] = calculate_traj(dx_all,dy_all);

figure;
plot(traj_x,traj_y)

%% Calculate velocity distributions

t_diff = diff(t_all);

v_x = dx_all(2:end) ./ t_diff;
v_y = dy_all(2:end) ./ t_diff;

v = sqrt(v_x.^2 + v_y.^2);

% hist(v,10000)
t_z = t_all-t_all(1);
figure;
plot(t_z(2:end),v)
axis tight;
xlabel('Time (s)');

%% 
figure;
hold on
plotyy([1:size(t_diff,2)],t_diff,[1:size(t_diff,2)],dy_all(2:end))

%%
%hist(t_diff,100000)
hist(v,100000)

%% Load all mat files from trials

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);

max_vsize = -1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d = data{trial_idx, j};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
        
        t_diff = diff(t);

        v_x = dx(2:end) ./ t_diff;
        v_y = dy(2:end) ./ t_diff;

        v = sqrt(v_x.^2 + v_y.^2);
        if( size(v,2) > max_vsize )
            max_vsize = size(v,2);
        end
        %all_data_v(trial_idx, j, 1:size(v,2)) = v;
    end
end

disp(['Max size: ' num2str(max_vsize)]);

clear all_data_v;
all_data = zeros(4*max(trial_type_cnt) * max_vsize,1 );

cur_idx = 1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d = data{trial_idx, j};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
        
        t_diff = diff(t);

        v_x = dx(2:end) ./ t_diff;
        v_y = dy(2:end) ./ t_diff;

        v = sqrt(v_x.^2 + v_y.^2);
        
        %all_data(cur_idx:(cur_idx+size(v,2))-1) = v;
        all_data(cur_idx:(cur_idx+size(v,2))-1) = t_diff;
        cur_idx = cur_idx + size(v,2);
    end
end

clear all_data_trim;
all_data_trim = all_data(1:cur_idx);
hist(all_data_trim,100000)