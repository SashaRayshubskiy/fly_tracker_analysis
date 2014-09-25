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

figure;
hist(v_f, 100000);
xlim([0 500]);
    
%% Load all mat files from trials

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

STIM = 15.0;
PRE_STIM = 5.0;

figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);

max_vsize = -1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d =  trial_data{ trial_idx, j }{2};
                
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

clear all_data_v all_data_vf all_data_vl;
all_data_vf = zeros(4*max(trial_type_cnt) * max_vsize,1 );
all_data_vl = zeros(4*max(trial_type_cnt) * max_vsize,1 );
all_data_v = zeros(4*max(trial_type_cnt) * max_vsize,1 );

cur_idx = 1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d =  trial_data{ trial_idx, j }{2};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
        
        t_diff = diff(t);

        v_x = dx(2:end) ./ t_diff;
        v_y = dy(2:end) ./ t_diff;

        v = sqrt(v_x.^2 + v_y.^2);
        
        all_data_vf(cur_idx:(cur_idx+size(v,2))-1) = v_y;
        all_data_vl(cur_idx:(cur_idx+size(v,2))-1) = v_x;
        all_data_v(cur_idx:(cur_idx+size(v,2))-1) = v;
        %all_data(cur_idx:(cur_idx+size(v,2))-1) = t_diff;
        cur_idx = cur_idx + size(v,2)
    end
end

clear all_data_vf_trim all_data_vl_trim all_data_v_trim;
all_data_vf_trim = all_data_vf(1:cur_idx);
all_data_vl_trim = all_data_vl(1:cur_idx);
all_data_v_trim  = all_data_v(1:cur_idx);

f = figure;
subplot(3,1,1);
hist(all_data_vf_trim,60);
xlabel('Velocity (au/s)','FontSize', 14);
ylabel('Count','FontSize', 14);
title('Forward velocity','FontSize', 16);

subplot(3,1,2);
hist(all_data_vl_trim,60);
xlabel('Velocity (au/s)','FontSize', 14);
ylabel('Count','FontSize', 14);
title('Lateral velocity','FontSize', 16);

subplot(3,1,3);
hist(all_data_v_trim,60);
xlabel('Velocity (au/s)','FontSize', 14);
ylabel('Count','FontSize', 14);
title('Velocity','FontSize', 16);

saveas(f, [basepath 'vel_hist_all_runs.png']);
saveas(f, [basepath 'vel_hist_all_runs.fig']);
saveas(f, [basepath 'vel_hist_all_runs.eps']);

%% 
trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

STIM = 15.0;
PRE_STIM = 5.0;

f = figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(f);

max_vsize = -1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d =  trial_data{ trial_idx, j }{2};
                
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

clear all_data_v_pre all_data_vf_pre all_data_vl_pre;
all_data_vf_pre = zeros(4*max(trial_type_cnt) * max_vsize,1 );
all_data_vl_pre = zeros(4*max(trial_type_cnt) * max_vsize,1 );
all_data_v_pre = zeros(4*max(trial_type_cnt) * max_vsize,1 );

clear all_data_v_stim all_data_vf_stim all_data_vl_stim;
all_data_vf_stim = zeros(4*max(trial_type_cnt) * max_vsize,1 );
all_data_vl_stim = zeros(4*max(trial_type_cnt) * max_vsize,1 );
all_data_v_stim = zeros(4*max(trial_type_cnt) * max_vsize,1 );

f = figure;
for trial_idx = 1:size(trial_type_cnt,1)
       
    cur_pre_idx = 1;
    cur_stim_idx = 1;
    for j=1:trial_type_cnt(trial_idx)
        d =  trial_data{ trial_idx, j }{2};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
        
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));
        
        pre_vel_x = dx(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        pre_vel_y = dy(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        pre_vel = sqrt(pre_vel_x.^2 + pre_vel_y.^2);        
        
        stim_vel_x = dx(stim_t(2:end)) ./ diff(t_z(stim_t));
        stim_vel_y = dy(stim_t(2:end)) ./ diff(t_z(stim_t));
        stim_vel = sqrt(stim_vel_x.^2 + stim_vel_y.^2);        
   
        cur_pre_size = size(pre_vel,2);
        cur_stim_size = size(stim_vel,2);
        
        all_data_vf_pre(cur_pre_idx:(cur_pre_idx+cur_pre_size-1)) = pre_vel_y;
        all_data_vl_pre(cur_pre_idx:(cur_pre_idx+cur_pre_size-1)) = pre_vel_x;
        all_data_v_pre(cur_pre_idx:(cur_pre_idx+cur_pre_size-1)) = pre_vel;
        cur_pre_idx = cur_pre_idx + cur_pre_size;
        
        all_data_vf_stim(cur_stim_idx:(cur_stim_idx+cur_stim_size-1)) = stim_vel_y;
        all_data_vl_stim(cur_stim_idx:(cur_stim_idx+cur_stim_size-1)) = stim_vel_x;
        all_data_v_stim(cur_stim_idx:(cur_stim_idx+cur_stim_size-1)) = stim_vel;
        cur_stim_idx = cur_stim_idx + cur_stim_size;
    end

    all_data_vl_pre_trim = all_data_vl_pre(1:cur_pre_idx);
    all_data_vf_pre_trim = all_data_vf_pre(1:cur_pre_idx);
    all_data_v_pre_trim = all_data_v_pre(1:cur_pre_idx);
    
    all_data_vl_stim_trim = all_data_vl_stim(1:cur_stim_idx);
    all_data_vf_stim_trim = all_data_vf_stim(1:cur_stim_idx);
    all_data_v_stim_trim = all_data_v_stim(1:cur_stim_idx);
 
    BIN_COUNT = 80;
    VEL_TYPE = 'Lat';
    subplot(2,4,2*trial_idx-1);
    hist(all_data_vl_pre_trim,BIN_COUNT);
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel pre-stim'], 'FontSize', 14);
    xlabel('Velocity (au/s)', 'FontSize', 14);
    ylabel('Count', 'FontSize', 14);
    xlim([-10000 10000]);
    ylim([0 1.8*16000]);
    
    subplot(2,4,2*trial_idx);
    hist(all_data_vl_stim_trim,BIN_COUNT );
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim'], 'FontSize', 14);
    xlabel('Velocity (au/s)', 'FontSize', 14);
    ylabel('Count', 'FontSize', 14);
    xlim([-10000 10000]);
    ylim([0 1.8*16000]);
end

saveas(f, [basepath '_' VEL_TYPE '_vel_hist_by_stim_type.png']);
saveas(f, [basepath '_' VEL_TYPE '_vel_hist_by_stim_type.eps']);
saveas(f, [basepath '_' VEL_TYPE '_vel_hist_by_stim_type.fig']);
