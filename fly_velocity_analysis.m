%% Load all mat files from trials

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

f = figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(f)

max_vsize = -1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)

        d = trial_data{ trial_idx }{j,3};

        % d =  trial_data{ trial_idx, j }{2};
                
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
        
        % if( length(t) <= 1 )
        %     continue;
        % end
        
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
        %d =  trial_data{ trial_idx, j }{2};
        d = trial_data{ trial_idx }{j,3};
        
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
        cur_idx = cur_idx + size(v,2);
    end
end

clear all_data_vf_trim all_data_vl_trim all_data_v_trim;
all_data_vf_trim = all_data_vf(1:cur_idx);
all_data_vl_trim = all_data_vl(1:cur_idx);
all_data_v_trim  = all_data_v(1:cur_idx);

NUMBINS = 50;
f = figure;
subplot(3,1,1);
hist(all_data_vf_trim,NUMBINS);
xlabel('Velocity (au/s)','FontSize', 14);
ylabel('Count','FontSize', 14);
title('Forward velocity','FontSize', 16);

subplot(3,1,2);
hist(all_data_vl_trim,NUMBINS);
xlabel('Velocity (au/s)','FontSize', 14);
ylabel('Count','FontSize', 14);
title('Lateral velocity','FontSize', 16);

subplot(3,1,3);
hist(all_data_v_trim,NUMBINS);
xlabel('Velocity (au/s)','FontSize', 14);
ylabel('Count','FontSize', 14);
title('Velocity','FontSize', 16);

saveas(f, [basepath 'vel_hist_all_runs.png']);
saveas(f, [basepath 'vel_hist_all_runs.fig']);
saveas(f, [basepath 'vel_hist_all_runs.eps']);

%% 
trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

f = figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(f);

max_vsize = -1;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        %d =  trial_data{ trial_idx, j }{2};
        d = trial_data{ trial_idx }{j,3};       
        
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

mean_vf_stim = zeros(4,1);
mean_vl_stim = zeros(4,1);
mean_v_stim = zeros(4,1);
std_vf_stim = zeros(4,1);
std_vl_stim = zeros(4,1);
std_v_stim = zeros(4,1);

mean_vf_pre = zeros(4,1);
mean_vl_pre = zeros(4,1);
mean_v_pre = zeros(4,1);
std_vf_pre = zeros(4,1);
std_vl_pre = zeros(4,1);
std_v_pre = zeros(4,1);

correct_with_prestim = 1;

f = figure;
for trial_idx = 1:size(trial_type_cnt,1)
       
    cur_pre_idx = 1;
    cur_stim_idx = 1;
    for j=1:trial_type_cnt(trial_idx)
        %d =  trial_data{ trial_idx, j }{2};
        d = trial_data{ trial_idx }{j,3};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);

        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));

        if( size(pre_stim_t,2) <= 1 || (size(stim_t,2) <= 1 ))
            continue;
        end
        
        
        if (correct_with_prestim == 1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Rotate the trial run by the direction of the pre_stim
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dir_pre_x = sum(dx(pre_stim_t));
            dir_pre_y = sum(dy(pre_stim_t));
            pre_angle_rad = atan2( dir_pre_y, dir_pre_x );
            
            rot_rad = pre_angle_rad - pi/2.0;
            R = [cos(rot_rad) -sin(rot_rad); sin(rot_rad) cos(rot_rad)];
            
            v = double([dx; dy]');
            vR = v * R;
            dx = vR(:,1)';
            dy = vR(:,2)';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
        disp(['pre_stim_t: ' num2str(size(pre_stim_t,2))]);        
        
        pre_vel_x = dx(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        pre_vel_y = dy(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        pre_vel = sqrt(pre_vel_x.^2 + pre_vel_y.^2);        
        
        stim_vel_x = dx(stim_t(2:end)) ./ diff(t_z(stim_t));
        stim_vel_y = dy(stim_t(2:end)) ./ diff(t_z(stim_t));
        stim_vel = sqrt(stim_vel_x.^2 + stim_vel_y.^2);        
   
        cur_pre_size = size(pre_vel,2)
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

    disp(['cur_pre_idx ' num2str(cur_pre_idx)])
    
    all_data_vl_pre_trim = all_data_vl_pre(1:cur_pre_idx);
    all_data_vf_pre_trim = all_data_vf_pre(1:cur_pre_idx);
    all_data_v_pre_trim = all_data_v_pre(1:cur_pre_idx);
    
    mean_vf_pre(trial_idx) = mean(all_data_vf_pre_trim);
    mean_vl_pre(trial_idx) = mean(all_data_vl_pre_trim);
    mean_v_pre(trial_idx) = mean(all_data_v_pre_trim);
    sem_vf_pre(trial_idx) = std(all_data_vf_pre_trim,1)./sqrt(size(all_data_vf_pre_trim,1));
    sem_vl_pre(trial_idx) = std(all_data_vl_pre_trim,1)./sqrt(size(all_data_vl_pre_trim,1));
    sem_v_pre(trial_idx) = std(all_data_v_pre_trim,1)./sqrt(size(all_data_v_pre_trim,1));
    
    all_data_vl_stim_trim = all_data_vl_stim(1:cur_stim_idx);
    all_data_vf_stim_trim = all_data_vf_stim(1:cur_stim_idx);
    all_data_v_stim_trim = all_data_v_stim(1:cur_stim_idx);

    mean_vf_stim(trial_idx) = mean(all_data_vf_stim_trim);
    mean_vl_stim(trial_idx) = mean(all_data_vl_stim_trim);
    mean_v_stim(trial_idx) = mean(all_data_v_stim_trim);
    sem_vf_stim(trial_idx) = std(all_data_vf_stim_trim,1)./sqrt(size(all_data_vf_stim_trim,1));;
    sem_vl_stim(trial_idx) = std(all_data_vl_stim_trim,1)./sqrt(size(all_data_vl_stim_trim,1));;
    sem_v_stim(trial_idx) = std(all_data_v_stim_trim,1)./sqrt(size(all_data_v_stim_trim,1));;
 
    BIN_COUNT = 50
    VEL_TYPE = 'Lat';
    subplot(2,4,2*trial_idx-1);
    hist(all_data_vl_pre_trim,BIN_COUNT);
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel pre-stim'], 'FontSize', 14);
    xlabel('Velocity (au/s)', 'FontSize', 14);
    ylabel('Count', 'FontSize', 14);
    xlim([-10000 10000]);
    ylim([0 8*16000]);
    %ylim([0 11000]);
    
    subplot(2,4,2*trial_idx);
    hist(all_data_vl_stim_trim,BIN_COUNT );
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim'], 'FontSize', 14);
    xlabel('Velocity (au/s)', 'FontSize', 14);
    ylabel('Count', 'FontSize', 14);
    xlim([-10000 10000]);
    ylim([0 8*16000]);
    %ylim([0 11000]);
end

saveas(f, [basepath '_' VEL_TYPE '_vel_hist_by_stim_type.png']);
saveas(f, [basepath '_' VEL_TYPE '_vel_hist_by_stim_type.eps']);
saveas(f, [basepath '_' VEL_TYPE '_vel_hist_by_stim_type.fig']);

f = figure;

subplot(3,1,1);
mean_vf_bars(1,:) = mean_vf_pre;
mean_vf_bars(2,:) = mean_vf_stim;
sem_vf_bars(1,:) = sem_vf_pre;
sem_vf_bars(2,:) = sem_vf_stim;
barwitherr(sem_vf_bars', mean_vf_bars');
set(gca,'Xtick',1:4,'XTickLabel', ...
{['Both Air (' num2str(trial_type_cnt(1)) ')']; ['Both Odor (' num2str(trial_type_cnt(2)) ')']; ['Left Odor (' num2str(trial_type_cnt(3)) ')']; ['Right Odor (' num2str(trial_type_cnt(4)) ')'];});

ylabel('Velocity (au/s)');
title(['Avg forward velocity']);
lh = legend('Pre-stim', 'Stim');
set(lh,'location','northeastoutside');
ylim([0 4500]);

subplot(3,1,2);
mean_vl_bars(1,:) = mean_vl_pre;
mean_vl_bars(2,:) = mean_vl_stim;
sem_vl_bars(1,:) = sem_vl_pre;
sem_vl_bars(2,:) = sem_vl_stim;
barwitherr(sem_vl_bars', mean_vl_bars');
set(gca,'Xtick',1:4,'XTickLabel', ...
{['Both Air (' num2str(trial_type_cnt(1)) ')']; ['Both Odor (' num2str(trial_type_cnt(2)) ')']; ['Left Odor (' num2str(trial_type_cnt(3)) ')']; ['Right Odor (' num2str(trial_type_cnt(4)) ')'];});
ylabel('Velocity (au/s)');
title(['Avg lateral velocity']);
lh = legend('Pre-stim', 'Stim');
set(lh,'location','northeastoutside');
yl = ylim;
ylim([yl(1) 750]);

subplot(3,1,3);
mean_v_bars(1,:) = mean_v_pre;
mean_v_bars(2,:) = mean_v_stim;
sem_v_bars(1,:) = sem_v_pre;
sem_v_bars(2,:) = sem_v_stim;
barwitherr(sem_v_bars', mean_v_bars');
set(gca,'Xtick',1:4,'XTickLabel', ...
{['Both Air (' num2str(trial_type_cnt(1)) ')']; ['Both Odor (' num2str(trial_type_cnt(2)) ')']; ['Left Odor (' num2str(trial_type_cnt(3)) ')']; ['Right Odor (' num2str(trial_type_cnt(4)) ')'];});
ylabel('Velocity (au/s)');
title(['Avg velocity']);
lh = legend('Pre-stim', 'Stim');
set(lh,'location','northeastoutside');
ylim([0 5000]);

% subplot(1,1,1);
% vf_bars(1,:) = mean_vf_pre;
% vf_bars(2,:) = mean_vf_stim;
% bar3(vf_bars);
% set(gca,'Xtick',1:4,'XTickLabel',{'Both Air'; 'Both Odor'; 'Left Odor'; 'Right Odor'});
% set(gca,'Ytick',1:2,'YTickLabel',{'Pre-stim'; 'Stim'});
% title(['Forward velocity']);
% zlabel('Velocity (au/s)');

figname = '133407_avg_vel';
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps']);

%% Plot average velocity to figure out a turning index

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

avg_vel_f = zeros(4,max(trial_type_cnt));
avg_vel_l = zeros(4,max(trial_type_cnt));

correct_with_prestim = 1;
f = figure;
for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        %d =  trial_data{ trial_idx, j }{2};
        d = trial_data{ trial_idx }{j,3};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
        
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));
        
        if( size(pre_stim_t,2) <= 1 || (size(stim_t,2) <= 1 ))
            continue;
        end        
                
        if (correct_with_prestim == 1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Rotate the trial run by the direction of the pre_stim
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dir_pre_x = sum(dx(pre_stim_t));
            dir_pre_y = sum(dy(pre_stim_t));
            pre_angle_rad = atan2( dir_pre_y, dir_pre_x );
            
            rot_rad = pre_angle_rad - pi/2.0;
            R = [cos(rot_rad) -sin(rot_rad); sin(rot_rad) cos(rot_rad)];
            
            v = double([dx; dy]');
            vR = v * R;
            dx = vR(:,1)';
            dy = vR(:,2)';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        pre_vel_x = dx(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        pre_vel_y = dy(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        pre_vel = sqrt(pre_vel_x.^2 + pre_vel_y.^2);        
        
        stim_vel_x = dx(stim_t(2:end)) ./ diff(t_z(stim_t));
        stim_vel_y = dy(stim_t(2:end)) ./ diff(t_z(stim_t));
        stim_vel = sqrt(stim_vel_x.^2 + stim_vel_y.^2);        
   
        avg_vel_l(trial_idx, j) = mean(stim_vel_x); 
        avg_vel_f(trial_idx, j) = mean(stim_vel_y); 
    end

    VEL_TYPE = 'Lat';
    subplot(2,2,trial_idx);

    %plot(avg_vel_l(trial_idx,:));        
    %plot(avg_vel_f(trial_idx,:),'r');    
    colormap jet;
    bar([avg_vel_l(trial_idx,:)' avg_vel_f(trial_idx,:)'], 2.0, 'EdgeColor', 'none');
    
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim'], 'FontSize', 14);
    xlabel('Trial #', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    xlim([0 trial_type_cnt(trial_idx)]);
    ylim([-2500 5500]);
    legend('Lateral', 'Forward');
end

figname = 'all_vel_bars';
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps']);

%% Process average velocity time course with error bars for lateral and forward 

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

avg_vel_f = zeros(4,max(trial_type_cnt));
avg_vel_l = zeros(4,max(trial_type_cnt));

colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1),max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(gcf());

%f = figure;

clear time_grid;
TIME_GRID_SPACING = 50.0; % per second
TIME_GRID_SIZE    = 30;
time_grid = [0 : 1.0/TIME_GRID_SPACING : TIME_GRID_SIZE ];

time_grid_data = cell(4,size(time_grid,1));
correct_with_prestim = 1;

for trial_idx = 1:size(trial_type_cnt,1)
       
    for i = 1:size(time_grid,2)
        time_grid_data{trial_idx,i} = [];
    end

    for j=1:trial_type_cnt(trial_idx)
        %d =  trial_data{ trial_idx, j }{2};
        d = trial_data{ trial_idx }{j,3};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
                
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));
        
        if( size(pre_stim_t,2) <= 1 || (size(stim_t,2) <= 1 ))
            continue;
        end
        
        if (correct_with_prestim == 1)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Rotate the trial run by the direction of the pre_stim
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dir_pre_x = sum(dx(pre_stim_t));
            dir_pre_y = sum(dy(pre_stim_t));
            pre_angle_rad = atan2( dir_pre_y, dir_pre_x );
            
            rot_rad = pre_angle_rad - pi/2.0;
            R = [cos(rot_rad) -sin(rot_rad); sin(rot_rad) cos(rot_rad)];
            
            v = double([dx; dy]');
            vR = v * R;
            dx = vR(:,1)';
            dy = vR(:,2)';
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
                 
        t_diff = diff(t_z);
        stim_vel_x = dx(2:end) ./ t_diff;
        stim_vel_y = dy(2:end) ./ t_diff;
        
        time_grid_idx = 1;
        for t_i = 2:size(t,2)  
            while( time_grid(time_grid_idx) < t_z(t_i) )
                time_grid_idx = time_grid_idx  + 1;
            end
            
            % t(t_i) is => time_grid(time_grid_idx) here.
            %t_i
            time_grid_data{trial_idx,time_grid_idx} = cat(1, time_grid_data{trial_idx,time_grid_idx}, [j stim_vel_x(t_i-1) stim_vel_y(t_i-1)]);            
        end
    end
    
    if 0
    subplot(2,2,trial_idx);
    hold on;
    plot(t_z, dx, 'color', cs(j,:));
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    end
end

if 0
figname = 'all_vel_ts';
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps']);
end

%% Plot average velocity time course with error bars for lateral and forward 

avg_tc_lat = zeros(size(time_grid,2),1);
sem_tc_lat = zeros(size(time_grid,2),1);

avg_tc_fwd = zeros(size(time_grid,2),1);
sem_tc_fwd = zeros(size(time_grid,2),1);

f = figure;
f1 = figure;
for trial_idx = 1:size(trial_type_cnt,1)       
    for i = 1:size(time_grid,2)
        
        if(size(time_grid_data{trial_idx,i}) ~= 0 )
            mmm = mean(time_grid_data{trial_idx,i},1);
            avg_tc_lat( i ) = mmm( 2 );
            avg_tc_fwd( i ) = mmm( 3 );
            
            sss = std(time_grid_data{trial_idx,i},1);
            if(length(sss) == 1)
                sem_tc_lat( i ) = 0.0;
                sem_tc_fwd( i ) = 0.0;
            else            
                %sem_tc_lat( i ) = sss( 2 ) / sqrt(size(time_grid_data{trial_idx,i},1));
                %sem_tc_fwd( i ) = sss( 3 ) / sqrt(size(time_grid_data{trial_idx,i},1));
                sem_tc_lat( i ) = sss( 2 );
                sem_tc_fwd( i ) = sss( 3 );
            end
        end
    end
    
    figure(f);
    subplot(2,2,trial_idx);
    hold on;
    VEL_TYPE = 'Lat';
    %plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    %plot(time_grid, sem_tc_lat, 'color', rgb('Bisque'));
    % subplot(1,2,2);
    fh = fill([time_grid fliplr(time_grid)], ...
         [(avg_tc_lat+sem_tc_lat)' fliplr((avg_tc_lat-sem_tc_lat)')], ...
         rgb('Bisque'));
    set(fh, 'EdgeColor', 'None');
     
    plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim error: std'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    xlim([0 25]);
    ylim([-2000 2000]);

    figure(f1);
    subplot(2,2,trial_idx);
    hold on;
    VEL_TYPE = 'Fwd';
    %plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    %plot(time_grid, sem_tc_lat, 'color', rgb('Bisque'));
    % subplot(1,2,2);
    fh = fill([time_grid fliplr(time_grid)], ...
         [(avg_tc_fwd+sem_tc_fwd)' fliplr((avg_tc_fwd-sem_tc_fwd)')], ...
         rgb('Bisque'));
    set(fh, 'EdgeColor', 'None');
     
    plot(time_grid, avg_tc_fwd, 'color', rgb('Brown'));
    
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim error: std'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    xlim([0 25]);
    ylim([-2000 5000]);
end

figname = ['lat_vel_ts_per_type'];
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps']);

figname = ['fwd_vel_ts_per_type'];
saveas(f1, [basepath figname '.png']);
saveas(f1, [basepath figname '.fig']);
saveas(f1, [basepath figname '.eps']);
