%% Calculate turning index

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

STIM = 15.0;
PRE_STIM = 5.0;

f = figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(f)

f = figure;
f1 = figure;

pre_stim_sizes = [];
stim_sizes = [];

TIME_BEFORE_STIM_TO_EVALUATE_FOR_TURNING = 5.0; % seconds
TIME_AFTER_STIM_TO_EVALUATE_FOR_TURNING = 15.0; % seconds

for trial_idx = 1:size(trial_type_cnt,1)
       
    mean_turning_idx = zeros(trial_type_cnt(trial_idx),1);
    mean_speedup_idx = zeros(trial_type_cnt(trial_idx),1);
    
    cnt = 1;
    chosen_trials{trial_idx} = [];
        
    trial_cnt = trial_type_cnt(trial_idx);
    figure(f1);
    subplot(2,2,trial_idx);
    for j=1:trial_cnt

        d = trial_data{ trial_idx }{j,3};

        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
                        
        t_z = t-t(1);
        
        pre_stim_t = find((t_z > (PRE_STIM-TIME_BEFORE_STIM_TO_EVALUATE_FOR_TURNING)) & (t_z < PRE_STIM));
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+TIME_AFTER_STIM_TO_EVALUATE_FOR_TURNING)));
        
        % disp(['pre_stim size: ' num2str(size(pre_stim_t,2)) ' stim size: ' num2str(size(stim_t,2))]);
        
        if( size(pre_stim_t,2) <= 1 || (size(stim_t,2) <= 1 ))
            continue;
        end        
            
        pre_vel_x = dx(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        pre_vel_y = dy(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        % pre_vel = sqrt(pre_vel_x.^2 + pre_vel_y.^2);        
        
        stim_vel_x = dx(stim_t(2:end)) ./ diff(t_z(stim_t));
        stim_vel_y = dy(stim_t(2:end)) ./ diff(t_z(stim_t));
        % stim_vel = sqrt(stim_vel_x.^2 + stim_vel_y.^2);      
        
        mean_turning_idx(cnt) = mean(stim_vel_x) - mean(pre_vel_x);
        mean_speedup_idx(cnt) = mean(stim_vel_y) - mean(pre_vel_y);
        
        %cur_pre_stim_size(cnt) = size( pre_stim_t,2 );
        %cur_stim_size(cnt) = size( stim_t,2 );
        show_in_plot = 0;
        if((mean_turning_idx(cnt) < 0) & (trial_idx == 3))
            show_in_plot = 1;
        elseif(mean_turning_idx(cnt) > 0)
            show_in_plot = 1;
        end
                
        if (show_in_plot == 1 )
            pre_and_stim_t = horzcat( pre_stim_t, stim_t );
            [traj_x, traj_y] = calc_trial_trajectory( dx(pre_and_stim_t), dy(pre_and_stim_t) );
            hold on;
            plot(traj_x, traj_y, 'color', cs(j,:));
            plot( traj_x(size(pre_stim_t,2):end), traj_y(size(pre_stim_t,2):end), 'x', 'color', cs(j,:) );
            axis xy;
        end
              
        if(( trial_idx == 3 ) & (mean_turning_idx(cnt) < 0 ))
            chosen_trials{trial_idx} = [ chosen_trials{trial_idx} j];
        elseif(mean_turning_idx(cnt) > 0 )
            chosen_trials{trial_idx} = [ chosen_trials{trial_idx} j];
        end
        
        cnt = cnt + 1;
    end
    
    avg_turning_idx{trial_idx} = mean_turning_idx(1:cnt-1);
    avg_speedup_idx{trial_idx} = mean_speedup_idx(1:cnt-1);
    
    %pre_stim_sizes{ trial_idx } = cur_pre_stim_size(1:cnt-1);
    %stim_sizes{ trial_idx }     = cur_stim_size(1:cnt-1);
    
    
    if( trial_idx == 3 )
        turn_percent = size(find(avg_turning_idx{trial_idx} < 0),1) ./ size(avg_turning_idx{trial_idx},1);
    else
        turn_percent = size(find(avg_turning_idx{trial_idx} > 0),1) ./ size(avg_turning_idx{trial_idx},1);
    end
    
    
    speedup_percent = size(find(avg_speedup_idx{trial_idx} > 0),1) ./ size(avg_speedup_idx{trial_idx},1);
    
    figure(f)
    subplot(2,2,trial_idx);
    bar([avg_turning_idx{trial_idx} avg_speedup_idx{trial_idx}], 2.0, 'EdgeColor', 'none');  
        
    title( [trial_type_labels{trial_idx} ': turn\_idx: ' num2str(turn_percent) ' sd\_idx: ' num2str(speedup_percent) ], 'FontSize', 14);
    xlim([0 trial_cnt]);
    lh = legend('Turning index', 'Speed up index');
    %set(lh,'location','northeastoutside');
    xlabel('Trial #');
    ylabel('Index');

end

figname = 'turning_index';
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps']);

%% Process average velocity time course with error bars for lateral and forward 

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

STIM = 15.0;
PRE_STIM = 5.0;

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

for trial_idx = 1:size(trial_type_cnt,1)
       
    for i = 1:size(time_grid,2)
        time_grid_data{trial_idx,i} = [];
    end
    
    cur_chosen_trials = chosen_trials{trial_idx};
    for j=1:size(cur_chosen_trials,1);
        
        % WARNING: This only works if there are no skipped trials when
        % avg_turning_idx is generated
        jj = cur_chosen_trials(j); 
        %d =  trial_data{ trial_idx, j }{2};
        d = trial_data{ trial_idx }{jj,3};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
                
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));
        
        if( size(pre_stim_t,2) <= 1 || (size(stim_t,2) <= 1 ))
            continue;
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
            time_grid_data{trial_idx,time_grid_idx} = cat(1, time_grid_data{trial_idx,time_grid_idx}, [jj stim_vel_x(t_i-1) stim_vel_y(t_i-1)]);            
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
        
    cur_chosen_trials = chosen_trials{trial_idx};
    for j=1:size(cur_chosen_trials,1);
        
        % WARNING: This only works if there are no skipped trials when
        % avg_turning_idx is generated
        jj = cur_chosen_trials(j); 
        %d =  trial_data{ trial_idx, j }{2};
        d = trial_data{ trial_idx }{jj,3};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
                
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));
        
        if( size(pre_stim_t,2) <= 1 || (size(stim_t,2) <= 1 ))
            continue;
        end        
                 
        t_diff = diff(t_z);
        stim_vel_x = dx(2:end) ./ t_diff;
        stim_vel_y = dy(2:end) ./ t_diff;
        hold on;
        plot( t_z(2:end), stim_vel_x, 'color', rgb('Bisque') );
    end    
    
    hold on;
    VEL_TYPE = 'Lat';
    %plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    %plot(time_grid, sem_tc_lat, 'color', rgb('Bisque'));
    % subplot(1,2,2);
    %fh = fill([time_grid fliplr(time_grid)], ...
    %     [(avg_tc_lat+sem_tc_lat)' fliplr((avg_tc_lat-sem_tc_lat)')], ...
    %     rgb('Bisque'));
    %set(fh, 'EdgeColor', 'None');
     
    plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim error: std'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    xlim([0 25]);
    %ylim([-2000 2000]);

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
