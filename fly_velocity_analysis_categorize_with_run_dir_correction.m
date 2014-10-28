%% Plot per trial type runs

trial_types = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1),max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(gcf());

%figure('Name', ['Trial type: ' trial_types{i}]);    
f = figure('Name', ['Number of trials: ' num2str(size(files,1))]);    
f1 = figure();

for i = 1:size(trial_data,1)
        
    for j = 1:size(trial_data{i},1)
        
        % dx = trial_data{ i, j }{2}.dx;
        % dy = trial_data{ i, j }{2}.dy;
        % t  = trial_data{ i, j }{2}.t;
        dx = trial_data{ i }{j,3}.dx;
        dy = trial_data{ i }{j,3}.dy;
        t = trial_data{ i }{j,3}.t;
                               
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rotate the trial run by the direction of the pre_stim
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t_z = t-t(1);        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));
        dir_pre_x = sum(dx(pre_stim_t));
        dir_pre_y = sum(dy(pre_stim_t));
        pre_angle_rad = atan2( dir_pre_y, dir_pre_x );
        
        rot_rad = pre_angle_rad - pi/2.0;
        R = [cos(rot_rad) -sin(rot_rad); sin(rot_rad) cos(rot_rad)];

        v = double([dx; dy]');
        vR = v * R;
        dx_rot = vR(:,1);
        dy_rot = vR(:,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(f1)
        subplot(2,2,i);

        hold on;
        [traj_x_rot, traj_y_rot] = calc_trial_trajectory( dx_rot', dy_rot' );
        plot(traj_x_rot, traj_y_rot, 'color', cs(j,:));
        
        % label the start of stim with a 'X'       
        if(size(stim_t,1) > 0)
            plot( traj_x_rot(stim_t), traj_y_rot(stim_t), 'x', 'color', cs(j,:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(f);
        subplot(2,2,i);
        [traj_x, traj_y] = calc_trial_trajectory( dx, dy );
        plot(traj_x, traj_y, 'color', cs(j,:));
        hold on;
        
        % label the start of stim with a 'X'       
        if(size(stim_t,1) > 0)
            plot( traj_x(stim_t), traj_y(stim_t), 'x', 'color', cs(j,:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        xlabel('x distance (au)');
        ylabel('y distance (au)');
        %xlim([-10000 25000]); % fly 20
        %ylim([-1000 100000]); % fly 20
        %xlim([-5000 10000]); % fly 21
        %ylim([-1000 30000]); % fly 21
        %xlim([-6000 10000]); % fly 22
        %ylim([-1000 60000]); % fly 22
        %xlim([-3500 15000]); % fly 24
        %ylim([-1000 70000]); % fly 24
        %xlim([-10000 10000]); % fly 24
        %ylim([-1000 40000]); % fly 24

        %xlim([-6000 10000]); % fly 26
        %ylim([-1000 90000]); % fly 26

        xlim([-10000 10000]); % fly 26
        ylim([-1000 90000]); % fly 26

        title(['Trial type: ' trial_types{i} ' Trial #: ' num2str(trial_type_cnt(i))]);
        axis xy;
    end
end

search_dirs_for_output = strrep(search_dirs,'*','_wild_');

saveas(f, [basepath search_dirs_for_output '_stim_type_all_runs_pre_stim_corrected.png']);
saveas(f, [basepath search_dirs_for_output '_stim_type_all_runs_pre_stim_corrected.fig']);
saveas(f, [basepath search_dirs_for_output '_stim_type_all_runs_pre_stim_corrected.eps']);

%% Calculate turning index

trial_type_labels = { 'BA', 'BO', 'LO', 'RO', 'LA', 'RA' };

f = figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(f)

f = figure;

pre_stim_sizes = [];
stim_sizes = [];

PLOT_X_OFFSET_DELTA = 5000;

CHOOSE_PERCENT_CUTOFF = 0.0;
first_time = 1;

for trial_idx = 1:size(trial_type_cnt,1)
       
    mean_turning_idx = zeros(trial_type_cnt(trial_idx),1);
    mean_speedup_idx = zeros(trial_type_cnt(trial_idx),1);
    
    cnt = 1;
    correct_trials{trial_idx} = [];
    incorrect_trials{trial_idx} = [];
        
    trial_cnt = trial_type_cnt(trial_idx);
    
    subplot(2,3,trial_idx);
    
    cur_x_plot_offset = 0;
    
    qualified_trials_cnt = 0;
    for j=1:trial_cnt

        d = trial_data{ trial_idx }{j,3};

        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
                        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Rotate the trial run by the direction of the pre_stim
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t_z = t-t(1);        
                
        % pre_stim_t = find(t_z < PRE_STIM);
        LOOKBACK = PRE_STIM;
        % LOOKBACK = 5;
        pre_stim_t = find((t_z>(PRE_STIM-LOOKBACK)) & (t_z < PRE_STIM));
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));        
        
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
        
        if( size(pre_stim_t,2) <= 1 || (size(stim_t,2) <= 1 ))
            continue;
        end
        
        pre_vel_x = dx(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        pre_vel_y = dy(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t));
        % pre_vel = sqrt(pre_vel_x.^2 + pre_vel_y.^2);        
                
        stim_vel_x = dx(stim_t(2:end)) ./ diff(t_z(stim_t));
        stim_vel_y = dy(stim_t(2:end)) ./ diff(t_z(stim_t));
        % stim_vel = sqrt(stim_vel_x.^2 + stim_vel_y.^2);      
        
        QUALIFICATION_SPEED_LIMIT = 1800;
        if( (mean(pre_vel_y) < QUALIFICATION_SPEED_LIMIT) | (mean(stim_vel_y) < QUALIFICATION_SPEED_LIMIT ))
            continue;
        end

        mean_turning_idx(cnt) = (mean(stim_vel_x) - mean(pre_vel_x)) / abs(mean(pre_vel_x));
        
        mean_speedup_idx(cnt) = mean(stim_vel_y) - mean(pre_vel_y);

        choosen = 0;
        if( ((trial_idx == 3) | (trial_idx == 5)) & (mean_turning_idx(cnt) < -1.0*CHOOSE_PERCENT_CUTOFF) )
            choosen = 1;
        elseif ( ((trial_idx ~= 3) & (trial_idx ~= 5)) & (mean_turning_idx(cnt) > CHOOSE_PERCENT_CUTOFF) )
            choosen = 1;
        end
        
        if(choosen == 1)
            clr = rgb('Brown');
            correct_trials{trial_idx} = [ correct_trials{trial_idx} j];
            
            if(trial_idx == 3)
                cur_x_plot_offset = 0;
            else
                cur_x_plot_offset = 0;
            end
        else
            incorrect_trials{trial_idx} = [ incorrect_trials{trial_idx} j];
            clr = rgb('SandyBrown');
            
            if(trial_idx == 3)
                cur_x_plot_offset = 0;
            else
                cur_x_plot_offset = 0;
            end
        end
        
        pre_and_stim_t = horzcat( pre_stim_t, stim_t );
        [traj_x, traj_y] = calc_trial_trajectory( dx(pre_and_stim_t), dy(pre_and_stim_t) );
        hold on;
        ph = plot(traj_x  + cur_x_plot_offset, traj_y, 'color', clr);
        plot( traj_x(size(pre_stim_t,2):end) + cur_x_plot_offset, traj_y(size(pre_stim_t,2):end), 'x', 'color', clr );
        
        if(choosen)
            go_ph = ph;
        else
            no_go_ph = ph;
        end
        
        %text(cur_x_plot_offset, mean(traj_y(size(pre_stim_t,2):end)), ['TurnIdx: ' num2str(mean_turning_idx(cnt))]);
        %text(cur_x_plot_offset, mean(traj_y(size(pre_stim_t,2):end))+2000, ['PreVelX: ' num2str(mean(pre_vel_x))]);
        %text(cur_x_plot_offset, mean(traj_y(size(pre_stim_t,2):end))+4000, ['StimVelX: ' num2str(mean(stim_vel_x))]);
                
        cnt = cnt + 1;
        qualified_trials_cnt = qualified_trials_cnt + 1;
    end
    
    avg_turning_idx{trial_idx} = mean_turning_idx(1:cnt-1);
    avg_speedup_idx{trial_idx} = mean_speedup_idx(1:cnt-1);
    
    turn_percent = size(correct_trials{trial_idx}, 2) ./ qualified_trials_cnt;        
    speedup_percent = size(find(avg_speedup_idx{trial_idx} > 0),1) ./ size(avg_speedup_idx{trial_idx},1);
           
    title( [trial_type_labels{trial_idx} ': tid: ' num2str(turn_percent) ' (' num2str(size( correct_trials{trial_idx}, 2)) '/' ... 
        num2str(qualified_trials_cnt) ' : ' num2str(trial_cnt) ') sid: ' num2str(speedup_percent) ' (' ... 
        num2str(size(find(avg_speedup_idx{trial_idx} > 0),1)) '/' num2str(size(avg_speedup_idx{trial_idx},1)) ' : ' num2str(trial_cnt) ')'], ... 
        'FontSize', 14);
    
    %xlim([ -15000 15000 ]);
    %ylim([ 0 50000 ]);
    xlim([ -10000 10000 ]);
    ylim([ 0 40000 ]);
    lh = legend([go_ph; no_go_ph],{'Correct', 'Incorrect'});

    xlabel('X distance (au)', 'FontSize', 14);
    ylabel('Y distance (au)', 'FontSize', 14);    
end

figname = 'turning_index';
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps']);

%% Process average velocity time course with error bars for lateral and forward 

avg_vel_f = zeros(TRIAL_TYPE_CNT,max(trial_type_cnt));
avg_vel_l = zeros(TRIAL_TYPE_CNT,max(trial_type_cnt));

colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1),max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(gcf());

%f = figure;

clear time_grid;
TIME_GRID_SPACING = 50.0; % per second
TIME_GRID_SIZE    = PRE_STIM + STIM + FLUSH+5;
time_grid = [0 : 1.0/TIME_GRID_SPACING : TIME_GRID_SIZE ];

correct_time_grid_data = cell(TRIAL_TYPE_CNT,size(time_grid,2));

for trial_idx = 1:size(trial_type_cnt,1)
       
    for i = 1:size(time_grid,2)
        correct_time_grid_data{trial_idx,i} = [];
    end
    
    cur_correct_trials = correct_trials{trial_idx};
    for j=1:size(cur_correct_trials,2)
        
        % WARNING: This only works if there are no skipped trials when
        % avg_turning_idx is generated
        jj = cur_correct_trials(j); 
        %d =  trial_data{ trial_idx, j }{2};
        d = trial_data{ trial_idx }{jj,3};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);
                
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find((t_z >= PRE_STIM) & (t_z<(PRE_STIM+STIM)));
        
        % disp(['size(pre_stim_t): ' num2str(size(pre_stim_t))]);
        % disp(['size(stim_t): ' num2str(size(stim_t))]);
        
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
            correct_time_grid_data{trial_idx,time_grid_idx} = cat(1, correct_time_grid_data{trial_idx,time_grid_idx}, [jj stim_vel_x(t_i-1) stim_vel_y(t_i-1)]);            
        end
    end
end

% Process the incorrect trials
incorrect_time_grid_data = cell(TRIAL_TYPE_CNT,size(time_grid,2));
for trial_idx = 1:size(trial_type_cnt,1)
    for i = 1:size(time_grid,2)
        incorrect_time_grid_data{trial_idx,i} = [];
    end
    
    cur_incorrect_trials = incorrect_trials{trial_idx};
    for j=1:size(cur_incorrect_trials,2)
        
        % WARNING: This only works if there are no skipped trials when
        % avg_turning_idx is generated
        jj = cur_incorrect_trials(j); 
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
            
            incorrect_time_grid_data{trial_idx,time_grid_idx} = cat(1, incorrect_time_grid_data{trial_idx,time_grid_idx}, [jj stim_vel_x(t_i-1) stim_vel_y(t_i-1)]);
       end
   end
end

%% Plot average velocity time course with error bars for lateral and forward 

correct_avg_tc_lat = zeros(size(time_grid,2),1);
correct_sem_tc_lat = zeros(size(time_grid,2),1);
correct_avg_tc_fwd = zeros(size(time_grid,2),1);
correct_sem_tc_fwd = zeros(size(time_grid,2),1);

incorrect_avg_tc_lat = zeros(size(time_grid,2),1);
incorrect_sem_tc_lat = zeros(size(time_grid,2),1);
incorrect_avg_tc_fwd = zeros(size(time_grid,2),1);
incorrect_sem_tc_fwd = zeros(size(time_grid,2),1);

f = figure('units','normalized','outerposition',[0 0 1 1]);
f1 = figure('units','normalized','outerposition',[0 0 1 1]);
f2 = figure('units','normalized','outerposition',[0 0 1 1]);

for trial_idx = 1:size(trial_type_cnt,1)       
    
    for i = 1:size(time_grid,2)      

        disp(['size(stim_t): ' num2str(size(correct_time_grid_data{trial_idx,i}))]);

        if(size(correct_time_grid_data{trial_idx,i}) ~= 0 )
            mmm = mean(correct_time_grid_data{trial_idx,i},1);
            correct_avg_tc_lat( i ) = mmm( 2 );
            correct_avg_tc_fwd( i ) = mmm( 3 );
            
            sss = std(correct_time_grid_data{trial_idx,i},1);
            if(length(sss) == 1)
                correct_sem_tc_lat( i ) = 0.0;
                correct_sem_tc_fwd( i ) = 0.0;
            else            
                % correct_sem_tc_lat( i ) = sss( 2 ) / sqrt(size(correct_time_grid_data{trial_idx,i},1));
                % correct_sem_tc_fwd( i ) = sss( 3 ) / sqrt(size(correct_time_grid_data{trial_idx,i},1));
                correct_sem_tc_lat( i ) = sss( 2 );
                correct_sem_tc_fwd( i ) = sss( 3 );
            end
        end
    end
    
    for i = 1:size(time_grid,2)
        if(size(incorrect_time_grid_data{trial_idx,i}) ~= 0 )
            mmm = mean(incorrect_time_grid_data{trial_idx,i},1);
            incorrect_avg_tc_lat( i ) = mmm( 2 );
            incorrect_avg_tc_fwd( i ) = mmm( 3 );
            
            sss = std(incorrect_time_grid_data{trial_idx,i},1);
            if(length(sss) == 1)
                incorrect_sem_tc_lat( i ) = 0.0;
                incorrect_sem_tc_fwd( i ) = 0.0;
            else            
                % incorrect_sem_tc_lat( i ) = sss( 2 ) / sqrt(size(incorrect_time_grid_data{trial_idx,i},1));
                % incorrect_sem_tc_fwd( i ) = sss( 3 ) / sqrt(size(incorrect_time_grid_data{trial_idx,i},1));
                incorrect_sem_tc_lat( i ) = sss( 2 );
                incorrect_sem_tc_fwd( i ) = sss( 3 );
            end
        end
    end
    
    subplot_idx1 = -1;
    subplot_idx2 = -1;
    if (trial_idx == 1 )
        subplot_idx1 = 1;
        subplot_idx2 = 4;
    elseif (trial_idx == 2 )
        subplot_idx1 = 2;
        subplot_idx2 = 5;
    elseif (trial_idx == 3 )
        subplot_idx1 = 7;
        subplot_idx2 = 10;
    elseif (trial_idx == 4 )
        subplot_idx1 = 8;
        subplot_idx2 = 11;
    elseif (trial_idx == 5 )
        subplot_idx1 = 3;
        subplot_idx2 = 6;
    elseif (trial_idx == 6 )
        subplot_idx1 = 9;
        subplot_idx2 = 12;
    end
    
    figure(f);
    %subplot(4,2,subplot_idx);
    SPACING = 0.1;
    PADDING = 0;
    MARGIN = 0.05;
    subaxis(4,3,subplot_idx1, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);
    
    hold on;
    VEL_TYPE = 'Lat';
    %plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    %plot(time_grid, sem_tc_lat, 'color', rgb('Bisque'));
    % subplot(1,2,2);
    fh = fill([time_grid fliplr(time_grid)], ...
         [(correct_avg_tc_lat+correct_sem_tc_lat)' fliplr((correct_avg_tc_lat-correct_sem_tc_lat)')], ...
         rgb('Bisque'));
    set(fh, 'EdgeColor', 'None');
     
    plot(time_grid, correct_avg_tc_lat, 'color', rgb('Brown'));
    
    title(['Correct ' trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim error: std'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    xlim([0 PRE_STIM+STIM+FLUSH]);
    ylim([-3000 3000]);

    % subplot(4,2,subplot_idx+2);
    subaxis(4,3,subplot_idx2, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);

    hold on;
    VEL_TYPE = 'Lat';
    %plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    %plot(time_grid, sem_tc_lat, 'color', rgb('Bisque'));
    % subplot(1,2,2);
    fh = fill([time_grid fliplr(time_grid)], ...
         [(incorrect_avg_tc_lat + incorrect_sem_tc_lat)' fliplr((incorrect_avg_tc_lat-incorrect_sem_tc_lat)')], ...
         rgb('Bisque'));
    set(fh, 'EdgeColor', 'None');
     
    plot(time_grid, incorrect_avg_tc_lat, 'color', rgb('Brown'));
    
    title(['Incorrect ' trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim error: std'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    xlim([0 PRE_STIM+STIM+FLUSH]);
    ylim([-3000 3000]);    
    
    figure(f1);
    %subplot(4,2,subplot_idx);
    subaxis(4,3,subplot_idx1, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);

    hold on;
    VEL_TYPE = 'Fwd';
    %plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    %plot(time_grid, sem_tc_lat, 'color', rgb('Bisque'));
    % subplot(1,2,2);
    fh = fill([time_grid fliplr(time_grid)], ...
         [(correct_avg_tc_fwd+correct_sem_tc_fwd)' fliplr((correct_avg_tc_fwd-correct_sem_tc_fwd)')], ...
         rgb('Bisque'));
    set(fh, 'EdgeColor', 'None');
     
    plot(time_grid, correct_avg_tc_fwd, 'color', rgb('Brown'));
    
    title(['Correct ' trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim error: std'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    xlim([0 PRE_STIM+STIM+FLUSH]);
    ylim([-3000 6000]);

    %subplot(4,2,subplot_idx+2);
    subaxis(4,3,subplot_idx2, 'Spacing', SPACING, 'Padding', PADDING, 'Margin', MARGIN);

    hold on;
    VEL_TYPE = 'Fwd';
    %plot(time_grid, avg_tc_lat, 'color', rgb('Brown'));
    %plot(time_grid, sem_tc_lat, 'color', rgb('Bisque'));
    % subplot(1,2,2);
    fh = fill([time_grid fliplr(time_grid)], ...
         [(incorrect_avg_tc_fwd+incorrect_sem_tc_fwd)' fliplr((incorrect_avg_tc_fwd-incorrect_sem_tc_fwd)')], ...
         rgb('Bisque'));
    set(fh, 'EdgeColor', 'None');
     
    plot(time_grid, incorrect_avg_tc_fwd, 'color', rgb('Brown'));
    
    title(['Incorrect ' trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim error: std'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    xlim([0 PRE_STIM+STIM+FLUSH]);
    ylim([-3000 6000]);
    
    figure(f2);
    subplot(2,3,trial_idx);
    
    [ax1, h1, h2] = plotyy(time_grid, correct_avg_tc_lat, time_grid, correct_avg_tc_fwd);
    title([trial_type_labels{trial_idx} ': ' VEL_TYPE ' vel stim error: std'], 'FontSize', 14);
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Velocity (au/s)', 'FontSize', 14);
    
    xlim(ax1(1), [0 PRE_STIM+STIM+FLUSH]);
    xlim(ax1(2), [0 PRE_STIM+STIM+FLUSH]);
    ylim(ax1(1), [-3000 6000]);
    ylim(ax1(2), [-3000 6000]);
    legend([h1,h2], {'Lat', 'Fwd'});    

end

figname = ['lat_vel_ts_per_type'];
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps'], 'epsc');

figname = ['fwd_vel_ts_per_type'];
saveas(f1, [basepath figname '.png']);
saveas(f1, [basepath figname '.fig']);
saveas(f1, [basepath figname '.eps'], 'epsc');

figname = ['fwd_lat_vel_ts_per_type'];
saveas(f2, [basepath figname '.png']);
saveas(f2, [basepath figname '.fig']);
saveas(f2, [basepath figname '.eps'], 'epsc');
