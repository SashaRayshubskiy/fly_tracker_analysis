%% 

clear all;
close all;
datapath = '/Users/sasha/Documents/Wilson lab/Data/trackball/fly11';

search_path = '*_raw_trial_*.mat';
dirs = dir( [ datapath '/' search_path ]);

trial_type_cnt = zeros(4,1);

data = { [], [], [], [] };

for i = 1:length( dirs )
   
    filename = [datapath '/' dirs(i).name ];
    
    fs = strsplit(filename, '_');
    
    trial_type = [ fs{4} '_' fs{5} ];
    
    disp(['Trial type: ' trial_type]);
    
    d = load(filename);
    
    trial_type_idx = -1;
    if( strcmp(trial_type, 'Both_Air') )
        trial_type_idx = 1;
    elseif( strcmp(trial_type, 'Both_Odor') )
        trial_type_idx = 2;
    elseif( strcmp(trial_type, 'Left_Odor') )
        trial_type_idx = 3;
    elseif( strcmp(trial_type, 'Right_Odor') )
        trial_type_idx = 4;
    end
    
    trial_type_cnt(trial_type_idx) = trial_type_cnt(trial_type_idx) + 1;
    data{trial_type_idx, trial_type_cnt(trial_type_idx)} = d;
end

%% 
trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);

for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d = data{trial_idx, j};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);        

        subplot(2,2,trial_idx);
        
        [ traj_x, traj_y ] = calculate_traj(dx,dy);
        hold on;
        plot(traj_x, traj_y, 'color', cs(j,:));
        
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find(t_z >= PRE_STIM & t_z<(PRE_STIM+STIM));
        
        if(length(stim_t) > 0 )
            plot(traj_x(stim_t), traj_y(stim_t), 'x', 'color', cs(j,:));        
        end
    end 
    
    title(['Trial type: ' trial_type_labels{trial_idx} '     Count: ' num2str(trial_type_cnt(trial_idx))]);
end

%% Generate polar plots

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);

clear pre_dirs pre_vels stim_dirs stim_vels;
pre_dirs = {[],[],[],[]};
pre_vels = {[],[],[],[]};

stim_dirs = {[],[],[],[]};
stim_vels = {[],[],[],[]};

for trial_idx = 1:size(trial_type_cnt,1)
       
    for j=1:trial_type_cnt(trial_idx)
        d =  trial_data{ trial_idx, j }{2};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);        
        
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find(t_z >= PRE_STIM & t_z<(PRE_STIM+STIM));
        
        if(length(stim_t) > 0 )
            
            % Calculate avg direction and velocity
            dir_pre_x = sum(dx(pre_stim_t));
            dir_pre_y = sum(dy(pre_stim_t));
            pre_angle_rad = atan( dir_pre_y, dir_pre_x );
            pre_dirs{trial_idx,j} = pre_angle_rad;
            
            pre_vel_x = mean(dx(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t)));
            pre_vel_y = mean(dy(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t)));
            pre_vel = sqrt(pre_vel_x.^2 + pre_vel_y.^2);
            pre_vels{trial_idx,j} = pre_vel;
            
            dir_stim_x = sum(dx(stim_t));
            dir_stim_y = sum(dy(stim_t));
            stim_angle_rad = atan( dir_stim_y, dir_stim_x );
            stim_dirs{trial_idx,j} = stim_angle_rad;
            
            stim_vel_x = mean(dx(stim_t(2:end)) ./ diff(t_z(stim_t)));
            stim_vel_y = mean(dy(stim_t(2:end)) ./ diff(t_z(stim_t)));
            stim_vel = sqrt(stim_vel_x.^2 + stim_vel_y.^2);
            stim_vels{trial_idx,j} = stim_vel;

            %
           
            
            %[x,y] = pol2cart([pre_angle_rad],[pre_vel]);
            %compass(x,y);        
            % hold on
            %[x,y] = pol2cart([stim_angle_rad],[stim_vel]);
            %compass(x,y,'r');
            
            if 0
            subplot(2,2,trial_idx);
            hold on;
            polar([0 stim_angle_rad], [0 stim_vel], 'r');
            polar([0 pre_angle_rad], [0 pre_vel]);
            xlim([-1000 6000]);
            ylim([-1000 1500]);
            end
          
        end
    end 
        
    if 1
    subplot(2,2,trial_idx);
    
    pd = cell2mat(pre_dirs(trial_idx,:));
    pv = cell2mat(pre_vels(trial_idx,:));
    sd = cell2mat(stim_dirs(trial_idx,:));
    sv = cell2mat(stim_vels(trial_idx,:));

    if 0
    [x,y] = pol2cart(pd,pv);
    compass(x,y);    
    hold on;
    
    [x,y] = pol2cart(sd,sv);
    compass(x,y,'r');
    end
    
    if 1
    % Convert directions to a unit vector and take the angle of the average
    % vector
    [x,y] = pol2cart(pd,ones(1,length(pd)));  
    %compass(x,y); 
    avg_dir = atan2(mean(y), mean(x));
    [x,y] = pol2cart(avg_dir,mean(pv));
    compass(x,y, '--');
    
    [x,y] = pol2cart(sd,ones(1,length(sd)));    
    %compass(x,y,'r');
    hold on;
    avg_dir = atan2(mean(y),mean(x));
    [x,y] = pol2cart(avg_dir,mean(sv));
    compass(x,y, 'r--');
    end
    end
    title(['Trial type: ' trial_type_labels{trial_idx} '      Count: ' num2str(trial_type_cnt(trial_idx))]);
end

%% Generate pre-stim vs. stim polar plot diffs

trial_type_labels = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };

figure;
colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1), max(trial_type_cnt));
cs = cmap(floor(temp),:);

clear pre_dirs pre_vels stim_dirs stim_vels;
pre_dirs = {[],[],[],[]};
pre_vels = {[],[],[],[]};

stim_dirs = {[],[],[],[]};
stim_vels = {[],[],[],[]};

delta_dir_rad = {[],[],[],[]};

VEL_THRESHOLD = 3000;

for trial_idx = 1:size(trial_type_cnt,1)
    
    inc_idx = 1;
    for j=1:trial_type_cnt(trial_idx)
        d =  trial_data{ trial_idx, j }{2};
        
        t = d.t;
        dx = double(d.dx);
        dy = double(d.dy);        
        
        t_z = t-t(1);
        
        pre_stim_t = find(t_z < PRE_STIM);
        stim_t = find(t_z >= PRE_STIM & t_z<(PRE_STIM+STIM));
        
        if(length(stim_t) > 0 )
            
            % Calculate avg direction and velocity
            dir_pre_x = sum(dx(pre_stim_t));
            dir_pre_y = sum(dy(pre_stim_t));
            pre_angle_rad = atan2( dir_pre_y, dir_pre_x );
            
            pre_vel_x = mean(dx(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t)));
            pre_vel_y = mean(dy(pre_stim_t(2:end)) ./ diff(t_z(pre_stim_t)));
            pre_vel = pre_vel_y; 
            % pre_vel = sqrt(pre_vel_x.^2 + pre_vel_y.^2);
            
            dir_stim_x = sum(dx(stim_t));
            dir_stim_y = sum(dy(stim_t));
            stim_angle_rad = atan2( dir_stim_y , dir_stim_x );
            
            stim_vel_x = mean(dx(stim_t(2:end)) ./ diff(t_z(stim_t)));
            stim_vel_y = mean(dy(stim_t(2:end)) ./ diff(t_z(stim_t)));
            % stim_vel = sqrt(stim_vel_x.^2 + stim_vel_y.^2);
            stim_vel = stim_vel_y; 
            
            delta_pre_vs_stim_rad = atan2(sin(pre_angle_rad-stim_angle_rad), cos(stim_angle_rad-pre_angle_rad));

            %if( (pre_vel_y > VEL_THRESHOLD) && (stim_vel_y > VEL_THRESHOLD))                
            if( stim_vel_y > VEL_THRESHOLD)
                pre_dirs{trial_idx,inc_idx} = pre_angle_rad;
                pre_vels{trial_idx,inc_idx} = pre_vel;
                stim_dirs{trial_idx,inc_idx} = stim_angle_rad;
                stim_vels{trial_idx,inc_idx} = stim_vel;
                delta_dir_rad{trial_idx,inc_idx} = delta_pre_vs_stim_rad;
                inc_idx = inc_idx + 1;
            end
        end
    end 
        
    if 1
    subplot(2,2,trial_idx);
    
    pd = cell2mat(pre_dirs(trial_idx,:));
    pv = cell2mat(pre_vels(trial_idx,:));
    sd = cell2mat(stim_dirs(trial_idx,:));
    sv = cell2mat(stim_vels(trial_idx,:));
    dd = cell2mat(delta_dir_rad(trial_idx,:));
    
    if 0
        t = 0 : .01 : 2 * pi;
        P = polar(t, 5000 * ones(size(t)));
        set(P, 'Visible', 'off')
        
        hold on;
        [x,y] = pol2cart(-1*dd,sv-pv);
        compass(x,y);    
        hold on;
    
    %[x,y] = pol2cart(sd,sv);
    %compass(x,y,'r');
    end
    
    if 1
    % Convert directions to a unit vector and take the angle of the average
    % vector
    [x,y] = pol2cart(dd,ones(1,length(dd)));  
    %compass(x,y); 
    avg_dir = atan2(mean(y), mean(x));
    [x,y] = pol2cart(avg_dir,mean(sv-pv));
    
     t = 0 : .01 : 2 * pi;
     P = polar(t, 4000 * ones(size(t)));
     set(P, 'Visible', 'off')
     hold on;
        
    compass(x,y, '--');    
    end
    
    title(['Trial type: ' trial_type_labels{trial_idx} '      Count: ' num2str(trial_type_cnt(trial_idx))]);
end
end