%% Present runs per stim type

% Load trial data
clear all;
close all;
% basepath = '/home/sasha/fly_trackball_data/fly16/';
basepath = '/Users/sasha/Documents/Wilson lab/Data/trackball/fly20/';
cd(basepath);

search_dirs = '*raw_trial_*';
files = dir([search_dirs '.mat']);

clear trial_data;
trial_data = cell(4,1);

trial_type_cnt = ones(4,1);

for i=1:size(files,1)
    filename = files(i).name;
    fs = strsplit(filename, '_');
    tmp = strsplit(fs{8}, '.');
    
    trial_ord = tmp{1};    
    
    trial_type = [fs{4} '_' fs{5}];
    trial_id = [ fs{1} '_' fs{2} '_' fs{3} '_' trial_ord ];
       
    trial_type_idx = -1;
    
    if( strcmp(trial_type, 'Both_Air') )
        trial_type_idx = 1;
    elseif ( strcmp(trial_type, 'Both_Odor') )
        trial_type_idx = 2;        
    elseif ( strcmp(trial_type, 'Left_Odor') )
        trial_type_idx = 3;        
    elseif ( strcmp(trial_type, 'Right_Odor') )
        trial_type_idx = 4;        
    end
    
   disp(['Trial filename: ' filename]); 
   disp(['Trial id: ' trial_id]);    
   disp(['Trial type idx: ' num2str(trial_type_idx)]);
   disp(['Trial type: ' trial_type]);
   
   cur_idx = trial_type_cnt(trial_type_idx);
   
   trial_data{trial_type_idx,cur_idx} = { trial_id, load([basepath filename]) };
  
   trial_type_cnt(trial_type_idx) = cur_idx + 1;
   
   % disp(['Trial type cnt: ' num2str(trial_type_cnt(trial_type_idx))])
end

% subtract 1 from the indecies to get the count, because we started with
% base 1
trial_type_cnt = trial_type_cnt - 1;
%% Plot per trial type runs

trial_types = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };
STIM_ONSET = 5;
STIM_PERIOD = 15;

colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1),size(trial_data,2));
cs = cmap(floor(temp),:);
close(gcf());

%figure('Name', ['Trial type: ' trial_types{i}]);    
f = figure('Name', ['Number of trials: ' num2str(size(files,1))]);    
    
for i = 1:size(trial_data,1)
        
    subplot(2,2,i);
    for j = 1:size(trial_data,2)
       
        if (size(trial_data{ i, j },2) ~= 0)
            
            dx = trial_data{ i, j }{2}.dx;
            dy = trial_data{ i, j }{2}.dy;
            t  = trial_data{ i, j }{2}.t;
            
            [traj_x, traj_y] = calc_trial_trajectory( dx, dy );
            plot(traj_x, traj_y, 'color', cs(j,:));
            hold on;            

            % label the start of stim with a 'X'
            t_plot = t - t(1);
            
            t_plot_stim = find( (t_plot >= STIM_ONSET) & (t_plot <= (STIM_ONSET+STIM_PERIOD)));
             
            if(size(t_plot_stim,1) > 0)
                plot( traj_x(t_plot_stim), traj_y(t_plot_stim), 'x', 'color', cs(j,:));
            end
            
            xlabel('x distance (au)');
            ylabel('y distance (au)');
            xlim([-10000 25000]);
            ylim([-1000 100000]);
            title(['Trial type: ' trial_types{i} ' Trial #: ' num2str(trial_type_cnt(i))]);
            axis xy;
        end
    end
end

search_dirs_for_output = strrep(search_dirs,'*','_wild_');

saveas(f, [basepath search_dirs_for_output '_stim_type_all_runs.png']);
saveas(f, [basepath search_dirs_for_output '_stim_type_all_runs.fig']);
saveas(f, [basepath search_dirs_for_output '_stim_type_all_runs.eps']);

%% Plot avg distance traveled 

trial_types = { 'Both Air', 'Both Odor', 'Left Odor', 'Right Odor' };
STIM_ONSET = 5;
STIM_PERIOD = 15;

colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1),size(trial_data,2));
cs = cmap(floor(temp),:);
close(gcf());

%figure('Name', ['Trial type: ' trial_types{i}]);     
all_dist = zeros(size(trial_data,1),size(trial_data,2));
all_dist_stim = zeros(size(trial_data,1),size(trial_data,2));
all_dist_pre_stim = zeros(size(trial_data,1),size(trial_data,2));

avg_dist = zeros(size(trial_data,1),1);
std_dist = zeros(size(trial_data,1),1);

avg_dist_stim = zeros(size(trial_data,1),1);
std_dist_stim = zeros(size(trial_data,1),1);

avg_dist_pre_stim = zeros(size(trial_data,1),1);
std_dist_pre_stim = zeros(size(trial_data,1),1);

for i = 1:size(trial_data,1)
        
    for j = 1:size(trial_data,2)
       
        if (size(trial_data{ i, j },2) ~= 0)
            
            dx = trial_data{ i, j }{2}.dx;
            dy = trial_data{ i, j }{2}.dy;
            t  = trial_data{ i, j }{2}.t;
            
            % [traj_x, traj_y] = calc_trial_trajectory( dx, dy );
            
            all_dist(i,j) = sum(sqrt(double(dx.^2 + dy.^2)));          
                     
            t_plot = t - t(1);            
            t_plot_stim = find( (t_plot >= STIM_ONSET) & (t_plot <= (STIM_ONSET+STIM_PERIOD)));                        
            all_dist_stim(i,j) = sum(sqrt(double(dx(t_plot_stim).^2 + dy(t_plot_stim).^2)));

            t_plot_pre_stim = find( t_plot < STIM_ONSET);                        
            all_dist_pre_stim(i,j) = sum(sqrt(double(dx(t_plot_pre_stim).^2 + dy(t_plot_pre_stim).^2)));
                        
        end
    end
    
    cur_dist = nonzeros(all_dist(i,:));
    cur_dist_stim = nonzeros(all_dist_stim(i,:));
    cur_dist_pre_stim = nonzeros(all_dist_pre_stim(i,:));
    
    avg_dist( i ) = mean(cur_dist);
    std_dist( i ) = std(cur_dist);     

    avg_dist_stim( i ) = mean(cur_dist_stim);
    std_dist_stim( i ) = std(cur_dist_stim);     

    avg_dist_pre_stim( i ) = mean(cur_dist_pre_stim);
    std_dist_pre_stim( i ) = std(cur_dist_pre_stim);     
end

f = figure('Name', ['Number of trials: ' num2str(size(files,1))]);

subplot(1,3,1);
hold on;
bar(1:size(trial_data,1),avg_dist);
errorbar(1:size(trial_data,1),avg_dist, std_dist,'.');
set(gca,'Xtick',1:4,'XTickLabel',{'Both Air'; 'Both Odor'; 'Left Odor'; 'Right Odor'})
title(['Full trial']);
ylabel('Distance (au)');

subplot(1,3,2);
hold on;
bar(1:size(trial_data,1),avg_dist_stim);
errorbar(1:size(trial_data,1),avg_dist_stim, std_dist_stim,'.');
set(gca,'Xtick',1:4,'XTickLabel',{'Both Air'; 'Both Odor'; 'Left Odor'; 'Right Odor'})
title(['Stim period']);
ylabel('Distance (au)');

subplot(1,3,3);
hold on;
bar(1:size(trial_data,1),avg_dist_pre_stim);
errorbar(1:size(trial_data,1),avg_dist_pre_stim, std_dist_pre_stim,'.');
set(gca,'Xtick',1:4,'XTickLabel',{'Both Air'; 'Both Odor'; 'Left Odor'; 'Right Odor'})
title(['Pre stim period']);
ylabel('Distance (au)');

figname = '3_avg_distance';
saveas(f, [basepath figname '.png']);
saveas(f, [basepath figname '.fig']);
saveas(f, [basepath figname '.eps']);

