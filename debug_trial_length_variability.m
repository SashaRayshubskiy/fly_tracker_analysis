%% Debug variability in trial length, yikes!

f = figure;
f1 = figure;

colormap jet;
cmap = colormap;
temp = linspace(1,size(cmap,1),max(trial_type_cnt));
cs = cmap(floor(temp),:);
close(gcf())

%for i = 1:size(trial_data,1)        
for i = 1
    for j = 1:size(trial_data{i},1)
        
        % dx = trial_data{ i, j }{2}.dx;
        % dy = trial_data{ i, j }{2}.dy;
        % t  = trial_data{ i, j }{2}.t;
        dx = trial_data{ i }{j,3}.dx;
        dy = trial_data{ i }{j,3}.dy;
        t = trial_data{ i }{j,3}.t;
        
        t_0 = t-t(1);
        
        figure(f)
        hold off;
        plot(t_0,cumsum(dx), 'color', cs(j,:));
        hold on;
        plot([20 20], [min(cumsum(dx)) max(cumsum(dx))], 'k' )
        title('Lateral');

        figure(f1)
        hold off;
        plot(t_0,cumsum(dy), 'color', cs(j,:));
        hold on;
        plot([20 20], [min(cumsum(dy)) max(cumsum(dy))], 'k' )
        title('Forward');
        
        waitforbuttonpress;
    end
end