clear variables;
close all;
v = 0.33;       % Actin wave velocity % WT
%v = 0.20;      % Actin wave velocity % Shootin1b-KO
length = 41;    % Length of actin waves
time = 5700;    % Simulation time

pa.Base = 238;   % Base length of triangle (0.1 um/pixel)
pa.deg = 120;    % Vertex angle of triangle
pa.R = (90-(pa.deg/2))*(pi/180);    % Base angle of triangle (radian)
pa.Height = (pa.Base/2)*tan(pa.R);
pa.slope = pa.Height/(pa.Base/2);
pa.Area = (pa.Base*pa.Height)/2;
pa.Vert = (30725-pa.Area)/pa.Base;

C = -pi:0.01:0;     % Half circle
cir_x = (pa.Base/2)*cos(C)+(pa.Base/2);
cir_y = (pa.Base/2)*sin(C)+(pa.Base/2);

tri_x = [pa.Base, pa.Base, pa.Base/2, 0, 0];
tri_y= [pa.Base/2, pa.Base/2 + pa.Vert, round(pa.Height + pa.Vert + (pa.Base/2)), pa.Base/2 + pa.Vert, pa.Base/2];

pa.x = [cir_x tri_x];
pa.y = [cir_y tri_y];
pa.xl = [min(pa.x)-50 max(pa.x)+50];
pa.yl = [min(pa.y)-50 530];%max(pa.y)+80];
[pa.xgrid, pa.ygrid] = meshgrid(pa.xl(1):pa.xl(2), pa.yl(1):pa.yl(2));

f = figure;
f.Position = [100 100 900 400];

x = zeros(1, 1000);
y = zeros(1, 1000);
rad = zeros(1, 1000);
vec_x = zeros(1, 1000);
vec_y = zeros(1, 1000);
status = zeros(1,1000);
Lifetime = zeros(1,1000);
generation_time = zeros(1,1000);
x_dis = zeros(1, 1000);
y_dis = zeros(1, 1000);
h = zeros(1,1000);
i = 0;        % Initial number of actin waves
veerFlag = 0;
LifetimeFlag = 0;
d = 1;
f = 1;


sigma = 15;
kern = @(x) exp(- 0.5 * (x/sigma).^2);
total_I = zeros(size(pa.xgrid));

%%
current_time = datetime('now');
seed_value = current_time.Second;
% Set the seed of random generator
rng(seed_value);

outputVideo = VideoWriter('output_video.avi'); % Name of movie file
outputVideo.FrameRate = 60; % Framerate
open(outputVideo);

for t = 1:time  % Time

    generate = rand();

    if(generate > 0.9)  % Probability of actin wave generation
        i = i+1;    % Number of actin waves
        if(LifetimeFlag == 0)
            x(i) = rand()*pa.Base;      % x coordinate
            if(x(i) <= (pa.Base/2))     % Left side 
                C1 = -acos((x(i)-(pa.Base/2))/(pa.Base/2));
                y(i) = rand()*(pa.slope*x(i)+pa.Vert+(pa.Base/2)-((pa.Base/2)*sin(C1)+(pa.Base/2)))+((pa.Base/2)*sin(C1)+(pa.Base/2));     % y coordinate
            end
            if(x(i) > (pa.Base/2))      % Right side
                C1 = -acos((x(i)-(pa.Base/2))/(pa.Base/2));
                y(i) = rand()*(-pa.slope*x(i)+pa.Vert+(pa.Base/2)+(pa.Height*2)-((pa.Base/2)*sin(C1)+(pa.Base/2)))+((pa.Base/2)*sin(C1)+(pa.Base/2));    % y coordinate
            end
        elseif(LifetimeFlag == 1)
            x(i) = x_dis(d);
            y(i) = y_dis(d);
            d = d + 1;
            if(d >= f)
                LifetimeFlag = 0;
            end
        end
        rad(i) = 2*pi*rand();       % Direction of actin wave movement
        vec_x(i) = v*cos(rad(i));   % x vector
        vec_y(i) = v*sin(rad(i));   % y vector
        Lifetime(i) = gamrnd(3.016848, 118.9898);   % Actin wave lifetime % WT
        %Lifetime(i) = gamrnd(2.397247, 106.2886);  % Actin wave lifetime % Shootin1b-KO
        generation_time(i) = t; % Generated time
    end
    if(i>0)
        for num = 1:i  % Repeat for the number of actin waves
            if Lifetime(num) > 0
                status(num) = 0;
                if(mod((t - generation_time(num)), 50) == 0)
                    veer = (11.1138*randn()-0.614419)*(pi/180);     % Directional change of actin waves
                else
                    veer = 0;
                end

                rad(num) = rad(num) - veer;
                vec_x(num) = v*cos(rad(num));   % x vector
                vec_y(num) = v*sin(rad(num));   % y vector

                x(num) = x(num) + vec_x(num);   % Movement in x direction
                y(num) = y(num) + vec_y(num);   % Movement in y direction

                if(x(num)<=0)       % When actin waves pass through left base
                    x(num) = x(num) - vec_x(num);
                    status(num) = 1;
                end
                if(x(num)>=pa.Base) % When actin waves pass through right base
                    x(num) = x(num) - vec_x(num);
                    status(num) = 1;
                end

                % Processing of actin waves at the slope
                if(y(num)>=pa.slope*x(num)+pa.Vert+(pa.Base/2) && x(num)<=(pa.Base/2) && pa.deg~=180)   % Left slope
                    b = y(num) - (-1/pa.slope)*x(num); % A straight line perpendicular to the left slope and on x(num), y(num)
                    x(num) = (b-pa.Vert-(pa.Base/2))/(pa.slope-(-1/pa.slope));  % x coordinate of intersection of straight line and left slope
                    y(num) = pa.slope*x(num)+pa.Vert+(pa.Base/2);          % y coordinate of intersection of straight line and left slope
                    status(num) = 1;
                elseif(y(num) >= -pa.slope*x(num)+(pa.Height*2)+pa.Vert+(pa.Base/2) && x(num)>(pa.Base/2) && pa.deg~=180)  % Right slope
                    b = y(num) - (1/pa.slope)*x(num); % A straight line perpendicular to the right slope and on x(num), y(num)
                    x(num) = (b-pa.Height*2-pa.Vert-(pa.Base/2))/(-pa.slope-(1/pa.slope));  % x coordinate of intersection of straight line and right slope
                    y(num) = -pa.slope*x(num)+(pa.Height*2)+pa.Vert+(pa.Base/2);            % y coordinate of intersection of straight line and right slope
                    status(num) = 1;
                end

                if((y(num)<=(pa.Base/2)) && (x(num)<=(pa.Base/2)) && ((x(num)-(pa.Base/2))^2+(y(num)-(pa.Base/2))^2 >= (pa.Base/2)^2))   % Left half ciercle
                    theta = atan(((pa.Base/2)-y(num))/((pa.Base/2)-x(num)));
                    x(num) = (pa.Base/2)*cos(-pi+theta)+(pa.Base/2);  % x coordinate on circle
                    y(num) = (pa.Base/2)*sin(-pi+theta)+(pa.Base/2);  % y coordinate on circle
                    status(num) = 1;
                elseif((y(num)<=(pa.Base/2)) && (x(num)>=(pa.Base/2)) && ((x(num)-(pa.Base/2))^2+(y(num)-(pa.Base/2))^2 >= (pa.Base/2)^2))   % Right half ciercle
                    theta = atan(((pa.Base/2)-y(num))/(x(num)-(pa.Base/2)));
                    x(num) = (pa.Base/2)*cos(-theta)+(pa.Base/2);   % x coordinate on circle
                    y(num) = (pa.Base/2)*sin(-theta)+(pa.Base/2);   % y coordinate on circle
                    status(num) = 1;
                end

                if(y(num)>=pa.Height+pa.Vert+(pa.Base/2))
                    y(num) = y(num) - vec_y(num);
                    status(num) = 1;
                end

                % Remaining life
                Lifetime(num) = Lifetime(num) -1;
                if(Lifetime(num) <= 0)
                    x_dis(f) = x(num);
                    y_dis(f) = y(num);
                    LifetimeFlag = 1;
                    f = f + 1;
                    status(num) = 2;
                end

            end % Lifetime(num) > 0
        end % num = 1:i

        % Drawing actin wave
        if(mod(t,10)==0)

            idx = find(Lifetime > 0);

            tail_x = length*cos(rad)/2;
            tail_y = length*sin(rad)/2;

            nan_idx = find(Lifetime <= 0);
            x(nan_idx) = NaN;
            y(nan_idx) = NaN;
            tail_x(nan_idx) = NaN;
            tail_y(nan_idx) = NaN;

            % Actin wave concentration
            I = zeros(size(pa.xgrid));
            for num_ = 1:size(idx,2)
                num = idx(num_);
                dx = pa.xgrid - x(num);
                dy = pa.ygrid - y(num);
                I = I + kern(sqrt(dx.^2 + dy.^2));
            end

            subplot(1,2,1);
            pcolor(pa.xgrid, pa.ygrid, I);
            shading interp;
            colormap('jet');
            total_I = total_I + I;
            hold on;

            patch(pa.x, pa.y, 'w', 'EdgeColor', 'w', 'FaceColor', 'none');

            % Drawing actin wave
            x_ = x(idx);
            y_ = y(idx);

            h = quiver(x_, y_, tail_x(idx), tail_y(idx), 'AutoScale', 'off', 'LineWidth', 1.5, 'Color', 'magenta');
            [h.ShowArrowHead] = deal('off');

            h = quiver(x_, y_, -tail_x(idx), -tail_y(idx), 'AutoScale', 'off', 'LineWidth', 1.5, 'Color', 'magenta');
            [h.ShowArrowHead] = deal('off');

            % Text
            %text(-20, 510, sprintf("Timer: %d", t), 'FontSize', 15, 'Color', 'white');
            %text(max(pa.x)-50, 500, '50 um', 'FontSize', 15, 'Color', 'white');
            %plot([max(pa.x)-30 max(pa.x)+20], [520 520], 'w-', 'LineWidth', 2);

            xlim(pa.xl);
            ylim(pa.yl);
            axis off; xticks([]); yticks([]);
            daspect([1 1 1]);

            subplot(1,2,2);
            pcolor(pa.xgrid, pa.ygrid, total_I);
            shading interp;
            colormap('jet');
            xlim(pa.xl);
            ylim(pa.yl);
            axis off; xticks([]); yticks([]);
            title('Average intensity');
            daspect([1 1 1]);

            drawnow;

            frame = getframe(gcf);
            writeVideo(outputVideo, frame);

            clf;
        end % Drawing actin wave

    end % i > 0
end % t = 1:time

close(outputVideo);
