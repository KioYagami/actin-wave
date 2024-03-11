clear variables;
close all;
v = 0.33;       % Actin wave velocity % WT
%v = 0.20;      % Actin wave velocity % Shootin1b-KO
length = 41;    % Length of actin waves
time = 20000;   % Simulation time

pcell.Base = 350;  % Base length of triangle (0.1 μm/pixel)
pcell.deg = 60;    % Vertex angle of triangle
pcell.R = (90-(pcell.deg/2))*(pi/180);  % Base angle of triangle (radian)
pcell.Height = (pcell.Base/2)*tan(pcell.R);
pcell.slope = pcell.Height/(pcell.Base/2);

figure(1)
DrawCellRegion(pcell)

x = zeros(1, 1000);
y = zeros(1, 1000);
prex = zeros(1, 1000);
prey = zeros(1, 1000);
rad = zeros(1, 1000);
vec_x = zeros(1, 1000);
vec_y = zeros(1, 1000);
status = zeros(1,1000);
Lifetime = zeros(1,1000);
generation_time = zeros(1,1000);
h = zeros(1,1000);
x_dis = zeros(1, 1000);
y_dis = zeros(1, 1000);
i = 0;        % Initial number of actin waves
veerFlag = 0;
LifetimeFlag = 0;
d = 1;
f = 1;
map = zeros(int16(pcell.Height)+100, pcell.Base+100);

%%
current_time = datetime('now');
seed_value = current_time.Second;
% Set the seed of random generator
rng(seed_value);

for t = 1:time  % Time

  generate = rand();

  if(generate > 0.9)    % Probability of actin wave generation
    i = i+1;    % Number of actin waves

    if(LifetimeFlag == 0)
      p = rand();
      q = rand();
      m1 = min(p, q);
      m2 = 1.0 - max(p, q);
      m3 = max(p, q) - min(p, q);
      x(i) = m1*0 + m2*(pcell.Base/2) + m3*pcell.Base;  % x coordinate
      y(i) = m1*0 + m2*pcell.Height + m3*0;             % y coordinate
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
    Lifetime(i) = gamrnd(3.016848, 118.9898); % Actin wave lifetime % WT
    %Lifetime(i) = gamrnd(2.397247, 106.2886); % Actin wave lifetime % Shootin1b-KO
    generation_time(i) = t; % Generated time

  end

  if(i>0)
    for num = 1:i   % Repeat for the number of actin waves

      if Lifetime(num) > 0
        status(num) = 0;  % On edge = 1
        if(mod((t - generation_time(num)), 50) == 0)
          veer = (11.1138*randn()-0.614419)*(pi/180);     % Directional change of actin waves
        else
          veer = 0;
        end

        rad(num) = rad(num) + veer;
        vec_x(num) = v*cos(rad(num));   % x vector
        vec_y(num) = v*sin(rad(num));   % y vector

        x(num) = x(num) + vec_x(num);   % Movement in x direction
        y(num) = y(num) + vec_y(num);   % Movement in y direction

        % Processing when actin waves collide with the plasma membrane
        if(y(num) >= pcell.slope*x(num))
          b = y(num) - (-1/pcell.slope)*x(num);
          x(num) = b/(pcell.slope-(-1/pcell.slope));
          y(num) = pcell.slope*x(num);
          status(num) = 1;
        elseif(y(num) >= -pcell.slope*x(num)+(pcell.Height*2))
          b = y(num) - (1/pcell.slope)*x(num);
          x(num) = (b-pcell.Height*2)/(-pcell.slope-(1/pcell.slope));
          y(num) = -pcell.slope*x(num)+(pcell.Height*2);
          status(num) = 1;
        end

        if(y(num) <= 0)
          y(num) = y(num) - vec_y(num);
          status(num) = 1;
        end

        if(y(num) >= pcell.Height) 
          y(num) = y(num) - vec_y(num);
          status(num) = 1;
        end

        if(x(num) <= 0)
          x(num) = x(num) - vec_x(num);
          status(num) = 1;
        end
        
        if(x(num) >= pcell.Base)
          x(num) = x(num) - vec_x(num);
          status(num) = 1;
        end

        % Remaining life
        Lifetime(num) = Lifetime(num) - 1;
        if(Lifetime(num) <= 0)
          x_dis(f) = x(num);
          y_dis(f) = y(num);
          LifetimeFlag = 1;
          f = f + 1;
          status(num) = 2;
        end

        % F-Actin intensity map
        a = tan(rad(num));
        b = y(num) - (a*x(num));
        if(cos(rad(num)) >= 0)
            for p = x(num) - (length*cos(rad(num))) : 0.1 : x(num)
                q = a*p + b;
                x_map = int16(p) + 50;
                y_map = int16(q) + 50;
                map(y_map, x_map) = map(y_map, x_map) + 1;
            end
        end
        if(cos(rad(num)) < 0)
            for p = x(num)  : 0.1 : x(num) - (length*cos(rad(num)))
                q = a*p + b;
                x_map = int16(p) + 50;
                y_map = int16(q) + 50;
                map(y_map, x_map) = map(y_map, x_map) + 1;
            end
        end
      end
    end  % num = 1:i
  end
end  % t = 1:time

%%

figure(2)
imagesc(-50:pcell.Base+50, -50:pcell.Height+50, map);
colormap jet;
colorbar;
axis xy;
daspect([1 1 1]);

% F-actin intensity within a 2 μm-wide region from the cell periphery
edge_map = zeros(int16(pcell.Height)+100, pcell.Base+100);
imageSize = [round(pcell.Height)+100, round(pcell.Base)+100];

x_1 = [50, pcell.Base/2+50, pcell.Base+50];
y_1 = [50, pcell.Height+50, 50];

mask_1 = poly2mask(x_1, y_1, imageSize(1), imageSize(2));

x_2 = [(20*tan((pi/2)-(pcell.R/2)))+50, pcell.Base/2+50, pcell.Base-(20*tan((pi/2)-(pcell.R/2)))+50];
y_2 = [20+50, (pcell.Base-((20*tan((pi/2)-(pcell.R/2)))*2))/2*tan(pcell.R)+20+50, 20+50];

mask_2 = poly2mask(x_2, y_2, imageSize(1), imageSize(2));

for edge_x = 1:round(pcell.Base)+100
    for edge_y = 1:round(pcell.Height)+100
        if(mask_1(edge_y, edge_x) == 1 && mask_2(edge_y, edge_x) == 0)
            edge_map(edge_y, edge_x) = map(edge_y, edge_x);
        end
    end
end

figure(3)
imagesc(-50:pcell.Base+50, -50:pcell.Height+50, edge_map);
colormap jet;
colorbar;
axis xy;
daspect([1 1 1]);

d = 20/cos(pcell.R);
% F-actin intensity within a 2 μm-wide and 5 μm-long regions in both sides of the vertex
corner_map = zeros(int16(pcell.Height)+100, pcell.Base+100);
n_1 = (pcell.Base/2)-(50*cos(pcell.R));
m_1 = pcell.Height-(50*sin(pcell.R));
b_1 = m_1 - ((-1/pcell.slope)*n_1);
for x = (pcell.Base/2)-(50*cos(pcell.R))*2:(pcell.Base/2)
  for y = pcell.slope*x-d:pcell.slope*x
    if(y > ((-1/pcell.slope)*x + b_1))
      edge_x = int16(x)+50;
      edge_y = int16(y)+50;
    end
    if(edge_y > pcell.Height+100)
      edge_y = int16(pcell.Height)+100;
    end
    corner_map(edge_y, edge_x) = edge_map(edge_y, edge_x);
  end
end

n_2 = (pcell.Base/2)+(50*cos(pcell.R));
m_2 = pcell.Height-(50*sin(pcell.R));
b_2 = m_2 - ((1/pcell.slope)*n_2);
for x = (pcell.Base/2):(pcell.Base/2)+(50*cos(pcell.R))*2
  for y = -pcell.slope*x+(pcell.Height*2)-d:-pcell.slope*x+(pcell.Height*2)
    if(y > ((1/pcell.slope)*x + b_2))
      edge_x = int16(x)+50;
      edge_y = int16(y)+50;
    end
    if(edge_y > pcell.Height+100)
      edge_y = int16(pcell.Height)+100;
    end
    corner_map(edge_y, edge_x) = edge_map(edge_y, edge_x);
  end
end
figure(4)
imagesc(-50:pcell.Base+50, -50:pcell.Height+50, corner_map);
colormap jet;
colorbar;
axis xy;
daspect([1 1 1]);

% Area of ​​whole cell, periphery, and corner
body_area = (pcell.Base*pcell.Height)/2;

edge_base = (pcell.Base-((20*tan((pi/2)-(pcell.R/2)))*2));
edge_height = edge_base/2*tan(pcell.R);
edge_area = body_area - (edge_base * edge_height / 2);

corner_area = (20*50*2) - (d * (d/2)/tan(pcell.R)*2 / 2) - (20*20*tan((90-pcell.deg)*(pi/180)));

body_intensity = sum(map, "all")/body_area;
edge_intensity = sum(edge_map, "all")/edge_area;
corner_intensity = sum(corner_map, "all")/corner_area;

disp("degree: " + pcell.deg)
disp("edge/body: " + (edge_intensity/body_intensity))
disp("corner/edge: " + (corner_intensity/edge_intensity))


function DrawCellRegion(pcell)
X1 = 0:pcell.Base;    % Base of triangle
Y1 = 0*X1;
hold on
plot(X1, Y1);

X2 = 0:(pcell.Base/2);     % Left slope of triangle
Y2 = pcell.slope*X2;
hold on
plot(X2, Y2);

X3 = (pcell.Base/2):pcell.Base;      % Right slope of triangle
Y3 = -pcell.slope*X3+(pcell.Height*2);
hold on
plot(X3, Y3);

xlim([-50 pcell.Base+50])
ylim([-50 pcell.Height+50])

daspect([1 1 1])

end

