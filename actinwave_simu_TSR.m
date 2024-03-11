clear variables;
close all;
v = 0.33;       % Actin wave velocity % WT
%v = 0.20;      % Actin wave velocity % Shootin1b-KO#1
length = 41;    % Length of actin waves
time = 20000;   % Simulation time

Base = 238;  % Base length of triangle (0.1 μm/pixel)
deg = 60;    % Vertex angle of triangle

R = (90-(deg/2))*(pi/180); % Base angle of triangle (radian)
Height = (Base/2)*tan(R); 
Area = (Base*Height)/2;
slope = Height/(Base/2);

figure(1)
C = -pi:0.01:0;
plot((Base/2)*cos(C)+(Base/2), (Base/2)*sin(C)+(Base/2))

Vert = (30725-Area)/Base;
y = (Base/2):Vert+(Base/2);
x = 0*y;
hold on
plot(x, y);

y = (Base/2):Vert+(Base/2);
x = 0*y+Base;
hold on
plot(x, y);

x = 0:(Base/2);     % Left slope of triangle
y = slope*x+Vert+(Base/2);
hold on
plot(x, y);

x = (Base/2):Base;      % Right slope of triangle
y = -slope*x+(Height*2)+Vert+(Base/2);
hold on
plot(x, y);

xlim([-50 Base+50])
ylim([-50 Height+Vert+Base/2+50])
daspect([1 1 1])

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
map = zeros(int16(Height+Vert+(Base/2))+100, Base+100);

%%
current_time = datetime('now');
seed_value = current_time.Second;
% Set the seed of random generator
rng(seed_value);

for t = 1:time  % Time

    generate = rand();

    if(generate > 0.9)  % Probability of actin wave generation
        i = i+1;    % Number of actin waves
        if(LifetimeFlag == 0)
            x(i) = rand()*Base;     % x coordinate
            if(x(i) <= (Base/2))    % Left side
                C1 = -acos((x(i)-(Base/2))/(Base/2));
                y(i) = rand()*(slope*x(i)+Vert+(Base/2)-((Base/2)*sin(C1)+(Base/2)))+((Base/2)*sin(C1)+(Base/2));     % y coordinate
            end
            if(x(i) > (Base/2))     % Right side
                C1 = -acos((x(i)-(Base/2))/(Base/2));
                y(i) = rand()*(-slope*x(i)+Vert+(Base/2)+(Height*2)-((Base/2)*sin(C1)+(Base/2)))+((Base/2)*sin(C1)+(Base/2));    % y coordinate
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
        %Lifetime(i) = gamrnd(2.397247, 106.2886);  % Actin wave lifetime % Shootin1b-KO#1
        generation_time(i) = t; % Time of actin wave emergence
    end
    if(i>0)
        for num = 1:i   % Repeat for the number of actin waves
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
            
                % Processing when actin waves collide with the plasma membrane
                if(x(num)<=0)
                    x(num) = x(num) - vec_x(num);
                    status(num) = 1;
                end
                if(x(num)>=Base)
                    x(num) = x(num) - vec_x(num);
                    status(num) = 1;
                end
            
                if(y(num)>=slope*x(num)+Vert+(Base/2) && x(num)<=(Base/2) && deg~=180)
                    b = y(num) - (-1/slope)*x(num);
                    x(num) = (b-Vert-(Base/2))/(slope-(-1/slope));
                    y(num) = slope*x(num)+Vert+(Base/2);
                    status(num) = 1;
                elseif(y(num) >= -slope*x(num)+(Height*2)+Vert+(Base/2) && x(num)>(Base/2) && deg~=180)
                    b = y(num) - (1/slope)*x(num);
                    x(num) = (b-Height*2-Vert-(Base/2))/(-slope-(1/slope));
                    y(num) = -slope*x(num)+(Height*2)+Vert+(Base/2);
                    status(num) = 1;
                end
            
                if((y(num)<=(Base/2)) && (x(num)<=(Base/2)) && ((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2))
                    theta = atan(((Base/2)-y(num))/((Base/2)-x(num))); 
                    x(num) = (Base/2)*cos(-pi+theta)+(Base/2);
                    y(num) = (Base/2)*sin(-pi+theta)+(Base/2);
                    status(num) = 1;
                elseif((y(num)<=(Base/2)) && (x(num)>=(Base/2)) && ((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2))
                    theta = atan(((Base/2)-y(num))/(x(num)-(Base/2)));
                    x(num) = (Base/2)*cos(-theta)+(Base/2);
                    y(num) = (Base/2)*sin(-theta)+(Base/2);
                    status(num) = 1;
                end
            
                if(y(num)>=Height+Vert+(Base/2))
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
                
                % F-actin intensity map
                a = tan(rad(num));
                b = y(num) - (a*x(num));
                y_cir = 0;
               
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
        end
    end
end

disp("degree: " + deg)

figure(2)
imagesc(-50:Base+50, -50:Height+Vert+(Base/2)+50, map);
colormap jet;
colorbar;
axis xy;
daspect([1 1 1]);

% F-actin intensity within a 2 μm-wide region from the cell periphery
edge_map = zeros(int16(Height+Vert+(Base/2))+100, Base+100);
for y = (Base/2):Vert+(Base/2)  % Rectangle
    for x = 0:20
        edge_x = int16(x)+50;
        edge_y = int16(y)+50;
        edge_map(edge_y, edge_x) = map(edge_y, edge_x);
    end
    for x = Base-20:Base
        edge_x = int16(x)+50;
        edge_y = int16(y)+50;
        edge_map(edge_y, edge_x) = map(edge_y, edge_x);
    end
end
if(deg == 180)
    for x = 0:Base
        for y = Vert+(Base/2)-20:Vert+(Base/2)
            edge_x = int16(x)+50;
            edge_y = int16(y)+50;
            edge_map(edge_y, edge_x) = map(edge_y, edge_x);
        end
    end
end


d = 20/cos(R);
for x = 0:(Base/2)  % Left slope of triangle
    for y = slope*x+Vert+(Base/2)-d:slope*x+Vert+(Base/2)
        edge_x = int16(x)+50;
        edge_y = int16(y)+50;
        if(edge_y > Height+Vert+(Base/2)+100)
            edge_y = int16(Height+Vert+(Base/2))+100;
        end
        edge_map(edge_y, edge_x) = map(edge_y, edge_x);
    end
end
for x = (Base/2):Base % Right slope of triangle
    for y = -slope*x+(Height*2)+Vert+(Base/2)-d:-slope*x+(Height*2)+Vert+(Base/2)
        edge_x = int16(x)+50;
        edge_y = int16(y)+50;
        if(edge_y > Height+Vert+(Base/2)+100)
            edge_y = int16(Height+Vert+(Base/2))+100;
        end
        edge_map(edge_y, edge_x) = map(edge_y, edge_x);
    end
end

for C = -pi:0.01:0
    for j = 0:20*cos(C)
        edge_x = int16((Base/2)*cos(C)+(Base/2))+50-j;
        for k = 0:20*(-sin(C))
            edge_y = int16((Base/2)*sin(C)+(Base/2))+50+k;

            if(edge_x > Base+100)
                edge_x = int16(Base)+100;
            end
            if(edge_y > (Base/2)+50)
                edge_y = int16((Base/2))+50;
            end
            edge_map(edge_y, edge_x) = map(edge_y, edge_x);
        end
    end
    for j = 0:20*(-cos(C))
        edge_x = int16((Base/2)*cos(C)+(Base/2))+50+j;
        for k = 0:20*(-sin(C))
            edge_y = int16((Base/2)*sin(C)+(Base/2))+50+k;

            if(edge_x > Base+100)
                edge_x = int16(Base)+100;
            end
            if(edge_y > (Base/2)+50)
                edge_y = int16((Base/2))+50;
            end
            edge_map(edge_y, edge_x) = map(edge_y, edge_x);
        end
    end
end
figure(3)
imagesc(-50:Base+50, -50:int16(Height+Vert+(Base/2))+50, edge_map);
colormap jet;
colorbar;
axis xy;
daspect([1 1 1]);

% F-actin intensity within a 2 μm-wide and 5 μm-long regions in both sides of the vertex
corner_map = zeros(int16(Height+Vert+(Base/2))+100, Base+100);
n_1 = (Base/2)-(50*cos(R));
m_1 = Height+Vert+(Base/2)-(50*sin(R));
b_1 = m_1 - (-1/slope)*n_1;
for x = (Base/2)-(50*cos(R)):(Base/2)
    for y = slope*x+Vert+(Base/2)-d:slope*x+Vert+(Base/2)
        if(y > -1/slope*x+b_1)
            edge_x = int16(x)+50;
            edge_y = int16(y)+50;
        end
        if(edge_y > Height+Vert+(Base/2)+100)
            edge_y = int16(Height+Vert+(Base/2))+100;
        end
        corner_map(edge_y, edge_x) = edge_map(edge_y, edge_x);
    end
end
n_2 = (Base/2)+(50*cos(R));
m_2 = Height+Vert+(Base/2)-(50*sin(R));
b_2 = m_2 - (1/slope*n_2);
for x = (Base/2):(Base/2)+(50*cos(R))
    for y = -slope*x+(Height*2)+Vert+(Base/2)-d:-slope*x+(Height*2)+Vert+(Base/2)
        if(y > 1/slope*x+b_2)
            edge_x = int16(x)+50;
            edge_y = int16(y)+50;
        end
        if(edge_y > Height+Vert+(Base/2)+100)
            edge_y = int16(Height+Vert+(Base/2))+100;
        end
        corner_map(edge_y, edge_x) = edge_map(edge_y, edge_x);
    end
end
if(deg == 180)
    for x = (Base/2)-50:(Base/2)+50
        for y = (Base/2)+Vert-20:(Base/2)+Vert
            edge_x = int16(x)+50;
            edge_y = int16(y)+50;
            corner_map(edge_y, edge_x) = edge_map(edge_y, edge_x);
        end
    end
end
figure(4)
imagesc(-50:Base+50, -50:int16(Height+Vert+(Base/2))+50, corner_map);
colormap jet;
colorbar;
axis xy;
daspect([1 1 1]);

% Area of ​​whole cell, periphery, and corner
body_area = ((Base*Height)/2) + (Base*Vert) + ((Base/2)*(Base/2)*pi/2);
edge_base = Base - ((20/cos((pi/2)-(R)))*2);
edge_height = edge_base/2*tan(R);
edge_area = body_area - (edge_base*edge_height/2) - ((Base-40) * Vert) - (((Base/2)-20)*((Base/2)-20)*pi/2);
if(deg == 180)
    edge_area = body_area - ((Base-40)*(Vert-20)) - (((Base/2)-20)*((Base/2)-20)*pi/2);
end

corner_area = (20*50*2) - (d*(d/2)/tan(R)*2/2) - (20*20*tan((90-deg)*(pi/180)));
if(deg == 180)
    corner_area = 20*100;
end

body_intensity = sum(map, "all")/body_area;
edge_intensity = sum(edge_map, "all")/edge_area;
corner_intensity = sum(corner_map, "all")/corner_area;

disp("edge/body: " + (edge_intensity/body_intensity))
disp("corner/edge: " + (corner_intensity/edge_intensity))

