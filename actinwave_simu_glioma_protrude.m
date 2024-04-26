clear variables;
close all;
v = 0.33;       % Actin wave velocity % WT
%v = 0.20;      % Actin wave velocity % Shootin1b-KO#1
length = 41;    % Length of actin waves
time = 20000;   % Simulation time

Base = 400;  % Base length of triangle (0.1 μm/pixel)
deg = 180;    % Vertex angle of triangle

R = (90-(deg/2))*(pi/180); % Base angle of triangle (radian)
Height = (Base/2)*tan(R); 
Area = (Base*Height)/2;
slope = Height/(Base/2);

figure(1)
C = -pi:0.01:0;
plot((Base/2)*cos(C)+(Base/2), (Base/2)*sin(C)+(Base/2))

C1 = -pi:0.01:pi;
x_1 = Base/10;
theta = acos(((Base/2)-x_1)/(Base/2));
y_1 = (Base/2) + (Base/2)*sin(-theta);
hold on
plot(20*cos(C1)+(x_1), 20*sin(C1)+(y_1))

x_2 = Base/4;
theta = acos(((Base/2)-x_2)/(Base/2));
y_2 = (Base/2) + (Base/2)*sin(-theta);
hold on
plot(20*cos(C1)+(x_2), 20*sin(C1)+(y_2))

x_3 = Base/4*3;
theta = acos(((Base/2)-x_3)/(Base/2));
y_3 = (Base/2) + (Base/2)*sin(-theta);
hold on
plot(20*cos(C1)+(x_3), 20*sin(C1)+(y_3))

x_4 = Base/10*9;
theta = acos(((Base/2)-x_4)/(Base/2));
y_4 = (Base/2) + (Base/2)*sin(-theta);
hold on
plot(20*cos(C1)+(x_4), 20*sin(C1)+(y_4))

Vert = 0;
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
            
                if(x(num)>=(x_1-20) && x(num)<=x_1)
                    if(((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2) && ((x(num)-x_1)^2+(y(num)-y_1)^2 >= (20)^2))
                        theta = atan(((y_1)-y(num))/((x_1)-x(num))); 
                        x(num) = (20)*cos(-pi+theta)+(x_1);
                        y(num) = (20)*sin(-pi+theta)+(y_1);
                        status(num) = 1;
                    end
                elseif(x(num)>=x_1 && x(num)<=(x_1+20))
                    if(((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2) && ((x(num)-x_1)^2+(y(num)-y_1)^2 >= (20)^2))
                        theta = atan(((y_1)-y(num))/(x(num)-(x_1)));
                        x(num) = (20)*cos(-theta)+(x_1);
                        y(num) = (20)*sin(-theta)+(y_1);
                        status(num) = 1;
                    end
                elseif(x(num)>=(x_2-20) && x(num)<=x_2)
                    if(((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2) && ((x(num)-x_2)^2+(y(num)-y_2)^2 >= (20)^2))
                        theta = atan(((y_2)-y(num))/((x_2)-x(num))); 
                        x(num) = (20)*cos(-pi+theta)+(x_2);
                        y(num) = (20)*sin(-pi+theta)+(y_2);
                        status(num) = 1;
                    end
                elseif(x(num)>=x_2 && x(num)<=(x_2+20))
                    if(((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2) && ((x(num)-x_2)^2+(y(num)-y_2)^2 >= (20)^2))
                        theta = atan(((y_2)-y(num))/(x(num)-(x_2)));
                        x(num) = (20)*cos(-theta)+(x_2);
                        y(num) = (20)*sin(-theta)+(y_2);
                        status(num) = 1;
                    end
                elseif(x(num)>=(x_3-20) && x(num)<=x_3)
                    if(((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2) && ((x(num)-x_3)^2+(y(num)-y_3)^2 >= (20)^2))
                        theta = atan(((y_3)-y(num))/((x_3)-x(num))); 
                        x(num) = (20)*cos(-pi+theta)+(x_3);
                        y(num) = (20)*sin(-pi+theta)+(y_3);
                        status(num) = 1;
                    end
                elseif(x(num)>=x_3 && x(num)<=(x_3+20))
                    if(((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2) && ((x(num)-x_3)^2+(y(num)-y_3)^2 >= (20)^2))
                        theta = atan(((y_3)-y(num))/(x(num)-(x_3)));
                        x(num) = (20)*cos(-theta)+(x_3);
                        y(num) = (20)*sin(-theta)+(y_3);
                        status(num) = 1;
                    end
                elseif(x(num)>=(x_4-20) && x(num)<=x_4)
                    if(((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2) && ((x(num)-x_4)^2+(y(num)-y_4)^2 >= (20)^2))
                        theta = atan(((y_4)-y(num))/((x_4)-x(num))); 
                        x(num) = (20)*cos(-pi+theta)+(x_4);
                        y(num) = (20)*sin(-pi+theta)+(y_4);
                        status(num) = 1;
                    end
                elseif(x(num)>=x_4 && x(num)<=(x_4+20))
                    if(((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2) && ((x(num)-x_4)^2+(y(num)-y_4)^2 >= (20)^2))
                        theta = atan(((y_4)-y(num))/(x(num)-(x_4)));
                        x(num) = (20)*cos(-theta)+(x_4);
                        y(num) = (20)*sin(-theta)+(y_4);
                        status(num) = 1;
                    end
                else
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

% F-actin intensity 
protrusion_map = zeros(int16(Height+Vert+(Base/2))+100, Base+100);

for x = 0:Base+50
    for y = 0:(Base/2)+50
        if(x>=(x_1-20) && x<=(x_1+20))
            if(((x-(Base/2))^2+(y-(Base/2))^2 >= (Base/2)^2))
                edge_x = int16(x)+50;
                edge_y = int16(y)+50;
                protrusion_map(edge_y, edge_x) = map(edge_y, edge_x);
            end
        elseif(x>=(x_2-20) && x<=(x_2+20))
            if(((x-(Base/2))^2+(y-(Base/2))^2 >= (Base/2)^2))
                edge_x = int16(x)+50;
                edge_y = int16(y)+50;
                protrusion_map(edge_y, edge_x) = map(edge_y, edge_x);
            end
        elseif(x>=(x_3-20) && x<=(x_3+20))
            if(((x-(Base/2))^2+(y-(Base/2))^2 >= (Base/2)^2))
                edge_x = int16(x)+50;
                edge_y = int16(y)+50;
                protrusion_map(edge_y, edge_x) = map(edge_y, edge_x);
            end
        elseif(x>=(x_4-20) && x<=(x_4+20))
            if(((x-(Base/2))^2+(y-(Base/2))^2 >= (Base/2)^2))
                edge_x = int16(x)+50;
                edge_y = int16(y)+50;
                protrusion_map(edge_y, edge_x) = map(edge_y, edge_x);
            end
        end
    end
end

figure(3)
imagesc(-50:Base+50, -50:int16(Height+Vert+(Base/2))+50, protrusion_map);
colormap jet;
colorbar;
axis xy;
daspect([1 1 1]);

% Area of ​​whole cell, periphery, and corner
body_area = nnz(map);
protrusion_area = nnz(protrusion_map);

all_intensity = sum(map, "all");
protrusion_intensity = sum(protrusion_map, "all");

disp("protrusion/body: " + ((protrusion_intensity/protrusion_area)/(all_intensity/body_area)))
