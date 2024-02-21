clear variables;
close all;
v = 0.35;       %アクチン波の速度 %WT
%v = 0.20;      %アクチン波の速度 %Shootin1b-KO
length = 41;    %アクチン波の長さ
time = 20000;    %実行時間

Base = 238;  %三角形の底辺の長さ
deg = 60;    %二等辺三角形の頂点の角度

R = (90-(deg/2))*(pi/180); %二等辺三角形の2角
Height = (Base/2)*tan(R);  %三角形の高さ
Area = (Base*Height)/2;    %三角形の面積 
slope = Height/(Base/2);   %三角形の辺の傾き

figure(1)
C = -pi:0.01:0;     %半円を描く
plot((Base/2)*cos(C)+(Base/2), (Base/2)*sin(C)+(Base/2))

Vert = (30725-Area)/Base;    %四角形の縦の長さ
y = (Base/2):Vert+(Base/2);
x = 0*y;    %左側
hold on
plot(x, y);

Vert = (30725-Area)/Base;    %四角形の縦の長さ
y = (Base/2):Vert+(Base/2);
x = 0*y+Base;   %右側
hold on
plot(x, y);

x = 0:(Base/2);     %三角形の左斜辺
y = slope*x+Vert+(Base/2);
hold on
plot(x, y);

x = (Base/2):Base;      %三角形の右斜辺
y = -slope*x+(Height*2)+Vert+(Base/2);
hold on
plot(x, y);

xlim([-50 Base+50])
ylim([-50 Height+Vert+Base/2+50])
daspect([1 1 1])

x = zeros(1, 1000);      %事前割り当て
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
i = 0;        %アクチン波の初期個数は0個
veerFlag = 0;
LifetimeFlag = 0;
d = 1;
f = 1;
map = zeros(int16(Height+Vert+(Base/2))+100, Base+100);

for t = 1:time          %時間

    generate = rand();

    if(generate > 0.9)         %アクチン波の発生確率
        i = i+1;                %アクチン波の個数
        if(LifetimeFlag == 0)
            x(i) = rand()*Base;     %x座標
            if(x(i) <= (Base/2))                %左辺側
                C1 = -acos((x(i)-(Base/2))/(Base/2));
                y(i) = rand()*(slope*x(i)+Vert+(Base/2)-((Base/2)*sin(C1)+(Base/2)))+((Base/2)*sin(C1)+(Base/2));     %y座標
            end
            if(x(i) > (Base/2))                 %右辺側
                C1 = -acos((x(i)-(Base/2))/(Base/2));
                y(i) = rand()*(-slope*x(i)+Vert+(Base/2)+(Height*2)-((Base/2)*sin(C1)+(Base/2)))+((Base/2)*sin(C1)+(Base/2));    %y座標
            end
        elseif(LifetimeFlag == 1)
            x(i) = x_dis(d);
            y(i) = y_dis(d);
            d = d + 1;
            if(d >= f)
                LifetimeFlag = 0;
            end
        end
        rad(i) = 2*pi*rand();       %方向
        vec_x(i) = v*cos(rad(i));   %xベクトル
        vec_y(i) = v*sin(rad(i));   %yベクトル
        Lifetime(i) = gamrnd(3.016848, 118.9898);   %アクチン波の寿命 %WT
        %Lifetime(i) = gamrnd(2.397247, 106.2886);  %アクチン波の寿命 %Shootin1b-KO
        generation_time(i) = t; %生成された時の時間
    end
    if(i>0)
        for num = 1:i   %アクチン波の個数だけ繰り返す
            if Lifetime(num) > 0
                status(num) = 0;
                if(mod((t - generation_time(num)), 50) == 0)
                    veer = (11.1138*randn()-0.614419)*(pi/180);     %変針
                else
                    veer = 0;
                end
            
                rad(num) = rad(num) - veer;
                vec_x(num) = v*cos(rad(num));   %xベクトル
                vec_y(num) = v*sin(rad(num));   %yベクトル
                
                x(num) = x(num) + vec_x(num);   %x方向の移動
                y(num) = y(num) + vec_y(num);   %y方向の移動
            
                if(x(num)<=0)     %左直線に来た時
                    x(num) = x(num) - vec_x(num);
                    status(num) = 1;
                end
                if(x(num)>=Base)  %右直線に来た時
                    x(num) = x(num) - vec_x(num);
                    status(num) = 1;
                end
            
                if(y(num)>=slope*x(num)+Vert+(Base/2) && x(num)<=(Base/2) && deg~=180)   %左斜辺に来た時
                    b = y(num) - (-1/slope)*x(num); %左斜面に垂直でx(num),y(num)を通る直線
                    x(num) = (b-Vert-(Base/2))/(slope-(-1/slope));  %交点のx座標
                    y(num) = slope*x(num)+Vert+(Base/2);          %交点のy座標
                    status(num) = 1;
                elseif(y(num) >= -slope*x(num)+(Height*2)+Vert+(Base/2) && x(num)>(Base/2) && deg~=180)  %右斜辺に来た時
                    b = y(num) - (1/slope)*x(num); %右斜面に垂直でx(num),y(num)を通る直線
                    x(num) = (b-Height*2-Vert-(Base/2))/(-slope-(1/slope));  %交点のx座標
                    y(num) = -slope*x(num)+(Height*2)+Vert+(Base/2);         %交点のy座標
                    status(num) = 1;
                end
            
                if((y(num)<=(Base/2)) && (x(num)<=(Base/2)) && ((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2))   %左側の円上に来た時
                    theta = atan(((Base/2)-y(num))/((Base/2)-x(num)));  %円の中心に対する角度を求める
                    x(num) = (Base/2)*cos(-pi+theta)+(Base/2);  %円上に移動
                    y(num) = (Base/2)*sin(-pi+theta)+(Base/2);  %円上に移動
                    status(num) = 1;
                elseif((y(num)<=(Base/2)) && (x(num)>=(Base/2)) && ((x(num)-(Base/2))^2+(y(num)-(Base/2))^2 >= (Base/2)^2))   %右側の円上に来た時
                    theta = atan(((Base/2)-y(num))/(x(num)-(Base/2)));
                    x(num) = (Base/2)*cos(-theta)+(Base/2);
                    y(num) = (Base/2)*sin(-theta)+(Base/2);
                    status(num) = 1;
                end
            
                if(y(num)>=Height+Vert+(Base/2))
                    y(num) = y(num) - vec_y(num);
                    status(num) = 1;
                end
                
                Lifetime(num) = Lifetime(num) -1;
                if(Lifetime(num) <= 0)
                    x_dis(f) = x(num);
                    y_dis(f) = y(num);
                    LifetimeFlag = 1;
                    f = f + 1;
                    status(num) = 2;
                end
                
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

%Actin intensity map
%細胞外周幅20umのactin intensity
edge_map = zeros(int16(Height+Vert+(Base/2))+100, Base+100);
for y = (Base/2):Vert+(Base/2)  %四角形
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
for x = 0:(Base/2)  %三角形左側
    for y = slope*x+Vert+(Base/2)-d:slope*x+Vert+(Base/2)
        edge_x = int16(x)+50;
        edge_y = int16(y)+50;
        if(edge_y > Height+Vert+(Base/2)+100)
            edge_y = int16(Height+Vert+(Base/2))+100;
        end
        edge_map(edge_y, edge_x) = map(edge_y, edge_x);
    end
end
for x = (Base/2):Base %三角形右側
    for y = -slope*x+(Height*2)+Vert+(Base/2)-d:-slope*x+(Height*2)+Vert+(Base/2)
        edge_x = int16(x)+50;
        edge_y = int16(y)+50;
        if(edge_y > Height+Vert+(Base/2)+100)
            edge_y = int16(Height+Vert+(Base/2))+100;
        end
        edge_map(edge_y, edge_x) = map(edge_y, edge_x);
    end
end

for C = -pi:0.01:0     %半円を描く
    for j = 0:20*cos(C) %円右側
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
    for j = 0:20*(-cos(C))  %円左側
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

%頂点の角から両方向に50um, 幅20umの範囲のactin intensity
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

%細胞全体、外周、角の範囲の面積
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

