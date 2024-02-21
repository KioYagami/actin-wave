clear variables;
close all;
v = 0.35;       %アクチン波の速度 %WT
%v = 0.20;       %アクチン波の速度 %Shootin1b-KO
length = 41;    %アクチン波の長さ
time = 20000;   %実行時間

pcell.Base = 350;  %三角形の底辺の長さ
pcell.deg = 60;    %二等辺三角形の頂点の角度
pcell.R = (90-(pcell.deg/2))*(pi/180);  %底角の角度(ラジアン)
pcell.Height = (pcell.Base/2)*tan(pcell.R);  %三角形の高さ 
pcell.slope = pcell.Height/(pcell.Base/2);   %三角形の辺の傾き

figure(1)
DrawCellRegion(pcell)

x = zeros(1, 1000);      %事前割り当て
y = zeros(1, 1000);
prex = zeros(1, 1000);      %事前割り当て
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
i = 0;        %アクチン波の初期個数は0個
veerFlag = 0;
LifetimeFlag = 0;
d = 1;
f = 1;
map = zeros(int16(pcell.Height)+100, pcell.Base+100);

%%
current_time = datetime('now');
% 現在時刻をシード値として設定
seed_value = current_time.Second; % 例: 現在の秒をシード値として使用
% 乱数ジェネレーターのシードを設定
rng(seed_value);

for t = 1:time          %時間

  generate = rand();

  if(generate > 0.9)         %アクチン波の発生確率
    i = i+1;                %アクチン波の個数

    if(LifetimeFlag == 0)
      p = rand();
      q = rand();
      m1 = min(p, q);
      m2 = 1.0 - max(p, q);
      m3 = max(p, q) - min(p, q);
      x(i) = m1*0 + m2*(pcell.Base/2) + m3*pcell.Base; %x座標
      y(i) = m1*0 + m2*pcell.Height + m3*0;      %y座標
    elseif(LifetimeFlag == 1)
      x(i) = x_dis(d);
      y(i) = y_dis(d);
      d = d + 1;
      if(d >= f)
        LifetimeFlag = 0;
      end
    end

    rad(i) = 2*pi*rand();       %移動方向
    vec_x(i) = v*cos(rad(i));   %xベクトル
    vec_y(i) = v*sin(rad(i));   %yベクトル
    Lifetime(i) = gamrnd(3.016848, 118.9898); %アクチン波の寿命 %WT
    %Lifetime(i) = gamrnd(2.397247, 106.2886); %アクチン波の寿命 %Shootin1b-KO
    generation_time(i) = t; %生成された時の時間

  end

  if(i>0)
    for num = 1:i   %アクチン波の個数だけ繰り返す

      if Lifetime(num) > 0
        status(num) = 0;  % Edge 上 = 1
        if(mod((t - generation_time(num)), 50) == 0)
          veer = (11.1138*randn()-0.614419)*(pi/180);     %変針
        else
          veer = 0;
        end

        rad(num) = rad(num) + veer;
        vec_x(num) = v*cos(rad(num));   %xベクトル
        vec_y(num) = v*sin(rad(num));   %yベクトル

        x(num) = x(num) + vec_x(num);   %x方向の移動
        y(num) = y(num) + vec_y(num);   %y方向の移動

        % エッジに来たときの処理
        if(y(num) >= pcell.slope*x(num))   %左斜辺に来た時
          b = y(num) - (-1/pcell.slope)*x(num); %左斜面に垂直でx(num),y(num)を通る直線
          x(num) = b/(pcell.slope-(-1/pcell.slope));  %交点のx座標
          y(num) = pcell.slope*x(num);          %交点のy座標
          status(num) = 1;
        elseif(y(num) >= -pcell.slope*x(num)+(pcell.Height*2))  %右斜辺に来た時
          b = y(num) - (1/pcell.slope)*x(num); %右斜面に垂直でx(num),y(num)を通る直線
          x(num) = (b-pcell.Height*2)/(-pcell.slope-(1/pcell.slope));  %交点のx座標
          y(num) = -pcell.slope*x(num)+(pcell.Height*2);         %交点のy座標
          status(num) = 1;
        end

        if(y(num) <= 0) %底面に来た時
          y(num) = y(num) - vec_y(num);
          status(num) = 1;
        end

        % エッジに沿って動いて頂点にきたとき
        % 上
        if(y(num) >= pcell.Height) %頂点より上に来た時
          y(num) = y(num) - vec_y(num);
          status(num) = 1;
        end
        % 左下
        if(x(num) <= 0)
          x(num) = x(num) - vec_x(num);
          status(num) = 1;
        end
        % 右下
        if(x(num) >= pcell.Base)
          x(num) = x(num) - vec_x(num);
          status(num) = 1;
        end

        % 寿命減らし
        Lifetime(num) = Lifetime(num) - 1;  %Lifetime減算
        if(Lifetime(num) <= 0)  %消失時の座標を記録
          x_dis(f) = x(num);
          y_dis(f) = y(num);
          LifetimeFlag = 1;
          f = f + 1;
          status(num) = 2;
        end

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

%Actin intensity map
%細胞外周幅20umのactin intensity
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
%頂点の角から両方向に50um, 幅20umの範囲のactin intensity
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

%細胞全体、外周、角の範囲の面積
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
X1 = 0:pcell.Base;    %三角形の底辺
Y1 = 0*X1;
hold on
plot(X1, Y1);

X2 = 0:(pcell.Base/2);     %三角形の左斜辺
Y2 = pcell.slope*X2;
hold on
plot(X2, Y2);

X3 = (pcell.Base/2):pcell.Base;      %三角形の右斜辺
Y3 = -pcell.slope*X3+(pcell.Height*2);
hold on
plot(X3, Y3);

xlim([-50 pcell.Base+50])
ylim([-50 pcell.Height+50])

daspect([1 1 1])

end

