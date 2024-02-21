clear variables;
close all;
v = 0.33;       %アクチン波の速度 %WT
%v = 0.20;      %アクチン波の速度 %Shootin1b-KO
length = 41;    %アクチン波の長さ
time = 5700;    %実行時間

pa.Base = 238;  %三角形の底辺の長さ
pa.deg = 120;    %二等辺三角形の頂点の角度
pa.R = (90-(pa.deg/2))*(pi/180);  %底角の角度(ラジアン)
pa.Height = (pa.Base/2)*tan(pa.R);  %三角形の高さ
pa.slope = pa.Height/(pa.Base/2);   %三角形の辺の傾き
pa.Area = (pa.Base*pa.Height)/2;    %三角形の面積
pa.Vert = (30725-pa.Area)/pa.Base;    %四角形の縦の長さ

C = -pi:0.01:0;     %半円を描く
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

% 濃度表示の設定
%[x_grid, y_grid] = meshgrid(,5);
sigma = 15;
kern = @(x) exp(- 0.5 * (x/sigma).^2);
total_I = zeros(size(pa.xgrid));

current_time = datetime('now');
% 現在時刻をシード値として設定
seed_value = current_time.Second; % 例: 現在の秒をシード値として使用
% 乱数ジェネレーターのシードを設定
rng(seed_value);

outputVideo = VideoWriter('output_video.avi'); % 動画ファイル名を指定
outputVideo.FrameRate = 60; % フレームレートを設定
open(outputVideo);

for t = 1:time          %時間

    generate = rand();

    if(generate > 0.8)         %アクチン波の発生確率
        i = i+1;                %アクチン波の個数
        if(LifetimeFlag == 0)
            x(i) = rand()*pa.Base;     %x座標
            if(x(i) <= (pa.Base/2))                %左辺側
                C1 = -acos((x(i)-(pa.Base/2))/(pa.Base/2));
                y(i) = rand()*(pa.slope*x(i)+pa.Vert+(pa.Base/2)-((pa.Base/2)*sin(C1)+(pa.Base/2)))+((pa.Base/2)*sin(C1)+(pa.Base/2));     %y座標
            end
            if(x(i) > (pa.Base/2))                 %右辺側
                C1 = -acos((x(i)-(pa.Base/2))/(pa.Base/2));
                y(i) = rand()*(-pa.slope*x(i)+pa.Vert+(pa.Base/2)+(pa.Height*2)-((pa.Base/2)*sin(C1)+(pa.Base/2)))+((pa.Base/2)*sin(C1)+(pa.Base/2));    %y座標
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
                if(x(num)>=pa.Base)  %右直線に来た時
                    x(num) = x(num) - vec_x(num);
                    status(num) = 1;
                end

                if(y(num)>=pa.slope*x(num)+pa.Vert+(pa.Base/2) && x(num)<=(pa.Base/2) && pa.deg~=180)   %左斜辺に来た時
                    b = y(num) - (-1/pa.slope)*x(num); %左斜面に垂直でx(num),y(num)を通る直線
                    x(num) = (b-pa.Vert-(pa.Base/2))/(pa.slope-(-1/pa.slope));  %交点のx座標
                    y(num) = pa.slope*x(num)+pa.Vert+(pa.Base/2);          %交点のy座標
                    status(num) = 1;
                elseif(y(num) >= -pa.slope*x(num)+(pa.Height*2)+pa.Vert+(pa.Base/2) && x(num)>(pa.Base/2) && pa.deg~=180)  %右斜辺に来た時
                    b = y(num) - (1/pa.slope)*x(num); %右斜面に垂直でx(num),y(num)を通る直線
                    x(num) = (b-pa.Height*2-pa.Vert-(pa.Base/2))/(-pa.slope-(1/pa.slope));  %交点のx座標
                    y(num) = -pa.slope*x(num)+(pa.Height*2)+pa.Vert+(pa.Base/2);         %交点のy座標
                    status(num) = 1;
                end

                if((y(num)<=(pa.Base/2)) && (x(num)<=(pa.Base/2)) && ((x(num)-(pa.Base/2))^2+(y(num)-(pa.Base/2))^2 >= (pa.Base/2)^2))   %左側の円上に来た時
                    theta = atan(((pa.Base/2)-y(num))/((pa.Base/2)-x(num)));  %円の中心に対する角度を求める
                    x(num) = (pa.Base/2)*cos(-pi+theta)+(pa.Base/2);  %円上に移動
                    y(num) = (pa.Base/2)*sin(-pi+theta)+(pa.Base/2);  %円上に移動
                    status(num) = 1;
                elseif((y(num)<=(pa.Base/2)) && (x(num)>=(pa.Base/2)) && ((x(num)-(pa.Base/2))^2+(y(num)-(pa.Base/2))^2 >= (pa.Base/2)^2))   %右側の円上に来た時
                    theta = atan(((pa.Base/2)-y(num))/(x(num)-(pa.Base/2)));
                    x(num) = (pa.Base/2)*cos(-theta)+(pa.Base/2);
                    y(num) = (pa.Base/2)*sin(-theta)+(pa.Base/2);
                    status(num) = 1;
                end

                if(y(num)>=pa.Height+pa.Vert+(pa.Base/2))
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

            end % Lifetime(num) > 0
        end % num = 1:i

        % アクチン波の描画
        if(mod(t,10)==0)

            idx = find(Lifetime > 0);

            tail_x = length*cos(rad)/2;
            tail_y = length*sin(rad)/2;

            % Lifetime が 0 のアクチンの座標を NaN で埋める
            nan_idx = find(Lifetime <= 0);
            x(nan_idx) = NaN;
            y(nan_idx) = NaN;
            tail_x(nan_idx) = NaN;
            tail_y(nan_idx) = NaN;

            %濃度計算
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
            colormap('jet'); % カラーマップを'jet'に設定
            total_I = total_I + I;
            hold on;

            % エッジ表示
            patch(pa.x, pa.y, 'w', 'EdgeColor', 'w', 'FaceColor', 'none');

            % actin wave 描画
            x_ = x(idx);
            y_ = y(idx);

            h = quiver(x_, y_, tail_x(idx), tail_y(idx), 'AutoScale', 'off', 'LineWidth', 1.5, 'Color', 'magenta');
            [h.ShowArrowHead] = deal('off');

            h = quiver(x_, y_, -tail_x(idx), -tail_y(idx), 'AutoScale', 'off', 'LineWidth', 1.5, 'Color', 'magenta');
            [h.ShowArrowHead] = deal('off');

            % 時間とかの表示
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
            colormap('jet'); % カラーマップを'jet'に設定
            xlim(pa.xl);
            ylim(pa.yl);
            axis off; xticks([]); yticks([]);
            title('Average intensity');
            daspect([1 1 1]);

            drawnow;

            frame = getframe(gcf);
            writeVideo(outputVideo, frame);

            clf;
        end % 各アクチンの描画

    end % i > 0
end % t = 1:time

close(outputVideo);
