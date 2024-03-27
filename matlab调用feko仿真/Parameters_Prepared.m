%% 所有图片统一归一化！
%  本代码主要包含两部分：1、图片各自强度归一化 2、图片强度统一归一化
%  可以用subplot画连续的几个图，并且所有图片归一化了，但是也没有距离徙动校正。
%  在主程序中，调用的还是ImageOutput函数,该代码主要用于调试。

close all;

TorP='ETHETA';
c=3e8;%光速
C_freq=220e9;%中心频率

BW=10e9;%带宽
OriPhi=0;%自定义
Phi_Step=0.02;%方位角间隔

Phi_N=9000;
subplot_n=3; %画图时，默认凑成3*3的subplot,可以根据需求调整
a_num=101; %画图时，方位向的长度，默认和距离向频点数保持一致
B=BW;
f0=C_freq;
lambda=c/f0;
Omiga=Phi_Step*a_num;
theta=deg2rad(Omiga);
range_a=lambda/2/theta;
range_r=c/2/B;
% a_gap=102;
a_gap=floor(Phi_N/(subplot_n*subplot_n)); %默认凑成3*3的subplot,可以根据需求调整

k_start=1; %画图时，方位向的起点
k_end=Phi_N; %给最后一个子图保留位置
ColorRange=-40; %caxis函数的第一个参数值

% 1、图片各自归一化
figure;
for k=k_start:a_gap:k_end
    S=hamming(size(RCS_theta,1))*hamming(size(RCS_theta(:,k:k+a_num),2))'.*RCS_theta(:,k:k+a_num);
    
    Im=fty(ftx(S));
    G=mag2db(abs(Im)/max(max(abs(Im)))); 
    mm(k)=max(max(abs(Im)));%记录每个图片的最强幅度
    [N,M]=size(G);
    % fig=figure;
    % movegui(fig,[200 200]);
    subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
    imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
    colorbar
    colormap hot
    caxis([ColorRange,0])
    xlabel('Azimuth(m)');
    ylabel('Range(m)');
    title(['ETHETA  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
    % ftresize(18)
end

MM=max(mm);%全局最强幅度
figure;
% 2、图片统一归一化
for k=k_start:a_gap:k_end
    S=hamming(size(RCS_theta,1))*hamming(size(RCS_theta(:,k:k+a_num),2))'.*RCS_theta(:,k:k+a_num);
    
    Im=fty(ftx(S));
    G=mag2db(abs(Im)/MM); %所有成像图片强度统一归一化
    [N,M]=size(G);

    subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
    imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
    colormap hot
    caxis([ColorRange,0])
    xlabel('Azimuth(m)');
    ylabel('Range(m)');
    title(['ETHETA  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
    % ftresize(18)
end

% 添加统一的标题
axes('Position', [0, 0, 1, 1], 'Visible', 'off');
text(0.5, 0.97, '成像序列', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% 添加共用的colorbar
h = colorbar('location', 'eastoutside'); % 指定共用的colorbar位置，可根据需要调整
set(h, 'Position', [0.92 0.1 0.020 0.8]); % 调整共用colorbar的位置和尺寸


