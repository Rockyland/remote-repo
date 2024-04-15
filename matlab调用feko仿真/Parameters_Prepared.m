%% 所有图片统一归一化！
%  本代码主要包含两部分：1、图片各自强度归一化 2、图片强度统一归一化
%  可以用subplot画连续的几个图，并且所有图片归一化了，但是也没有距离徙动校正。
%  在主程序中，调用的还是ImageOutput函数,该代码主要用于调试。
% RCS_theta=RCS_theta(:,1:2000);
close all;clear all;
load RCS_theta_YZ02.mat

TorP='ETHETA';
c=3e8;%光速
C_freq=220e9;%中心频率
BW=10e9;%带宽
OriPhi=0;%自定义
Phi_Step=0.02;%方位角间隔
B=BW;
f0=C_freq;
lambda=c/f0;
a_num=size(RCS_theta,1); %画图时，方位向的长度，默认和距离向频点数保持一致
Phi_N=size(RCS_theta,2);
Omiga=Phi_Step*a_num;
theta=deg2rad(Omiga);
range_a=lambda/2/theta;
range_r=c/2/B;
kn=range_r;
counts=0;%for循环计数用的
fig = figure(2);%gif动图制作需要
% axis([0 101 0 101]);
% a_gap=102;\\

subplot_n=3; %画图时，默认凑成3*3的subplot,可以根据需求调整
subplot_nums=1;%用来判断是否用subplot来展示成像图片。因为subplot一般最多看9张。
if subplot_nums>1
    a_gap=floor(Phi_N/(subplot_n*subplot_n)); %默认凑成3*3的subplot,可以根据需求调整
    % a_gap=1500;
else
    a_gap=300;
end

k_start=1; %画图时，方位向的起点
k_end=Phi_N; %给最后一个子图保留位置
ColorRange=-30; %caxis函数的第一个参数值

flag1=max(a_gap,a_num);%仅仅为了防止索引超出数组边界
% 1、图片各自归一化

% axis([0 101 0 101])
for k=k_start:a_gap:k_end-flag1
    S=hamming(size(RCS_theta,1))*hamming(size(RCS_theta(:,k:k+a_num),2))'.*RCS_theta(:,k:k+a_num);
    
    Im=fty(ftx(S));
    G=mag2db(abs(Im)/max(max(abs(Im)))); 
    gg=G;%保存一个副本,用来寻找最大的n个像素的位置
    mm(k)=max(max(abs(Im)));%记录每个图片的最强幅度
    [N,M]=size(G);

    counts=counts+1;
    A=max(max(gg));%最大值
    a=min(min(gg));%最小值
    % a = sort(G(:));%按列排序
    % [row1,col1]=find(G>=a(end-5));
    %由于散射点会有展宽，因此，当确定了最大值后，其周围的八个点也都跟它一起被筛除

    % plot(row1,col1,'o-')
    % hold on
    
    %寻找最大的n个点 
    for ii=1:1
        [row(counts),col(counts)]=find(gg==A);
        % scatter(row(counts),col(counts),'red','filled');
        indx(ii,:)=[row(counts),col(counts)];
        % plot(row(counts),col(counts),'o-')
        % line(row(counts),col(counts))
        gg(row(counts)-3:row(counts)+3,col(counts)-3:col(counts)+3)=a;
        A=max(max(gg));%最大值
        % line(row1,col1)
        % axis([0 101 0 101])
        
    end
    hold on
    plot(indx(:,1),indx(:,2),'-o')
    disp(k*0.02)

    % if subplot_nums<=1
    %     figure
    % else
    %     subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
    % end
    % figure(2)
    % imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
    

    % imagesc(G)
    % colorbar
    % colormap hot
    % caxis([ColorRange,0])
    % xlabel('Azimuth(m)');
    % ylabel('Range(m)');
    title(['ETHETA  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
    % 
    pause(0.5)
    drawnow
    frame = getframe(fig);
    im{counts} = frame2im(frame);
    % ftresize(18)
end
filename = 'test222.gif'; % Specify the output file name
for idx = 1:counts
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
    end
end

%Keystone变换之后，再成像
figure
for k=k_start:a_gap:k_end-flag1

    S=hamming(size(RCS_theta,1))*hamming(size(RCS_theta(:,k:k+a_num),2))'.*RCS_theta(:,k:k+a_num);


    S_Keystone=Keystone(S,B,f0);

     % Keystone之后，做补偿，再成像
    hrrp=ftx(S_Keystone);
%         [N,M]=size(phi_rcs);
    w=deg2rad(Omiga)^2;
    y0=0;
    PRF=M;
    alpha=2*pi*f0/(c*PRF^2);
    phase=exp(1i*alpha*([-N/2:N/2-1].'*kn+y0)*w*([-M/2:M/2-1].^2));
    hrrp_phasecomp=hrrp.*phase;
%         hrrp_phasecomp=ftx(hamming(1024)*hamming(1024)'.*iftx(hrrp_phasecomp),4096);
%         IM=abs(fty(hrrp_phasecomp,4096));
    IM=abs(fty(hrrp_phasecomp));

    IM=mag2db(IM/max(max(IM)));

    if subplot_nums<=1
        figure
    else
        subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
    end

    imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,IM)
    % imagesc(G)
    colorbar
    colormap hot
    caxis([ColorRange,0])
    xlabel('Azimuth(m)');
    ylabel('Range(m)');
    % title([TorP,'  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°',' 相位补偿后']);
    title(['phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
end


MM=max(mm);%全局最强幅度
figure;
% 2、图片统一归一化
for k=k_start:a_gap:k_end-flag1
    S=hamming(size(RCS_theta,1))*hamming(size(RCS_theta(:,k:k+a_num),2))'.*RCS_theta(:,k:k+a_num);
    
    Im=fty(ftx(S));
    G=mag2db(abs(Im)/MM); %所有成像图片强度统一归一化
    [N,M]=size(G);
    if subplot_nums>1
        subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
    else
        fig=figure;
        movegui(fig,[200 200]);
    end
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
caxis([ColorRange,0])
set(h, 'Position', [0.92 0.1 0.020 0.8]); % 调整共用colorbar的位置和尺寸


figure;plot((1:length(RCS_theta(1,:)))*0.02,max(abs(RCS_theta)))
xlabel('方位角/°');ylabel('幅度');

% figure;plot((1:length(RCS_V(1,:)))*0.02,max(abs(RCS_V)))
% xlabel('方位角/°');ylabel('幅度');
% close all
% for jj=1:10:100
%     zz=RCS_theta(jj,:);

%     figure;plot((1:9000)*0.02,abs(zz))
%     xlabel('方位角/°')
%     ylabel('幅度')
%     % axis([-100 9100 -10 40])
% end
% close all
% for jj=1:1000:9000
%     yy=RCS_theta(:,jj);
%     figure;plot(1:length(yy),abs(yy),'-ro')
%     % xlabel('方位角/°')
%     % ylabel('幅度')
%     % axis([-10 120 -10 40])
% end


