
% FEKO软件的程序控制与序列程序
% 一定要打开cadfeko 才能运行 !!!
% 该版本新增了自动读取 频点数、初始theta角、初始phi角！！！由函数 Read_Original_Value实现。
% 确保该路径下只有一个.pre文件！！！

%% 要求本m文件与FEKO工程文件在同一路径下，然后在matlab中基于该路径运行m文件。
clear all
close all
clc
tic % 测试程序运行时间开始

%% 初始化，参数设置
Phi_N = 9000;% 方位角的个数A
Phi_Step = 0.02; %步进角度0.006°
c=3e8; % 光速

%% 确保该文件夹只有一个适用的.pre文件！！！

%% 第一部分：找到.pre文件、分别读取各自的关键初始值
filepath='E:\Satelites\zhangye_paowumian\目标特性\PaoWM\01';
cd(filepath);%切换到FEKO工程文件所在路径，这个一定要对应的替换哈！！！
feko_path = '"D:\Program Files\Altair\2021\feko\bin\runfeko"';%指定FEKO的安装路径，注意引号的使用，先是单引号，再加上双引号。这是处理路径存在空格的情况
A = dir(fullfile('','*.pre'));  %读取文件夹下所有的*.pre文件。A = dir(fullfile('','*.pre'));  A = dir(fullfile('','*.fek'))
prefile=strcat(filepath,'\',A.name);

[OriHz,S_freq,E_freq,OriTheta,OriPhi] = Read_Original_Value(prefile);

C_freq = (S_freq+E_freq)/2; %中心频率
BW = E_freq - S_freq; %带宽
RCS_theta=zeros(OriHz,Phi_N);  
RCS_phi=zeros(OriHz,Phi_N);
RCS_V=zeros(OriHz,Phi_N);

bar = waitbar(0,['初始phi角为',num2str(OriPhi),' 数据读取中...']);
ii = 1;
% ii=124
while ii<=Phi_N 

%% 对文件夹内的*.pre文件进行操作，默认只有一个.pre文件

    editfekoname=A.name;
    %%以下的并行计算与否以及并行核数的设置，下面默认的是调用本地全部核数，如果是希望2个核，即改np all为 np 2
%     d = strcat([feko_path,32,editfekoname]) ;% 不要用all 用2个或4个足矣
   d = strcat([feko_path,32,editfekoname,' -np 2 --parallel-authenticate localonly']);  % -np 4  d = strcat([feko_path,32,filename,' -np all']) 
dos(d); 

%% 第二部分：寻找、读取.out文件
lengthE = length(editfekoname);
lengthP = length('.pre');
tmp_name = editfekoname(1:lengthE-lengthP);
outname=strcat(tmp_name,'.out');

B = dir(fullfile('',outname));
filename_out=strcat(filepath,'\',B.name);

%读取.out文件分件
[RCS_complex_theta,RCS_complex_phi,RCS_value]=FekoReadFarfieldRCS_MonoStatic_OnlyForRCS(filename_out,OriHz);
%PlotRangeSequence(Phi_deg,Freq_Hz,RCS_complex_theta,1024);%用好用同极化数据(θ极化)


RCS_theta(:,ii)=RCS_complex_theta;
RCS_phi(:,ii)=RCS_complex_phi;
RCS_V(:,ii)=RCS_value;


%% 第三部分：修改.pre文件，并保存
% fid = fopen(editfekoname,'r');
CharPhi='A0: 0 :  : 1 : 1 : 1 : 1 : 0 : ';

Phi_new = OriPhi + Phi_Step*ii; % 待验证的角度

Phi_new_char = num2str(Phi_new); % 转换成字符串
new_contents_phi = strcat("A0: 0 :  : 1 : 1 : 1 : 1 : 0 : ",num2str(OriTheta)," : ",Phi_new_char," : 0 : 0 : ",num2str(Phi_Step)," : 0   ** PlaneWaveSource1");
new_contents_phi = new_contents_phi + newline;


text_modify(editfekoname,editfekoname,new_contents_phi,CharPhi)%读取文件夹中的文件，用新的内容new_contents替换 指定行的旧的内容old_contents

str_bar = ['初始phi角为',num2str(OriPhi),' 计算中...',num2str(100*ii/Phi_N),'%'];
waitbar(ii/Phi_N,bar,str_bar);
ii = ii + 1
end % while ii 结束
t = toc % 测试程序运行时间结尾，用于查看程序运行时间。

%% 保存回波数据
save RCS_theta.mat RCS_theta
save RCS_phi.mat RCS_phi
save RCS_V.mat RCS_V

%% 第四部分： 成像
%% 成像相关参数 初始化(为第四部分的成像做准备，与前三部分无关）
subplot_n=3; %画图时，默认凑成3*3的subplot,可以根据需求调整
a_num=OriHz; %画图时，方位向的长度，默认和距离向频点数保持一致
a_gap=floor(Phi_N/(subplot_n*subplot_n)); %默认凑成3*3的subplot,可以根据需求调整
k_start=1; %画图时，方位向的起点
k_end=Phi_N; %给最后一个子图保留位置
% k_end=Phi_N;
ColorRange=-50; %caxis函数的第一个参数值

load RCS_theta.mat
ImageOutput(subplot_n,a_num,a_gap,k_start,k_end,...
   RCS_theta,BW,C_freq,c,Phi_Step,OriPhi,ColorRange,'ETHETA')
% TorP='ETHETA' or 'EPHI', RCS_TorP 中的TorP必须与 TorP 对应
load RCS_phi.mat
ImageOutput(subplot_n,a_num,a_gap,k_start,k_end,...
   RCS_phi,BW,C_freq,c,Phi_Step,OriPhi,ColorRange,'EPHI')




% figure;
% for k=k_start:a_gap:k_end
% S=hamming(size(RCS_theta,1))*hamming(size(RCS_theta(:,k:k+a_num),2))'.*RCS_theta(:,k:k+a_num);
% 
% Im=fty(ftx(S));
% G=mag2db(abs(Im)/max(max(abs(Im))));
% [N,M]=size(G);
% B=BW;
% f0=C_freq;
% lambda=c/f0;
% Omiga=Phi_Step*a_num;
% theta=deg2rad(Omiga);
% range_a=lambda/2/theta;
% range_r=c/2/B;
% subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
% imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
% % imagesc(G)
% colorbar
% colormap hot
% caxis([-10,0])
% xlabel('Azimuth(m)');
% ylabel('Range(m)');
% % xlim([-0.3,0.3]);
% % ylim([-0.3,0.3]);
% title(['ETHETA  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
% % ftresize(18)
% end
% 
% 
% load RCS_phi.mat
% 
% 
% figure;
% for k=k_start:a_gap:k_end
% S=hamming(size(RCS_phi,1))*hamming(size(RCS_phi(:,k:k+a_num),2))'.*RCS_phi(:,k:k+a_num);
% 
% Im=fty(ftx(S));
% G=mag2db(abs(Im)/max(max(abs(Im))));
% [N,M]=size(G);
% B=BW;
% f0=C_freq;
% lambda=c/f0;
% Omiga=Phi_Step*a_num;
% theta=deg2rad(Omiga);
% range_a=lambda/2/theta;
% range_r=c/2/B;
% subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
% imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
% % imagesc(G)
% colorbar
% colormap hot
% caxis([-10,0])
% xlabel('Azimuth(m)');
% ylabel('Range(m)');
% % xlim([-0.3,0.3]);
% % ylim([-0.3,0.3]);
% title(['EPHI  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
% % ftresize(18)
% end





% S=hamming(size(RCS_phi,1))*hamming(size(RCS_phi,2))'.*RCS_phi;
% 
% Im=fty(ftx(S,size(RCS_phi,1)),size(RCS_phi,2));
% G=mag2db(abs(Im)/max(max(abs(Im))));
% figure;imagesc(G)
% colorbar
% colormap hot
% caxis([-50,0])

% [Im,y,x,Polar2Cartesian_data,fx,fy]=TurnImaging(Phi_deg,Freq_Hz,RCS_complex_theta,4,0,2048);

% RCS_complex_phi=load("zzz_RCS_complex_phi.mat");
% RCS_complex_phi=getfield(RCS_complex_phi,'RCS_complex_phi');

% Im=ftx(fty(RCS_phi(:,1:50)));
% G=20*log10(abs(Im));
% gm=max(max(G));
% gn=gm-60;
% G=255/(gm-gn)*(G-gn).*(G>gn);
% figure
% imagesc(G)
% colormap hot
% set( gca, 'XTick', [], 'YTick', [] );%不显示坐标轴XY
% colorbar;
% title('方位角-6°到6°');
% xlabel('方位向');ylabel('距离向');
% PlotImage(Im,-y,x,[],[],50);