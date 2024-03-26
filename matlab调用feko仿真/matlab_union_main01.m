
% FEKO����ĳ�����������г���
% һ��Ҫ��cadfeko �������� !!!
% �ð汾�������Զ���ȡ Ƶ��������ʼtheta�ǡ���ʼphi�ǣ������ɺ��� Read_Original_Valueʵ�֡�
% ȷ����·����ֻ��һ��.pre�ļ�������

%% Ҫ��m�ļ���FEKO�����ļ���ͬһ·���£�Ȼ����matlab�л��ڸ�·������m�ļ���
clear all
close all
clc
tic % ���Գ�������ʱ�俪ʼ

%% ��ʼ������������
Phi_N = 9000;% ��λ�ǵĸ���A
Phi_Step = 0.02; %�����Ƕ�0.006��
c=3e8; % ����

%% ȷ�����ļ���ֻ��һ�����õ�.pre�ļ�������

%% ��һ���֣��ҵ�.pre�ļ����ֱ��ȡ���ԵĹؼ���ʼֵ
filepath='E:\Satelites\zhangye_paowumian\Ŀ������\PaoWM\01';
cd(filepath);%�л���FEKO�����ļ�����·�������һ��Ҫ��Ӧ���滻��������
feko_path = '"D:\Program Files\Altair\2021\feko\bin\runfeko"';%ָ��FEKO�İ�װ·����ע�����ŵ�ʹ�ã����ǵ����ţ��ټ���˫���š����Ǵ���·�����ڿո�����
A = dir(fullfile('','*.pre'));  %��ȡ�ļ��������е�*.pre�ļ���A = dir(fullfile('','*.pre'));  A = dir(fullfile('','*.fek'))
prefile=strcat(filepath,'\',A.name);

[OriHz,S_freq,E_freq,OriTheta,OriPhi] = Read_Original_Value(prefile);

C_freq = (S_freq+E_freq)/2; %����Ƶ��
BW = E_freq - S_freq; %����
RCS_theta=zeros(OriHz,Phi_N);  
RCS_phi=zeros(OriHz,Phi_N);
RCS_V=zeros(OriHz,Phi_N);

bar = waitbar(0,['��ʼphi��Ϊ',num2str(OriPhi),' ���ݶ�ȡ��...']);
ii = 1;
% ii=124
while ii<=Phi_N 

%% ���ļ����ڵ�*.pre�ļ����в�����Ĭ��ֻ��һ��.pre�ļ�

    editfekoname=A.name;
    %%���µĲ��м�������Լ����к��������ã�����Ĭ�ϵ��ǵ��ñ���ȫ�������������ϣ��2���ˣ�����np allΪ np 2
%     d = strcat([feko_path,32,editfekoname]) ;% ��Ҫ��all ��2����4������
   d = strcat([feko_path,32,editfekoname,' -np 2 --parallel-authenticate localonly']);  % -np 4  d = strcat([feko_path,32,filename,' -np all']) 
dos(d); 

%% �ڶ����֣�Ѱ�ҡ���ȡ.out�ļ�
lengthE = length(editfekoname);
lengthP = length('.pre');
tmp_name = editfekoname(1:lengthE-lengthP);
outname=strcat(tmp_name,'.out');

B = dir(fullfile('',outname));
filename_out=strcat(filepath,'\',B.name);

%��ȡ.out�ļ��ּ�
[RCS_complex_theta,RCS_complex_phi,RCS_value]=FekoReadFarfieldRCS_MonoStatic_OnlyForRCS(filename_out,OriHz);
%PlotRangeSequence(Phi_deg,Freq_Hz,RCS_complex_theta,1024);%�ú���ͬ��������(�ȼ���)


RCS_theta(:,ii)=RCS_complex_theta;
RCS_phi(:,ii)=RCS_complex_phi;
RCS_V(:,ii)=RCS_value;


%% �������֣��޸�.pre�ļ���������
% fid = fopen(editfekoname,'r');
CharPhi='A0: 0 :  : 1 : 1 : 1 : 1 : 0 : ';

Phi_new = OriPhi + Phi_Step*ii; % ����֤�ĽǶ�

Phi_new_char = num2str(Phi_new); % ת�����ַ���
new_contents_phi = strcat("A0: 0 :  : 1 : 1 : 1 : 1 : 0 : ",num2str(OriTheta)," : ",Phi_new_char," : 0 : 0 : ",num2str(Phi_Step)," : 0   ** PlaneWaveSource1");
new_contents_phi = new_contents_phi + newline;


text_modify(editfekoname,editfekoname,new_contents_phi,CharPhi)%��ȡ�ļ����е��ļ������µ�����new_contents�滻 ָ���еľɵ�����old_contents

str_bar = ['��ʼphi��Ϊ',num2str(OriPhi),' ������...',num2str(100*ii/Phi_N),'%'];
waitbar(ii/Phi_N,bar,str_bar);
ii = ii + 1
end % while ii ����
t = toc % ���Գ�������ʱ���β�����ڲ鿴��������ʱ�䡣

%% ����ز�����
save RCS_theta.mat RCS_theta
save RCS_phi.mat RCS_phi
save RCS_V.mat RCS_V

%% ���Ĳ��֣� ����
%% ������ز��� ��ʼ��(Ϊ���Ĳ��ֵĳ�����׼������ǰ�������޹أ�
subplot_n=3; %��ͼʱ��Ĭ�ϴճ�3*3��subplot,���Ը����������
a_num=OriHz; %��ͼʱ����λ��ĳ��ȣ�Ĭ�Ϻ;�����Ƶ��������һ��
a_gap=floor(Phi_N/(subplot_n*subplot_n)); %Ĭ�ϴճ�3*3��subplot,���Ը����������
k_start=1; %��ͼʱ����λ������
k_end=Phi_N; %�����һ����ͼ����λ��
% k_end=Phi_N;
ColorRange=-50; %caxis�����ĵ�һ������ֵ

load RCS_theta.mat
ImageOutput(subplot_n,a_num,a_gap,k_start,k_end,...
   RCS_theta,BW,C_freq,c,Phi_Step,OriPhi,ColorRange,'ETHETA')
% TorP='ETHETA' or 'EPHI', RCS_TorP �е�TorP������ TorP ��Ӧ
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
% title(['ETHETA  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'��']);
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
% title(['EPHI  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'��']);
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
% set( gca, 'XTick', [], 'YTick', [] );%����ʾ������XY
% colorbar;
% title('��λ��-6�㵽6��');
% xlabel('��λ��');ylabel('������');
% PlotImage(Im,-y,x,[],[],50);