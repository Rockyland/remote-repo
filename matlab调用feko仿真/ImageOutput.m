
% 对用matlab_union_main01生成的回波进行成像
function ImageOutput(subplot_n,a_num,a_gap,k_start,k_end,...
   RCS_TorP,BW,C_freq,c,Phi_Step,OriPhi,ColorRange,TorP)

    [~,n]=size(RCS_TorP);
    B=BW;
    f0=C_freq;
    lambda=c/f0;
    Omiga=Phi_Step*a_num;
    theta=deg2rad(Omiga);
    range_a=lambda/2/theta;
    range_r=c/2/B;
    kn=range_r;
    
for k=k_start:a_gap:k_end
    S=hamming(size(RCS_TorP,1))*hamming(size(RCS_TorP(:,k:k+a_num),2))'.*RCS_TorP(:,k:k+a_num);
    
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
    S=hamming(size(RCS_TorP,1))*hamming(size(RCS_TorP(:,k:k+a_num),2))'.*RCS_TorP(:,k:k+a_num);
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%小目标，暂时不需要。
%     %Keystone变换之后，再成像
%     figure
%     for k=k_start:a_gap:k_end
% 
%         if k+a_num<n
%             S=hamming(size(RCS_TorP,1))*hamming(size(RCS_TorP(:,k:k+a_num),2))'.*RCS_TorP(:,k:k+a_num);
%         else
%             S=hamming(size(RCS_TorP,1))*hamming(size(RCS_TorP(:,k:n),2))'.*RCS_TorP(:,k:n);
%         end
% 
%         S_Keystone=Keystone(S,B,f0);
% 
%          % Keystone之后，做补偿，再成像
%         hrrp=ftx(S_Keystone);
% %         [N,M]=size(phi_rcs);
%         w=deg2rad(Omiga)^2;
%         y0=0;
%         PRF=M;
%         alpha=2*pi*f0/(c*PRF^2);
%         phase=exp(1i*alpha*([-N/2:N/2-1].'*kn+y0)*w*([-M/2:M/2-1].^2));
%         hrrp_phasecomp=hrrp.*phase;
% %         hrrp_phasecomp=ftx(hamming(1024)*hamming(1024)'.*iftx(hrrp_phasecomp),4096);
% %         IM=abs(fty(hrrp_phasecomp,4096));
%         IM=abs(fty(hrrp_phasecomp));
% 
%         IM=mag2db(IM/max(max(IM)));
%         subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
%         imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,IM)
%         % imagesc(G)
%         colorbar
%         colormap hot
%         caxis([ColorRange,0])
%         xlabel('Azimuth(m)');
%         ylabel('Range(m)');
%         title([TorP,'  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°',' 相位补偿后']);
%     end
    
   
    
    
    
    
end