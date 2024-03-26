
% 对用matlab_union_main01生成的回波进行成像
function ImageOutput(subplot_n,a_num,a_gap,k_start,k_end,...
   RCS_TorP,BW,C_freq,c,Phi_Step,OriPhi,ColorRange,TorP)

    [~,n]=size(RCS_TorP);

    figure;%直接成像，不做补偿
    for k=k_start:a_gap:k_end

        if k+a_num<n
            S=hamming(size(RCS_TorP,1))*hamming(size(RCS_TorP(:,k:k+a_num),2))'.*RCS_TorP(:,k:k+a_num);
        else
            S=hamming(size(RCS_TorP,1))*hamming(size(RCS_TorP(:,k:n),2))'.*RCS_TorP(:,k:n);
        end

        Im=fty(ftx(S));
        G=mag2db(abs(Im)/max(max(abs(Im))));
        [N,M]=size(G);
        B=BW;
        f0=C_freq;
        lambda=c/f0;
        Omiga=Phi_Step*a_num;
        theta=deg2rad(Omiga);
        range_a=lambda/2/theta;
        range_r=c/2/B;
        kn=range_r;
        PRF=M;
        alpha=2*pi*f0/(c*PRF^2);
        subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
        imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
        % imagesc(G)
        colorbar
        colormap hot
        caxis([ColorRange,0])
        xlabel('Azimuth(m)');
        ylabel('Range(m)');
        title([TorP,'  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
    end
    
    %Keystone变换之后，再成像
    figure
    for k=k_start:a_gap:k_end
         
        if k+a_num<n
            S=hamming(size(RCS_TorP,1))*hamming(size(RCS_TorP(:,k:k+a_num),2))'.*RCS_TorP(:,k:k+a_num);
        else
            S=hamming(size(RCS_TorP,1))*hamming(size(RCS_TorP(:,k:n),2))'.*RCS_TorP(:,k:n);
        end

        S_Keystone=Keystone(S,B,f0);
%         IMRDK=abs(fty(ftx(S_Keystone)));
%         IMRDK=mag2db(IMRDK/max(max(IMRDK)));
%         subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
%         imagesc(range_a*linspace(-0.5,0.5,M)*M,N*linspace(-0.5,0.5,N)*range_r,IMRDK)
%         xlabel('Azimuth (m)');
%         ylabel('Range (m)')
    %     ftresize(18)
%         colorbar
%         colormap hot
%         caxis([ColorRange,0])
        
         % Keystone之后，做补偿，再成像
        hrrp=ftx(S_Keystone);
%         [N,M]=size(phi_rcs);
        w=deg2rad(Omiga)^2;
        y0=0;
        phase=exp(1i*alpha*([-N/2:N/2-1].'*kn+y0)*w*([-M/2:M/2-1].^2));
        hrrp_phasecomp=hrrp.*phase;
%         hrrp_phasecomp=ftx(hamming(1024)*hamming(1024)'.*iftx(hrrp_phasecomp),4096);
%         IM=abs(fty(hrrp_phasecomp,4096));
        IM=abs(fty(hrrp_phasecomp));
     
        IM=mag2db(IM/max(max(IM)));
        subplot(subplot_n,subplot_n,floor(k/a_gap)+1)
        imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,IM)
        % imagesc(G)
        colorbar
        colormap hot
        caxis([ColorRange,0])
        xlabel('Azimuth(m)');
        ylabel('Range(m)');
        title([TorP,'  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°',' 相位补偿后']);
    end
    
   
    
    
    
    
end