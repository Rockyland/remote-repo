load RCS_theta.mat

a_num=256;
a_gap=300;
k_start=1;
k_end=4600;

figure;
for k=k_start:a_gap:k_end
S=hamming(size(RCS_theta,1))*hamming(size(RCS_theta(:,k:k+a_num),2))'.*RCS_theta(:,k:k+a_num);

Im=fty(ftx(S));
G=mag2db(abs(Im)/max(max(abs(Im))));
[N,M]=size(G);
B=20E9;
c=3e8;
f0=216e9;
lambda=c/f0;
Omiga=Phi_Step*a_num;
theta=deg2rad(Omiga);
range_a=lambda/2/theta;
range_r=c/2/B;
subplot(4,4,floor(k/a_gap)+1)
imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
% imagesc(G)
colorbar
colormap hot
caxis([-10,0])
xlabel('Azimuth(m)');
ylabel('Range(m)');
% xlim([-0.3,0.3]);
% ylim([-0.3,0.3]);
title(['ETHETA  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
% ftresize(18)
end


load RCS_phi.mat


figure;
for k=k_start:a_gap:k_end
S=hamming(size(RCS_phi,1))*hamming(size(RCS_phi(:,k:k+a_num),2))'.*RCS_phi(:,k:k+a_num);

Im=fty(ftx(S));
G=mag2db(abs(Im)/max(max(abs(Im))));
[N,M]=size(G);
B=20E9;
c=3e8;
f0=216e9;
lambda=c/f0;
Omiga=Phi_Step*a_num;
theta=deg2rad(Omiga);
range_a=lambda/2/theta;
range_r=c/2/B;
subplot(4,4,floor(k/a_gap)+1)
imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
% imagesc(G)
colorbar
colormap hot
caxis([-10,0])
xlabel('Azimuth(m)');
ylabel('Range(m)');
% xlim([-0.3,0.3]);
% ylim([-0.3,0.3]);
title(['EPHI  ','phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
% ftresize(18)
end