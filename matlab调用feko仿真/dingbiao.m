for k=100:300:4800
S=hamming(size(RCS_phi,1))*hamming(size(RCS_phi(:,k:k+256),2))'.*RCS_phi(:,k:k+256);

Im=fty(ftx(S));
G=mag2db(abs(Im)/max(max(abs(Im))));
[N,M]=size(G);
B=20E9;
c=3e8;
f0=216e9;
lambda=c/f0;
Omiga=Phi_Step*256;
theta=deg2rad(Omiga);
range_a=lambda/2/theta;
range_r=c/2/B;
% figure;
subplot(4,4,floor(k/300)+1)
imagesc(range_a*linspace(-0.5,0.5,M)*M ,N*linspace(-0.5,0.5,N)*range_r,G)
% imagesc(G)
colorbar
colormap hot
caxis([-50,0])
xlabel('Azimuth(m)');
ylabel('Range(m)');
xlim([-0.3,0.3]);
ylim([-0.3,0.3]);
title(['phi = ',num2str(round(OriPhi+Phi_Step*k,1)),'°']);
% ftresize(18)
end
