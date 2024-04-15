clear;clc
lim=30;
[suf, ~, ~, c] = stlread("tank.stl");
sufPoints(:,3)=suf.Points(:,3);
sufPoints(:,2)=suf.Points(:,2);
sufPoints(:,1)=suf.Points(:,1);
figure
p = patch('Faces',suf.ConnectivityList,'Vertices',sufPoints(:,:));
p.FaceVertexCData=abs(sufPoints(:,1));
shading interp
daspect([1 1 1]);
camlight left;lighting gouraud
axis([-lim,lim,-lim,lim,-lim/2,lim/2]);
xl=xlabel('Range (cm)');yl=ylabel('Azimuth (cm)');zlabel('Elevation (cm)');
view(90,90)