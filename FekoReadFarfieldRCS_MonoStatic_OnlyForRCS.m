function [RCS_complex_theta,RCS_complex_phi,RCS_value]=FekoReadFarfieldRCS_MonoStatic_OnlyForRCS(filename,OriHz)
% 

RCS_complex_theta=zeros(OriHz,1);
RCS_complex_phi=zeros(OriHz,1);
RCS_value=zeros(OriHz,1);
k = 1;

fid = fopen(filename,'r');

Char1='   THETA    PHI      magn.    phase     magn.    phase        in m*m                     axial r. angle   direction';
%   THETA    PHI      magn.    phase     magn.    phase        in m*m                     axial r. angle   direction
while ~feof(fid)
    tline = fgets(fid);%必须用fgets不能用fgetl
    logic = strncmp(Char1,tline,length(Char1));
    if logic
        tmp1 =fscanf(fid,'%*g%*g %g%g %g%g',[1,4]);
%         RCS_complex_theta=[RCS_complex_theta;2*sqrt(pi)*tmp1(1)*exp(1j*deg2rad(tmp1(2)))];
        RCS_complex_theta(k)=2*sqrt(pi)*tmp1(1)*exp(1j*deg2rad(tmp1(2)));
%         RCS_complex_phi=[RCS_complex_phi;2*sqrt(pi)*tmp1(3)*exp(1j*deg2rad(tmp1(4)))];
        RCS_complex_phi(k)=2*sqrt(pi)*tmp1(3)*exp(1j*deg2rad(tmp1(4)));
        RCS_value(k)=tmp1(5);
        k = k+1; 
    end
end

if 1==0  %相当于注释掉，因为在Y9000P 上总是报错
    RCS_complex_theta=reshape(RCS_complex_theta,1,1,1200);
    % RCS_complex_theta=reshape(RCS_complex_theta,length(Theta_deg),501,1200);
    RCS_complex_theta=permute(RCS_complex_theta,[3,2,1]);
    RCS_complex_theta=squeeze(RCS_complex_theta);

    RCS_complex_phi=reshape(RCS_complex_phi,1,1,1200);
    % RCS_complex_phi=reshape(RCS_complex_phi,length(Theta_deg),501,1200);
    RCS_complex_phi=permute(RCS_complex_phi,[3,2,1]);
    RCS_complex_phi=squeeze(RCS_complex_phi);
end
fclose(fid);
