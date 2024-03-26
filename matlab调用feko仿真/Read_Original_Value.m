% �ú���������ȡ .pre�ļ���3����ʼֵ��Ƶ������theta��ʼ�ǡ�phi��ʼ��
function [OriHz,S_freq,E_freq,OriTheta,OriPhi] =  Read_Original_Value(filename)
fid = fopen(filename,'r');

Char1='** Set frequency';
Char2='** Sources';
%   THETA    PHI      magn.    phase     magn.    phase        in m*m                     axial r. angle   direction
while ~feof(fid)
    tline = fgets(fid);%������fgets������fgetl
    logic1 = strncmp(Char1,tline,length(Char1));
    logic2 = strncmp(Char2,tline,length(Char2));
    if logic1
        tmp1 =fscanf(fid,'%*s%s',[1,1]); % tmp1��Ƶ�ʵ���
        tmp4 =fscanf(fid,'%*s%*s%*s%*s%*s%*s%s',[1,1]); % tmp4�ǵ�Ų����Ƶ��
        tmp5 =fscanf(fid,'%*s%*s%s',[1,1]); % tmp5�ǵ�Ų��յ�Ƶ��
        OriHz=str2double(tmp1);
        S_freq=str2double(tmp4);
        E_freq=str2double(tmp5);
    end
    if logic2
        tmp2 =fscanf(fid,'%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s %s',[1,1]);%tmp2��theta��ʼ��
        tmp3 =fscanf(fid,'%*s %s',[1,1]); % tmp3��phi��ʼ��
        OriTheta=str2double(tmp2);
        OriPhi=str2double(tmp3);
    end
end
fclose(fid);
end

