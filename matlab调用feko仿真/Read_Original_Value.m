% 该函数用来读取 .pre文件的3个初始值：频点数、theta初始角、phi初始角
function [OriHz,S_freq,E_freq,OriTheta,OriPhi] =  Read_Original_Value(filename)
fid = fopen(filename,'r');

Char1='** Set frequency';
Char2='** Sources';
%   THETA    PHI      magn.    phase     magn.    phase        in m*m                     axial r. angle   direction
while ~feof(fid)
    tline = fgets(fid);%必须用fgets不能用fgetl
    logic1 = strncmp(Char1,tline,length(Char1));
    logic2 = strncmp(Char2,tline,length(Char2));
    if logic1
        tmp1 =fscanf(fid,'%*s%s',[1,1]); % tmp1是频率点数
        tmp4 =fscanf(fid,'%*s%*s%*s%*s%*s%*s%s',[1,1]); % tmp4是电磁波起点频率
        tmp5 =fscanf(fid,'%*s%*s%s',[1,1]); % tmp5是电磁波终点频率
        OriHz=str2double(tmp1);
        S_freq=str2double(tmp4);
        E_freq=str2double(tmp5);
    end
    if logic2
        tmp2 =fscanf(fid,'%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s %s',[1,1]);%tmp2是theta初始角
        tmp3 =fscanf(fid,'%*s %s',[1,1]); % tmp3是phi初始角
        OriTheta=str2double(tmp2);
        OriPhi=str2double(tmp3);
    end
end
fclose(fid);
end

