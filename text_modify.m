function text_modify(filename_ini,filename_new,new_contents_phi,CharPhi)%读取文件夹中的文件，用新的内容new_contents替换line行的旧的内容old_contents

%% 读取文件并修改
 fileID = fopen(filename_ini,'r+');                    %以可读写的方式打开待修改的文件
 i=0;
 while ~feof(fileID)
        tline=fgets(fileID);
        logic = strncmp(CharPhi,tline,length(CharPhi));%逐行读取原始文件
        i=i+1;
        newrow{i} = tline;                               %创建新对象接受原始文件每行数据
        if logic                                        %判断是否到达待修改的行
            newrow{i} = strrep(tline,tline,new_contents_phi);%替换函数，将待修改的行全部替换为新内容
        end
 end
    fclose(fileID);                                       %关闭文件
    
%% 输出文本
    fileID = fopen(filename_new,'w+');                    %以可读写的方式打开输出文件，文件若存在则清空文件内容从文件头部开始写，若不存在则根据文件名创建新文件并只写打开
    for k=1:i
        fprintf(fileID,'%s',newrow{k});              %将newline内的内容逐行写出
    end
    fclose(fileID);                                       %关闭文件
end