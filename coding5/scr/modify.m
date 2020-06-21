function modify(DATA,IGV_angle)

[fname,location]=uigetfile({'*.geomTurbo';'*.txt';'*.dat';'*.*'},'r');
%copyfile('DTS-02.geomTurbo', 'modify6.geomTurbo');  %先拷贝一份作为备份,modify作为修改文件
    fid=fopen(fname,'r');
    fd=fopen([['New_angle_',num2str(IGV_angle),'_'] ,fname],'w');

%words=input('Give the string you want to delete:\n','s');
%words='47';%标志物

     mark=0;%标记mark=1（IGV_suction）；mark=2（IGV_pressure）；mark=3（rotor_suction）；mark=4（rotor_pressure）；mark=5（stator_suction）；mark=6（stator_pressure）；
    
    for j=1:10000  %第j行
        if ~feof(fid)  %只要文档还开着
           
            line=fgetl(fid);%就一行一行复制下去
             if (~isempty(strfind(line,'suction')))|(~isempty(strfind(line,'pressure'))) %当复制到标志点时  |strfind(line,'pressure')
                         mark=mark+1;
                fprintf(fd,'%s\r\n',line);  %但还在复制该信息 
                line=fgetl(fid);
                fprintf(fd,'%s\r\n',line); 
                section = fscanf(fid, '%d',1);  %获取某些参量信息
                fprintf(fd,'%s',num2str(section));  %但还在复制该信息    %%整个函数结构可以大幅度提升
                continue;
             elseif isempty(strfind(line,['#','   section ']))    % else if isempty(strfind(line,['#','   section ',num2str(k)]))
                 temp=0;  %标记状态
                 fprintf(fd,'%s\r\n',line);
             else if ~isempty(strfind(line,['#','   section ']))
                        temp=1;
                        fprintf(fd,'%s\r\n',line);
                        k=str2num(line(isstrprop(line,'digit')));  %提取section number'
                        line=fgetl(fid);
                        fprintf(fd,'%s\r\n',line);
                        numberOfPoints(k) =fscanf(fid, '%d',1);
                        fprintf(fd,'%s\r',num2str(numberOfPoints(k)));   
                        line=fgetl(fid);
                        fprintf(fd,'%s\r\n',line);
                        
                        if mark==1
                             matrix=DATA.IGV_suction_rotate{1,k};  %替换新数据
                        elseif mark==2
                             matrix=DATA.IGV_pressure_rotate{1,k};  %替换新数据
                        elseif mark==3
                             matrix=DATA.rotor_suction{1,k};  %替换新数据
                       elseif mark==4
                             matrix=DATA.rotor_pressure{1,k};  %替换新数据
                       elseif mark==5
                             matrix=DATA.stator_suction{1,k};  %替换新数据
                       elseif mark==6
                             matrix=DATA.stator_pressure{1,k};  %替换新数据
                        end
                        
                        [m,n]=size(matrix);
                            if m~=numberOfPoints(k)
                                disp('新数据有误')
                                break;
                             end
                        for kk=1:numberOfPoints(k)
                            temp=2;
                            line=fgetl(fid);
                            %fprintf(fd,'%s\r\n',line);  %替换！
                            
                                fprintf(fd,'%f %f %f\n',matrix(kk,:)); 
                            
                        end    
                                
                               
                            
                        end
                            
                      
               %%函数完整性得以保证，但是还无法想到有效加入数据的结构
                           
                         end
                     
                   
                 
             

             
        else
            break;
        end
    end
     fclose(fd);
    if(temp>10)
        delete(['G' cell2mat(fname(i))]);
    end
    fclose(fid);
    disp('导出成功！！')
   
end
