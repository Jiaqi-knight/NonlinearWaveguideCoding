function Data=ReadgeomTurbo(fname)
      ZR_mark=0;%标志hub，shround线
      Z_mark=0;%标志IGV_leading=1;IGV_trailing=2;rotor_leading=3;rotor_trailing=4;stator_leading=5;stator_trailing=6;
      mark=0;%标志IGV_suction-1;IGV_pressure-2;rotor_suction-3;rotor_pressure-4;stator_suction-5;stator_pressure-6;
      fid=fopen(fname,'r');     
      
    for j=1:10000  %第j行
        if ~feof(fid)  %只要文档还开着
            
            line=fgetl(fid);%就一行一行读下去
               if (~isempty(strfind(line,'NI_BEGIN zrcurve')))  %读取hub和shround线
                ZR_mark=ZR_mark+1;              
                line=fgetl(fid);
                ZR_numberOfPoints(ZR_mark) = fscanf(fid, '%d',1);  %获取某些参量信息
                for n = 1:ZR_numberOfPoints(ZR_mark)
                         fscanf(fid, '%c',1);
                         z=fscanf(fid, '%f',1);
                         r=fscanf(fid, '%f',1);   
                         fscanf(fid, '%c',1);    
                         if ZR_mark==1
                         Z_hub(n,1)= [z];
                         R_hub(n,1)= [r];
                         else
                         Z_shround(n,1)= [z];
                         R_shround(n,1)= [r];
                         end
                             
                 end
               
                                    switch ZR_mark
                            case 1
                                Data.hub=[Z_hub,R_hub];
                                clear Z_hub
                                clear R_hub
                            case 2
                                Data.shround=[Z_shround,R_shround];
                                clear Z_shround
                                clear Z_shround
                                    end
                                     
                 end                    
%% 读取前缘和尾缘数据

             
                if ((~isempty(strfind(line,'NI_BEGIN nisolid_angle_at_leading_edge')))|(~isempty(strfind(line,'NI_BEGIN nisolid_angle_at_trailing_edge'))))  %读取hub和shround线
                  Z_mark=Z_mark+1;              
                  %line=fgetl(fid);
                  Z_numberOfPoints(Z_mark) = fscanf(fid, '%d',1);  %获取某些参量信息
                  for n = 1:Z_numberOfPoints(Z_mark)
                         fscanf(fid, '%c',1);
                         x=fscanf(fid, '%f',1);
                         y=fscanf(fid, '%f',1); 
                         z=fscanf(fid, '%f',1); 
                         fscanf(fid, '%c',1);    
                         X1(n,:)= [x,y,z];
                         end
                             
                
            
                                    switch Z_mark
                                                                    case 1
                                Data.IGV_leading=[X1];
                            case 2
                                Data.IGV_trailing=[X1];
                            case 3
                                Data.rotor_leading=[X1]; 
                            case 4
                                Data.rotor_trailing=[X1];
                            case 5
                                Data.stator_leading=[X1];
                            case 6
                                Data.stator_trailing=[X1];
                            end
                                    clear X1
                end 
                
%% 读取 叶型数据    
             if (~isempty(strfind(line,'suction')))|(~isempty(strfind(line,'pressure'))) %当复制到标志点时  |strfind(line,'pressure')
                mark=mark+1;              
                line=fgetl(fid);
                section = fscanf(fid, '%d',1);  %获取某些参量信息
                continue;
             elseif isempty(strfind(line,['#','   section ']))    % else if isempty(strfind(line,['#','   section ',num2str(k)]))
                 temp=0;  %标记状态
             elseif ~isempty(strfind(line,['#','   section ']))
                        temp=1;
                        k=str2num(line(isstrprop(line,'digit')));  %提取文本中的数字
                        line=fgetl(fid);
                        numberOfPoints(mark,k) =fscanf(fid, '%d',1);
                        
                       for n = 1:numberOfPoints(mark,k)
                         fscanf(fid, '%c',1);
                         x=fscanf(fid, '%f',1);
                         y=fscanf(fid, '%f',1);
                         z=fscanf(fid, '%f',1);
                         fscanf(fid, '%c',1);
                        
                         X(n,:)= [x,y,z]; 
                       end
                       
                        switch mark
                            case 1
                                Data.IGV_suction{k}=[X];    
                            case 2
                                Data.IGV_pressure{k}=[X];
                            case 3
                                Data.rotor_suction{k}=[X];
                            case 4
                                Data.rotor_pressure{k}=[X];
                            case 5
                                Data.stator_suction{k}=[X];
                            case 6
                                Data.stator_pressure{k}=[X];
                        end
                        clear X  %不能加封号
                        
                 end
             end
       
    end
      Data.numberOfPoints=numberOfPoints; 
    end
