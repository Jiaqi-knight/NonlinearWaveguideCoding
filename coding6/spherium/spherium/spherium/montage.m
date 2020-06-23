%montage
% Function that takes a set of images of common filenames and arranges them
% in a 'near square' format, the saves this montaged image in PNG format.
%
% LAST UPDATED by Andy French 30-November-2007

function montage(file_directory,common_filename,delete_files)

%Start clock
t0=clock;

%Get filenames in file_directory which include common_filename
d=dir(file_directory);
num_files=length(d);
disp(['    Montaging files in directory ',file_directory,...
    ' with common filename ',common_filename])

%Step through dir and remove any file names with '.' or '..'
i=0;
filenames=cell(1,1);
for n=1:num_files
    if strcmp(d(n).name,'.') | strcmp(d(n).name,'..')
    else
        i=i+1;
        filenames{i} = d(n).name;
    end
end

%Step through files in filenames and only load those with a valid image
%file extension and contain filename "common_filename"
num_files=length(filenames);
if ~isempty(filenames{1})
    i=0;
    ok_files=cell(1,1);
    for n=1:num_files

        %Obtain file extension
        index_of_dot = strfind(filenames{n},'.');
        fname = filenames{n};
        extension = fname(index_of_dot+1:end);

        %Check that filenames contain common_filename and are recognised image
        %formats
        if ( ~isempty( strmatch( common_filename, fname  ) ) ) &...
                ( ~isempty( imformats( extension ) ) )
            i=i+1;
            ok_files{i}=fname;
        end
    end

    %Only proceed if there are some OK files!
    if ~isempty( ok_files{1} )

        %Determine number of OK files and compute rows and columns of the
        %montage
        num_files = length(ok_files);
        rows = ceil( sqrt(num_files) );
        cols = rows - 1;
        if num_files>rows*cols
            cols=cols+1;
        end

        %Sequentially read files into cell array I and form arrays of X and Y
        %pixel size
        I = cell(1,num_files);
        for n=1:num_files
            disp(['    Reading file ',file_directory,'\',ok_files{n},' .....'])
            I{n} = imread( [ file_directory,'\',ok_files{n} ] );
            dim = size(I{n}) ;
            Xpixels(n) = dim(1);
            Ypixels(n) = dim(2);
        end
        max_Xpixels = max(Xpixels);
        max_Ypixels = max(Ypixels);

        %Sequentially add data to matrix M, padding with zeros (white) if image
        %size is less than largest in X and Y directions
        M=[];
        file_num=1;
        for c=1:cols
            R=[];
            for r=1:rows
                if file_num>0
                    
                    disp(['    Adding file ',file_directory,'\',ok_files{file_num},' to montage .....'])

                    %Pad X and Y pixel ones. (White) Attempt to center images in
                    %space of largest
                    start_X_pixels = round( (max_Xpixels-Xpixels(file_num) )/2 );
                    end_X_pixels = max_Xpixels-Xpixels(file_num) - start_X_pixels ;
                    start_Y_pixels = round( (max_Ypixels-Ypixels(file_num) )/2 );
                    end_Y_pixels = max_Ypixels-Ypixels(file_num) - start_Y_pixels ;
                    
                    I{file_num} = [ones(start_X_pixels,Ypixels(file_num),3);...
                        I{file_num};...
                        ones(end_X_pixels,Ypixels(file_num),3)];
                    I{file_num} = [ones(max_Xpixels,start_Y_pixels,3),...
                        I{file_num},...
                        ones(max_Xpixels,end_Y_pixels,3)];

                    %Add to row image and then remove data from I{file_num}
                    %to save memory
                    R = [R,I{file_num}];
                    I{file_num} = [];

                    %Update file number
                    file_num = file_num+1;
                    if file_num>num_files
                        file_num=0;
                    end
                end
            end
            
            %Check columns of R match M, if not pad with ones. Note
            %dim_M(2) will always be >= than dim_R(2)
            dim_M = size(M);
            dim_R = size(R);
            if dim_R(2)~=dim_M(2)
                R = [R,ones(dim_R(1),dim_M(2)-dim_R(2),3)];
            end
            M=[M ; R];
        end

        %Save montage
        disp('... Saving montage')
        imwrite(M,[file_directory,'\','Montage ',common_filename,'.png'],'png')
        
        %Delete composite files
        if delete_files==1
            for n=1:num_files
                delete([file_directory,'\',ok_files{n}]);
            end
        end
    end
end

%Stop clock and display ending message with elapsed time
disp(['... Montage complete in ',num2str( etime(clock,t0) ),' seconds.'])

%End of code

