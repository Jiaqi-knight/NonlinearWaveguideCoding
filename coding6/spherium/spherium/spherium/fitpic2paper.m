%fitpic2paper
% Fits a bitmap image to a defined paper size aspect ratio.

function fitpic2paper( filename, aspect, paper_color )

%Load image
I = imread(filename);

%Find updated pixel dimensions of cropped image
dim = size(I);
width_I = dim(2);
height_I = dim(1);

%Determine aspect ratio
aspect_I = width_I/height_I;

%Depending on aspect ratios of initial_image and desired image, choose
%image forming method
if aspect<aspect_I
    
    %Output image is 'taller' than I. Need to pad in height, fit by width
    num_blanks = round( 0.5 * abs(height_I*aspect_I/aspect - height_I) );
    I = addblanks( I , num_blanks , 'height' , paper_color );
    
else
    
    %Output image is 'wider' than I. Need to pad in width, fit by height
    num_blanks = round( 0.5 * abs(width_I*aspect/aspect_I - width_I) );
    I = addblanks( I , num_blanks , 'width' , paper_color );
end

%Save image
imwrite( uint8(I) , filename );

%%

%Add blanks function
function II = addblanks( I , num_blanks , height_or_width , blank_RGB )

%Get dimensions of I
[heightI,widthI,num_cols]  = size(I);

%Create R,G,B matricies from image matrix I.
R = ones(heightI,widthI);
G = ones(heightI,widthI);
B = ones(heightI,widthI);
R(:,:) = I(:,:,1);
G(:,:) = I(:,:,2);
B(:,:) = I(:,:,3);

%Construct new RGB array II with required number of 'blanks' added
if strcmp(height_or_width,'width')
    II = ones( heightI,widthI+2*num_blanks,3);
    II(:,:,1) = [ blank_RGB(1)*ones(heightI,num_blanks) ,...
        R , blank_RGB(1)*ones(heightI,num_blanks) ];
    II(:,:,2) = [ blank_RGB(2)*ones(heightI,num_blanks),...
        G , blank_RGB(2)*ones(heightI,num_blanks) ];
    II(:,:,3) = [ blank_RGB(3)*ones(heightI,num_blanks) ,...
        B , blank_RGB(3)*ones(heightI,num_blanks) ];
else
    II = ones( heightI+2*num_blanks,widthI,3);
    II(:,:,1) = [ blank_RGB(1)*ones(num_blanks,widthI) ;...
        R ; blank_RGB(1)*ones(num_blanks,widthI) ];
    II(:,:,2) = [ blank_RGB(2)*ones(num_blanks,widthI) ;...
        G ; blank_RGB(2)*ones(num_blanks,widthI) ];
    II(:,:,3) = [ blank_RGB(3)*ones(num_blanks,widthI) ;...
        B ; blank_RGB(3)*ones(num_blanks,widthI) ];
end

%End of code