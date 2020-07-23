%interp_image
% Bitmap image interpolation function.
%
% LAST UPDATED by Andy French 13th March 2011

function II = interp_image( I , II_width , II_height , fit_by )

%Message to command window
disp('.... Reducing image size via interpolation')

%Get dimensions of I
[heightI,widthI,num_cols]  = size(I);

%Define index vectors for height and width of 'initial_image'
w = 1:widthI;
h = 1:heightI;

%Define index vectors for height and width of output image
if strcmp(fit_by,'height')
    W = linspace(1,widthI,II_height);
    H = linspace(1,heightI,II_height );
elseif strcmp(fit_by,'width')
    W = linspace(1,widthI,II_width);
    H = linspace(1,heightI,II_width );
else
    W = linspace(1,widthI,II_width);
    H = linspace(1,heightI,II_height );
end

%Meshgrid these
[h,w] = meshgrid(h,w);
[H,W] = meshgrid(H,W);

%If I is a grayscale image, then replicate for R,G,B components
if num_cols ==1
    II = interp2(h,w,double(I.'),H,W).';
else
    %Create R,G,B matricies from image matrix I.
    R = ones(heightI,widthI);
    G = ones(heightI,widthI);
    B = ones(heightI,widthI);
    R(:,:) = I(:,:,1);
    G(:,:) = I(:,:,2);
    B(:,:) = I(:,:,3);
    
    %Interpolate these colour matricies
    R = interp2(h,w,R.',H,W);
    G = interp2(h,w,G.',H,W);
    B = interp2(h,w,B.',H,W);
    
    %Initialise output array
    dim = size(R.');
    II = zeros(dim(1),dim(2),3);
    II(:,:,1) = R.';
    II(:,:,2) = G.';
    II(:,:,3) = B.';
end

%End of code