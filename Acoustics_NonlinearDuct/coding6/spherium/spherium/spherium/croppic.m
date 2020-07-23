%croppic
% Automatic cropping of a bitmap image. This function determines the
% largest rectangular bounding box around pixels in the image which are NOT
% the background colour. The latter is specified as a [R,G,B] vector with
% values 0.... 255 for red, green and blue.

function croppic( filename, background_colour )

%Load image
I = imread(filename);

%Determine maximum and minimum indices of pixels which are NOT the
%paper colour
[row_ind,col_ind] =...
    find( (( I(:,:,1)==background_colour(1) ) .*...
    ( I(:,:,2)==background_colour(2) ) .*...
    ( I(:,:,3)==background_colour(3) ))==0 );
min_row = min(row_ind);
max_row = max(row_ind);
min_col = min(col_ind);
max_col = max(col_ind);

%Crop image
I = I( min_row : max_row, min_col:max_col, : );

%Save image
imwrite( uint8(I) , filename );

%End of code