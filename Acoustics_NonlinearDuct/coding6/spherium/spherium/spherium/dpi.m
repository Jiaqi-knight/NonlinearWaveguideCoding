%dpi
% Function which returns height and width of an image of x dots per inch
% with defined width and height in cm.
%
% LAST UPDATED by Andy French 29-Dec-2007
%
% Syntax: [width,height] = dpi( x , width_cm , height_cm)
% OR dpi( x , paper)
%
% x is DPI of image
% a is width / height aspect ratio. i.e. 6/4 for a standard photo.
%
% paper is a string i.e. 'A4'
%
% width, height are dimensions in pixels of desired image.
%
% Example: [width,height] = dpi( 300 , 'A4' )

function [width,height] = dpi( varargin)

%Assign (variable) inputs
if nargin==2
    x = varargin{1};
    paper = varargin{2};
    [width_cm,height_cm] = papersize(paper);
elseif nargin==3
    x = varargin{1};
    width_cm = varargin{2};
    height_cm = varargin{3};
else
    width=[];
    height=[];
end

%Number of centimeters per inch
cm_per_inch = 1/0.3937;

%Pixels per centimeter
pixels_per_cm = x / cm_per_inch ;

%Let x be the required number of pixels per inch. Convert width_cm and
%height_cm to numbers of pixels
width = round( width_cm * pixels_per_cm );
height = round( height_cm * pixels_per_cm );

%End of code