%papersize
% Function which gives standard paper sizes in cm.
%
% References: 
% http://www.cl.cam.ac.uk/~mgk25/iso-paper.html
% http://en.wikipedia.org/wiki/Photo_print_sizes

function [width_cm,height_cm] = papersize(paper)

if strcmp( paper,'A6 landscape')==1
    [height_cm,width_cm] = ABCsize('A6');

elseif strcmp( paper,'A6 portrait')==1
    [width_cm,height_cm] = ABCsize('A6');

elseif strcmp( paper,'A5 landscape')==1
    [height_cm,width_cm] = ABCsize('A5');

elseif strcmp( paper,'A5 portrait')==1
    [width_cm,height_cm] = ABCsize('A5');

elseif strcmp( paper,'A4 landscape')==1
    [height_cm,width_cm] = ABCsize('A4');

elseif strcmp( paper,'A4 portrait')==1
    [width_cm,height_cm] = ABCsize('A4');

elseif strcmp( paper,'A3 landscape')==1
    [height_cm,width_cm] = ABCsize('A3');

elseif strcmp( paper,'A3 portrait')==1
    [width_cm,height_cm] = ABCsize('A3');

elseif strcmp( paper,'A2 landscape')==1
    [height_cm,width_cm] = ABCsize('A2');

elseif strcmp( paper,'A2 portrait')==1
    [width_cm,height_cm] = ABCsize('A2');

elseif strcmp( paper,'A1 landscape')==1
    [height_cm,width_cm] = ABCsize('A1');

elseif strcmp( paper,'A1 portrait')==1
    [width_cm,height_cm] = ABCsize('A1');

elseif strcmp( paper,'A0 landscape')==1
    [height_cm,width_cm] = ABCsize('A0');

elseif strcmp( paper,'A0 portrait')==1
    [width_cm,height_cm] = ABCsize('A0');

elseif strcmp( paper,'6" x 4"')==1
    width_cm = 15.2;
    height_cm = 10.2;

elseif strcmp( paper,'7" x 5"')==1
    width_cm = 17.8;
    height_cm = 12.7;

elseif strcmp( paper,'8" x 6"')==1
    width_cm = 20.3;
    height_cm = 15.2;

elseif strcmp( paper,'10" x 12"')==1
    width_cm = 30.5;
    height_cm = 25.4;

elseif strcmp( paper,'10" x 15"')==1
    width_cm = 38.1;
    height_cm = 25.4;
else
    width_cm = NaN;
    height_cm = NaN;
end
%%

% Portrait sizes of A,B,C names

function [width_cm,height_cm] = ABCsize(XN)

%Break up paper size string XN
X = XN(1);
N = str2num( XN(2:end) );

%Compute size of paper /m
if strcmp(X,'A')
    width=2^(-0.25-N/2);
    height=2^(0.25-N/2);
elseif strcmp(X,'B')
    width=2^(-N/2);
    height=2^(0.5-N/2);
elseif strcmp(X,'C')
    width=2^(-1/8-N/2);
    height=2^(3/8-N/2);
else
    width=NaN;
    height=NaN;
end

%Convert to centimeters
width_cm=width*100;
height_cm=height*100;

%End of code
