function[format1,format2]=get_currentdate()

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 2nd, 2014
%    Last update: July 3rd, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

c = clock;

switch c(2)
    case 1
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Jan_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' January 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Jan_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' January ',num2str(c(3)));
        end
    case 2
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Feb_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' February 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Feb_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' February ',num2str(c(3)));
        end
    case 3
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Mar_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' March 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Mar_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' March ',num2str(c(3)));
        end
    case 4
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Apr_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' April 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Apr_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' April ',num2str(c(3)));
        end
    case 5
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_May_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' May 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_May_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' May ',num2str(c(3)));
        end
    case 6
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Jun_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' June 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Jun_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' June ',num2str(c(3)));
        end
    case 7
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Jul_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' July 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Jul_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' July ',num2str(c(3)));
        end
    case 8
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Aug_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' August 0',num2str(c(3)));
            
        else
            format1 = strcat(num2str(c(1)),'_Aug_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' August ',num2str(c(3)));
        end
    case 9
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Sep_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' September 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Sep_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' September ',num2str(c(3)));
        end
    case 10
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Oct_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' October 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Oct_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' October ',num2str(c(3)));
        end
    case 11
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Nov_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' November 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Nov_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' November ',num2str(c(3)));
        end
    case 12
        if c(3)<10
            format1 = strcat(num2str(c(1)),'_Dec_0',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' December 0',num2str(c(3)));
        else
            format1 = strcat(num2str(c(1)),'_Dec_',num2str(c(3)));
            format2 = strcat(num2str(c(1)),' December ',num2str(c(3)));
        end
end

return