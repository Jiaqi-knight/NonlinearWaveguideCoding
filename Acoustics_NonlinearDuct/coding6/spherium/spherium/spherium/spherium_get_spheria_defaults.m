%spherium_get_spheria_defaults
% Function which looks in a (defined) directory and constructs (i) a cell
% array of spheria generatig plugin names from the directory titles, (ii)
% stores the default parameters for each option in structured array d

function [spheria_list, d] = spherium_get_spheria_defaults( spheria_dir )

%Get directory names
D = dir( spheria_dir );
N = length(D)-2;
spheria_list = cell(1,N); %Ignore . and .. fieldnames
if N>0
    for n=1:N
        spheria_list{n} = D(n+2).name;
        eval( ['d.',spheria_list{n},'=',spheria_list{n},'_defaults;']);
    end
else
    d = [];
end

%End of code