function[basepath,folder,configname]=get_inputpath()

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Zürich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 9th, 2014
%    Last update: July 9th, 2014
%
%    Description: 
%          Input: 
%         Output: 

%%

basepath = '/home/lucadistasio/Documents/ETH/Research_material/Numerics';
folder = 'Numerical_experiments';
configname = 'config';

fprintf('The basepath is set to %s. Do you want to change it? (yes/no) ',basepath)
basepathchange = input('','s');
fprintf('\n')
if isempty(basepathchange)
    basepathchange = 'no';
end
if strcmp(basepathchange,'yes') || strcmp(basepathchange,'Yes') || strcmp(basepathchange,'YES') || strcmp(basepathchange,'yep') || strcmp(basepathchange,'Y') || strcmp(basepathchange,'y') || strcmp(basepathchange,'ye') || strcmp(basepathchange,'yeah')  || strcmp(basepathchange,'ok') || strcmp(basepathchange,'OK')  || strcmp(basepathchange,'Ok')
    fprintf('Input the new basepath: ')
    basepath = input('','s');
    fprintf('\n')
end

fprintf('The current folder is %s. Do you want to change it? (yes/no) ',folder)
folderchange = input('','s');
fprintf('\n')
if isempty(folderchange)
    folderchange = 'no';
end
if strcmp(folderchange,'yes') || strcmp(folderchange,'Yes') || strcmp(folderchange,'YES') || strcmp(folderchange,'yep') || strcmp(folderchange,'Y') || strcmp(folderchange,'y') || strcmp(folderchange,'ye') || strcmp(folderchange,'yeah')  || strcmp(folderchange,'ok') || strcmp(folderchange,'OK')  || strcmp(folderchange,'Ok')
    fprintf('Input the new folder: ')
    folder = input('','s');
    fprintf('\n')
end

fprintf('The current configuration file is %s. Do you want to change it? (yes/no) ',configname)
confignamechange = input('','s');
fprintf('\n')
if isempty(confignamechange)
    confignamechange = 'no';
end
if strcmp(confignamechange,'yes') || strcmp(confignamechange,'Yes') || strcmp(confignamechange,'YES') || strcmp(confignamechange,'yep') || strcmp(confignamechange,'Y') || strcmp(confignamechange,'y') || strcmp(confignamechange,'ye') || strcmp(confignamechange,'yeah')  || strcmp(confignamechange,'ok') || strcmp(confignamechange,'OK')  || strcmp(confignamechange,'Ok')
    fprintf('Input the new configuration file: ')
    configname = input('','s');
    fprintf('\n')
end

return