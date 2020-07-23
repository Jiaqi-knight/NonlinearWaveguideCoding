%getax
% Get non read only axes properties. Use this function to set another axes
% to look exactly the same as another via set( new_ax, f, c);
%
% f is a cell aray fo non-read only field names
% c is the corresponding cell aray of properties.

function [f,c] = getax(ax)

a = get(ax);
f = fieldnames(a);
c = get(ax,f);

%Find indices of fields that should not be cloned
do_not = [];
for n=1:length(f)
    if strcmp( f{n},'BeingDeleted' )==1
        do_not = [do_not,n];
    elseif strcmp( f{n},'Type' )==1
        do_not = [do_not,n];
    elseif strcmp( f{n},'CurrentPoint' )==1
        do_not = [do_not,n];
    elseif strcmp( f{n},'TightInset' )==1
        do_not = [do_not,n];
    elseif strcmp( f{n},'Children' )==1
        do_not = [do_not,n];
    elseif strcmp( f{n},'Tag' )==1
        do_not = [do_not,n];
    elseif strcmp( f{n},'Parent' )==1
        do_not = [do_not,n];
    end
end

%Remove these fields from f and c
f(do_not) = [];
c(do_not) = [];

%Only choose camera axes properties
ok = [];
for n=1:length(f)
    if strcmp(f{n},'CameraPosition')==1
        ok = [ok,n];
    elseif strcmp(f{n},'CameraPositionMode')==1
        ok = [ok,n];
    elseif strcmp(f{n},'CameraTarget')==1
        ok = [ok,n];
    elseif strcmp(f{n},'CameraTargetMode')==1
        ok = [ok,n];
    elseif strcmp(f{n},'CameraUpVector')==1
        ok = [ok,n];
    elseif strcmp(f{n},'CameraUpVectorMode')==1
        ok = [ok,n];
    elseif strcmp(f{n},'CameraViewAngle')==1
        ok = [ok,n];
    elseif strcmp(f{n},'CameraViewAngleMode')==1
        ok = [ok,n];
    end
end
f = f(ok);
c = c(ok);

%End of code