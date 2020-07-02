%%
%==============================================================================
% Copyright (c) 2016 Universit? de Lorraine & Lule? tekniska universitet
% Author: Luca Di Stasio <luca.distasio@gmail.com>
%                        <luca.distasio@ingpec.eu>
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 
% Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
% Neither the name of the Universit? de Lorraine or Lule? tekniska universitet
% nor the names of its contributors may be used to endorse or promote products
% derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================
%
%  DESCRIPTION
%  
%  A function to perform 
%
%  Input: 
%  Output: 
%
%%

clear all
close all
clc
subfunction_path1=genpath('C:\Users\wjq\Desktop\workspace\mesh_generation-master\matlab\Structured');
subfunction_path2=genpath('C:\Users\wjq\Desktop\workspace\interpolation-master\matlab');
subfunction_path3=genpath('C:\Users\wjq\Desktop\differential_geometry-master\differential_geometry-master\matlab');
subfunction_path4=genpath('C:\Users\wjq\Desktop\workspace\io_manager\matlab');

addpath(subfunction_path1);
addpath(subfunction_path2);
addpath(subfunction_path3);
addpath(subfunction_path4);

formatOut = 'mm-dd-yy-HH-MM-SS';
logfullfile=[datestr(now,formatOut),'.log'];




D = 3;
dim1min = -5; 
dim1max = 5;
Ndim1 = 10;
dim2min = 0;
dim2max = 5;
Ndim2 = 10;
dim3min = 0;
dim3max = 10;
Ndim3 = 10;
Nlinesdim1 = 1;
Nlinesdim2 = 1;
Nlinesdim3 = 1;
edgeflag = 0;
faceflag = 0;
cellflag = 0;
printflag = 0;

mesh = generate_regularcuboids_mesh(D,dim1min,dim1max,Ndim1,dim2min,dim2max,Ndim2,dim3min,dim3max,Ndim3,Nlinesdim1,Nlinesdim2,Nlinesdim3,edgeflag,faceflag,cellflag,printflag);

changeddomain = changedomain(D,mesh,8,10,0,2*pi,0,10);

funcs = {'r*cos(theta)','r*sin(theta)','z'};
args = {'r','theta','z'};


Dshow = 3;
pointcolor = 'r';
pointdim = 2;
linecolor = 'b';
linedim = 1;
titlestring = 'Regular mesh';
xstring = 'x';
ystring = 'y';
zstring = 'z';
shownodes = true;
showlines = false;
shownodelabels = false;
showedgelabels = false;
showfacelabels = false;
showcelllabels = false;
nodelabelcolor = [0 0 0]; % black in rgb
edgelabelcolor = [1 0 0]; % red in rgb
facelabelcolor = [0 1 0]; % green in rgb
celllabelcolor = [0 0 1]; % blue in rgb
labelsize = 12;
xfigsize = 100;
yfigsize = 100;
xaxismin = -15;
xaxismax = 15;
yaxismin = -15;
yaxismax = 15;
zaxismin = -15;
zaxismax = 15;

f1 = show_mesh(Dshow,mesh,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

f3 = show_mesh(Dshow,changeddomain,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

analyticalmorph = analytictransformation(D,changeddomain,funcs);%,args

f5 = show_mesh(Dshow,analyticalmorph,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);


D = 2;
dim1min = -5; 
dim1max = 5;
Ndim1 = 100;
dim2min = 0;
dim2max = 5;
Ndim2 = 100;
dim3min = 0;
dim3max = 10;
Ndim3 = 100;
Nlinesdim1 = 1;
Nlinesdim2 = 1;
Nlinesdim3 = 1;
edgeflag = 0;
faceflag = 0;
cellflag = 0;
printflag = 0;

mesh = generate_regularcuboids_mesh(D,dim1min,dim1max,Ndim1,dim2min,dim2max,Ndim2,dim3min,dim3max,Ndim3,Nlinesdim1,Nlinesdim2,Nlinesdim3,edgeflag,faceflag,cellflag,printflag);

changeddomain = changedomain(D,mesh,8,10,0,2*pi,0,10);

funcs = {'r*cos(theta)','r*sin(theta)'};
args = {'r','theta'};


Dshow = 2;
pointcolor = 'r';
pointdim = 2;
linecolor = 'b';
linedim = 1;
titlestring = 'Regular mesh';
xstring = 'x';
ystring = 'y';
zstring = 'z';
shownodes = true;
showlines = false;
shownodelabels = false;
showedgelabels = false;
showfacelabels = false;
showcelllabels = false;
nodelabelcolor = [0 0 0]; % black in rgb
edgelabelcolor = [1 0 0]; % red in rgb
facelabelcolor = [0 1 0]; % green in rgb
celllabelcolor = [0 0 1]; % blue in rgb
labelsize = 12;
xfigsize = 1000;
yfigsize = 1000;
xaxismin = dim1min - 5;
xaxismax = dim1max + 5;
yaxismin = dim2min - 5;
yaxismax = dim2max + 5;
zaxismin = dim3min - 5;
zaxismax = dim3max + 5;

f2 = show_mesh(Dshow,mesh,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

f4 = show_mesh(Dshow,changeddomain,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

analyticalmorph = analytictransformation(D,changeddomain,funcs,args);

f6 = show_mesh(Dshow,analyticalmorph,shownodes,showlines,shownodelabels,showedgelabels,showfacelabels,showcelllabels,pointcolor,pointdim,linecolor,linedim,nodelabelcolor,edgelabelcolor,facelabelcolor,celllabelcolor,labelsize,titlestring,xstring,ystring,zstring,xfigsize,yfigsize,xaxismin,xaxismax,yaxismin,yaxismax,zaxismin,zaxismax);

