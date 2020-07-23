function[section]=homogeneousDilation2D(section,f)
%%
%==============================================================================
% Copyright (c) 2016 Université de Lorraine & Luleå tekniska universitet
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
% Neither the name of the Université de Lorraine or Luleå tekniska universitet
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
%  A function to generate a homogeneous dilation of a section with a generic
%  curvilinear boundary (output is the new dilated boundary)
%
%  f is the multiplicative factor for the dilation
%
%%

if ~isfield(section,'boundingCircle')
    section = computeBoundingCircle(section);
end

xC = section.boundingCircle.center(1);
yC = section.boundingCircle.center(2);

section.dilatedBoundary.c1 = [xC+f*cos(atan2(section.c1(2)-yC,section.c1(1)-xC))*computeDistance2D(section.c1,section.boundingCircle.center) ...
                              yC+f*sin(atan2(section.c1(2)-yC,section.c1(1)-xC))*computeDistance2D(section.c1,section.boundingCircle.center)];
section.dilatedBoundary.c2 = [xC+f*cos(atan2(section.c2(2)-yC,section.c2(1)-xC))*computeDistance2D(section.c2,section.boundingCircle.center) ...
                              yC+f*sin(atan2(section.c2(2)-yC,section.c2(1)-xC))*computeDistance2D(section.c2,section.boundingCircle.center)];
section.dilatedBoundary.c3 = [xC+f*cos(atan2(section.c3(2)-yC,section.c3(1)-xC))*computeDistance2D(section.c3,section.boundingCircle.center) ...
                              yC+f*sin(atan2(section.c3(2)-yC,section.c3(1)-xC))*computeDistance2D(section.c3,section.boundingCircle.center)];
section.dilatedBoundary.c4 = [xC+f*cos(atan2(section.c4(2)-yC,section.c4(1)-xC))*computeDistance2D(section.c4,section.boundingCircle.center) ...
                              yC+f*sin(atan2(section.c4(2)-yC,section.c4(1)-xC))*computeDistance2D(section.c4,section.boundingCircle.center)];

section.dilatedBoundary.e1 = zeros(length(section.e1),2);
for i=1:length(section.e1)
    section.dilatedBoundary.e1(i,:) = [xC+f*cos(atan2(section.c1(2)-yC,section.c1(1)-xC))*computeDistance2D(section.c1,section.boundingCircle.center) ...
                              yC+f*sin(atan2(section.c1(2)-yC,section.c1(1)-xC))*computeDistance2D(section.c1,section.boundingCircle.center)];
end
section.dilatedBoundary.e2 = zeros(length(section.e2),2);
for i=1:length(section.e2)
    
end
section.dilatedBoundary.e3 = zeros(length(section.e3),2);
for i=1:length(section.e3)
    
end
section.dilatedBoundary.e4 = zeros(length(section.e4),2);
for i=1:length(section.e4)
    
end

return