function[region]=computeBoundingCircle(region)
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
%  A function to compute the bounding circle of an 2D region
%
%  Output: 
%
%%

meanX = region.c1(1) + region.c2(1) + region.c3(1) + region.c4(1);
meanY = region.c1(2) + region.c2(2) + region.c3(2) + region.c4(2);
for i=1:length(region.e1)
   meanX = meanX + region.e1(i,1);
   meanY = meanY + region.e1(i,2);
end
for i=1:length(region.e2)
   meanX = meanX + region.e2(i,1);
   meanY = meanY + region.e2(i,2);
end
for i=1:length(region.e3)
   meanX = meanX + region.e3(i,1);
   meanY = meanY + region.e3(i,2);
end
for i=1:length(region.e4)
   meanX = meanX + region.e4(i,1);
   meanY = meanY + region.e4(i,2);
end
meanX = meanX/(4+length(region.e1)+length(region.e2)+length(region.e3)+length(region.e4));
meanY = meanY/(4+length(region.e1)+length(region.e2)+length(region.e3)+length(region.e4));

region.boundingCircle.center = [meanX meanY];
radius = 0;

if computeDistance2D(region.boundingCircle.center,region.c1)>radius
    radius = computeDistance2D(region.boundingCircle.center,region.c1);
end
if computeDistance2D(region.boundingCircle.center,region.c2)>radius
    radius = computeDistance2D(region.boundingCircle.center,region.c2);
end
if computeDistance2D(region.boundingCircle.center,region.c3)>radius
    radius = computeDistance2D(region.boundingCircle.center,region.c3);
end
if computeDistance2D(region.boundingCircle.center,region.c4)>radius
    radius = computeDistance2D(region.boundingCircle.center,region.c4);
end
for i=1:length(region.e1)
    if computeDistance2D(region.boundingCircle.center,region.e1(i))>radius
        radius = computeDistance2D(region.boundingCircle.center,region.e1(i));
    end
end
for i=1:length(region.e2)
    if computeDistance2D(region.boundingCircle.center,region.e2(i))>radius
        radius = computeDistance2D(region.boundingCircle.center,region.e2(i));
    end
end
for i=1:length(region.e3)
    if computeDistance2D(region.boundingCircle.center,region.e3(i))>radius
        radius = computeDistance2D(region.boundingCircle.center,region.e3(i));
    end
end
for i=1:length(region.e4)
    if computeDistance2D(region.boundingCircle.center,region.e4(i))>radius
        radius = computeDistance2D(region.boundingCircle.center,region.e4(i));
    end
end

region.boundingCircle.radius = radius;

return