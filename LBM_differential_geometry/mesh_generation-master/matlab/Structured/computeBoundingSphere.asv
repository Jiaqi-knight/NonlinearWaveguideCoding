function[region]=computeBoundingSphere(region)
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

meanX = region.c1(1) + region.c2(1) + region.c3(1) + region.c4(1) + region.c5(1) + region.c6(1) + region.c7(1) + region.c8(1);
meanY = region.c1(2) + region.c2(2) + region.c3(2) + region.c4(2) + region.c5(2) + region.c6(2) + region.c7(2) + region.c8(2);
meanZ = region.c1(3) + region.c2(3) + region.c3(3) + region.c4(3) + region.c5(3) + region.c6(3) + region.c7(3) + region.c8(3);
for i=1:length(region.e1)
   meanX = meanX + region.e1(i,1);
   meanY = meanY + region.e1(i,2);
   meanZ = meanZ + region.e1(i,3);
end
for i=1:length(region.e2)
   meanX = meanX + region.e2(i,1);
   meanY = meanY + region.e2(i,2);
   meanZ = meanZ + region.e2(i,3);
end
for i=1:length(region.e3)
   meanX = meanX + region.e3(i,1);
   meanY = meanY + region.e3(i,2);
   meanZ = meanZ + region.e3(i,3);
end
for i=1:length(region.e4)
   meanX = meanX + region.e4(i,1);
   meanY = meanY + region.e4(i,2);
   meanZ = meanZ + region.e4(i,3);
end
for i=1:length(region.e5)
   meanX = meanX + region.e5(i,1);
   meanY = meanY + region.e5(i,2);
   meanZ = meanZ + region.e5(i,3);
end
for i=1:length(region.e6)
   meanX = meanX + region.e6(i,1);
   meanY = meanY + region.e6(i,2);
   meanZ = meanZ + region.e6(i,3);
end
for i=1:length(region.e7)
   meanX = meanX + region.e7(i,1);
   meanY = meanY + region.e7(i,2);
   meanZ = meanZ + region.e7(i,3);
end
for i=1:length(region.e8)
   meanX = meanX + region.e8(i,1);
   meanY = meanY + region.e8(i,2);
   meanZ = meanZ + region.e8(i,3);
end
for i=1:length(region.e9)
   meanX = meanX + region.e9(i,1);
   meanY = meanY + region.e9(i,2);
   meanZ = meanZ + region.e9(i,3);
end
for i=1:length(region.e10)
   meanX = meanX + region.e10(i,1);
   meanY = meanY + region.e10(i,2);
   meanZ = meanZ + region.e10(i,3);
end
for i=1:length(region.e11)
   meanX = meanX + region.e11(i,1);
   meanY = meanY + region.e11(i,2);
   meanZ = meanZ + region.e11(i,3);
end
for i=1:length(region.e12)
   meanX = meanX + region.e12(i,1);
   meanY = meanY + region.e12(i,2);
   meanZ = meanZ + region.e12(i,3);
end
for i=1:length(region.f1)
   meanX = meanX + region.f1(i,1);
   meanY = meanY + region.f1(i,2);
   meanZ = meanZ + region.f1(i,3);
end
for i=1:length(region.f2)
   meanX = meanX + region.f2(i,1);
   meanY = meanY + region.f2(i,2);
   meanZ = meanZ + region.f2(i,3);
end
for i=1:length(region.f3)
   meanX = meanX + region.f3(i,1);
   meanY = meanY + region.f3(i,2);
   meanZ = meanZ + region.f3(i,3);
end
for i=1:length(region.f4)
   meanX = meanX + region.f4(i,1);
   meanY = meanY + region.f4(i,2);
   meanZ = meanZ + region.f4(i,3);
end
for i=1:length(region.f5)
   meanX = meanX + region.f5(i,1);
   meanY = meanY + region.f5(i,2);
   meanZ = meanZ + region.f5(i,3);
end
for i=1:length(region.f6)
   meanX = meanX + region.f6(i,1);
   meanY = meanY + region.f6(i,2);
   meanZ = meanZ + region.f6(i,3);
end
meanX = meanX/(8+length(region.e1)+length(region.e2)+length(region.e3)+length(region.e4)+length(region.e5)+length(region.e6)+length(region.e7)+length(region.e8)+length(region.e9)+length(region.e10)+length(region.e11)+length(region.e12)+length(region.f1)+length(region.f2)+length(region.f3)+length(region.f4)+length(region.f5)+length(region.f6));
meanY = meanY/(8+length(region.e1)+length(region.e2)+length(region.e3)+length(region.e4)+length(region.e5)+length(region.e6)+length(region.e7)+length(region.e8)+length(region.e9)+length(region.e10)+length(region.e11)+length(region.e12)+length(region.f1)+length(region.f2)+length(region.f3)+length(region.f4)+length(region.f5)+length(region.f6));
meanZ = meanZ/(8+length(region.e1)+length(region.e2)+length(region.e3)+length(region.e4)+length(region.e5)+length(region.e6)+length(region.e7)+length(region.e8)+length(region.e9)+length(region.e10)+length(region.e11)+length(region.e12)+length(region.f1)+length(region.f2)+length(region.f3)+length(region.f4)+length(region.f5)+length(region.f6));

region.boundingSphere.center = [meanX meanY meanZ];


radius = 0;

if computeDistance3D(region.boundingSphere.center,region.c1)>radius
    radius = computeDistance3D(region.boundingSphere.center,region.c1);
end
if computeDistance3D(region.boundingSphere.center,region.c2)>radius
    radius = computeDistance3D(region.boundingSphere.center,region.c2);
end
if computeDistance3D(region.boundingSphere.center,region.c3)>radius
    radius = computeDistance3D(region.boundingSphere.center,region.c3);
end
if computeDistance3D(region.boundingSphere.center,region.c4)>radius
    radius = computeDistance3D(region.boundingSphere.center,region.c4);
end
if computeDistance3D(region.boundingSphere.center,region.c5)>radius
    radius = computeDistance3D(region.boundingSphere.center,region.c5);
end
if computeDistance3D(region.boundingSphere.center,region.c6)>radius
    radius = computeDistance3D(region.boundingSphere.center,region.c6);
end
if computeDistance3D(region.boundingSphere.center,region.c7)>radius
    radius = computeDistance3D(region.boundingSphere.center,region.c7);
end
if computeDistance3D(region.boundingSphere.center,region.c8)>radius
    radius = computeDistance3D(region.boundingSphere.center,region.c8);
end
for i=1:length(region.e1)
    if computeDistance3D(region.boundingSphere.center,region.e1(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e1(i));
    end
end
for i=1:length(region.e2)
    if computeDistance3D(region.boundingSphere.center,region.e2(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e2(i));
    end
end
for i=1:length(region.e3)
    if computeDistance3D(region.boundingSphere.center,region.e3(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e3(i));
    end
end
for i=1:length(region.e4)
    if computeDistance3D(region.boundingSphere.center,region.e4(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e4(i));
    end
end
for i=1:length(region.e5)
    if computeDistance3D(region.boundingSphere.center,region.e5(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e5(i));
    end
end
for i=1:length(region.e6)
    if computeDistance3D(region.boundingSphere.center,region.e6(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e6(i));
    end
end
for i=1:length(region.e7)
    if computeDistance3D(region.boundingSphere.center,region.e7(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e7(i));
    end
end
for i=1:length(region.e8)
    if computeDistance3D(region.boundingSphere.center,region.e8(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e8(i));
    end
end
for i=1:length(region.e9)
    if computeDistance3D(region.boundingSphere.center,region.e9(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e9(i));
    end
end
for i=1:length(region.e10)
    if computeDistance3D(region.boundingSphere.center,region.e10(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e10(i));
    end
end
for i=1:length(region.e11)
    if computeDistance3D(region.boundingSphere.center,region.e11(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e11(i));
    end
end
for i=1:length(region.e12)
    if computeDistance3D(region.boundingSphere.center,region.e12(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.e12(i));
    end
end
for i=1:length(region.f1)
    if computeDistance3D(region.boundingSphere.center,region.f1(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.f1(i));
    end
end
for i=1:length(region.f2)
    if computeDistance3D(region.boundingSphere.center,region.f2(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.f2(i));
    end
end
for i=1:length(region.f3)
    if computeDistance3D(region.boundingSphere.center,region.f3(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.f3(i));
    end
end
for i=1:length(region.f4)
    if computeDistance3D(region.boundingSphere.center,region.f4(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.f4(i));
    end
end
for i=1:length(region.f5)
    if computeDistance3D(region.boundingSphere.center,region.f5(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.f5(i));
    end
end
for i=1:length(region.f6)
    if computeDistance3D(region.boundingSphere.center,region.f6(i))>radius
        radius = computeDistance3D(region.boundingSphere.center,region.f6(i));
    end
end

region.boundingSphere.radius = radius;

return