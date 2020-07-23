function[region]=computeBoundingHex(region)
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

region.boundingHex.center = [meanX meanY meanZ];

lx = 0;
ly = 0;
lz = 0;

if abs(region.boundingHex.center(1)-region.c1(1))>lx
    lx = abs(region.boundingHex.center(1)-region.c1(1));
end
if abs(region.boundingHex.center(1)-region.c2(1))>lx
    lx = abs(region.boundingHex.center(1)-region.c2(1));
end
if abs(region.boundingHex.center(1)-region.c3(1))>lx
    lx = abs(region.boundingHex.center(1)-region.c3(1));
end
if abs(region.boundingHex.center(1)-region.c4(1))>lx
    lx = abs(region.boundingHex.center(1)-region.c4(1));
end
if abs(region.boundingHex.center(1)-region.c5(1))>lx
    lx = abs(region.boundingHex.center(1)-region.c5(1));
end
if abs(region.boundingHex.center(1)-region.c6(1))>lx
    lx = abs(region.boundingHex.center(1)-region.c6(1));
end
if abs(region.boundingHex.center(1)-region.c7(1))>lx
    lx = abs(region.boundingHex.center(1)-region.c7(1));
end
if abs(region.boundingHex.center(1)-region.c8(1))>lx
    lx = abs(region.boundingHex.center(1)-region.c8(1));
end

if abs(region.boundingHex.center(2)-region.c1(2))>ly
    ly = abs(region.boundingHex.center(2)-region.c1(2));
end
if abs(region.boundingHex.center(2)-region.c2(2))>ly
    ly = abs(region.boundingHex.center(2)-region.c2(2));
end
if abs(region.boundingHex.center(2)-region.c3(2))>ly
    ly = abs(region.boundingHex.center(2)-region.c3(2));
end
if abs(region.boundingHex.center(2)-region.c4(2))>ly
    ly = abs(region.boundingHex.center(2)-region.c4(2));
end
if abs(region.boundingHex.center(2)-region.c5(2))>ly
    ly = abs(region.boundingHex.center(2)-region.c5(2));
end
if abs(region.boundingHex.center(2)-region.c6(2))>ly
    ly = abs(region.boundingHex.center(2)-region.c6(2));
end
if abs(region.boundingHex.center(2)-region.c7(2))>ly
    ly = abs(region.boundingHex.center(2)-region.c7(2));
end
if abs(region.boundingHex.center(2)-region.c8(2))>ly
    ly = abs(region.boundingHex.center(2)-region.c8(2));
end

if abs(region.boundingHex.center(3)-region.c1(3))>lz
    lz = abs(region.boundingHex.center(3)-region.c1(3));
end
if abs(region.boundingHex.center(3)-region.c2(3))>lz
    lz = abs(region.boundingHex.center(3)-region.c2(3));
end
if abs(region.boundingHex.center(3)-region.c3(3))>lz
    lz = abs(region.boundingHex.center(3)-region.c3(3));
end
if abs(region.boundingHex.center(3)-region.c4(3))>lz
    lz = abs(region.boundingHex.center(3)-region.c4(3));
end
if abs(region.boundingHex.center(3)-region.c5(3))>lz
    lz = abs(region.boundingHex.center(3)-region.c5(3));
end
if abs(region.boundingHex.center(3)-region.c6(3))>lz
    lz = abs(region.boundingHex.center(3)-region.c6(3));
end
if abs(region.boundingHex.center(3)-region.c7(3))>lz
    lz = abs(region.boundingHex.center(3)-region.c7(3));
end
if abs(region.boundingHex.center(3)-region.c8(3))>lz
    lz = abs(region.boundingHex.center(3)-region.c8(3));
end

for i=1:length(region.e1)
    if abs(region.boundingHex.center(1)-region.e1(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e1(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e1(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e1(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e1(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.e1(i,3));
    end
end
for i=1:length(region.e2)
    if abs(region.boundingHex.center(1)-region.e2(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e2(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e2(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e2(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e2(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.e2(i,3));
    end
end
for i=1:length(region.e3)
    if abs(region.boundingHex.center(1)-region.e3(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e3(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e3(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e3(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e3(i,3))>ly
        ly = abs(region.boundingHex.center(3)-region.e3(i,3));
    end
end
for i=1:length(region.e4)
    if abs(region.boundingHex.center(1)-region.e4(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e4(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e4(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e4(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e4(i,3))>ly
        ly = abs(region.boundingHex.center(3)-region.e4(i,3));
    end
end
for i=1:length(region.e5)
    if abs(region.boundingHex.center(1)-region.e5(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e5(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e5(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e5(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e5(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.e5(i,3));
    end
end
for i=1:length(region.e6)
    if abs(region.boundingHex.center(1)-region.e6(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e6(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e6(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e6(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e6(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.e6(i,3));
    end
end
for i=1:length(region.e7)
    if abs(region.boundingHex.center(1)-region.e7(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e7(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e7(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e7(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e7(i,3))>ly
        ly = abs(region.boundingHex.center(3)-region.e7(i,3));
    end
end
for i=1:length(region.e8)
    if abs(region.boundingHex.center(1)-region.e8(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e8(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e8(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e8(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e8(i,3))>ly
        ly = abs(region.boundingHex.center(3)-region.e8(i,3));
    end
end
for i=1:length(region.e9)
    if abs(region.boundingHex.center(1)-region.e9(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e9(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e9(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e9(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e9(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.e9(i,3));
    end
end
for i=1:length(region.e10)
    if abs(region.boundingHex.center(1)-region.e10(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e10(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e10(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e10(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e10(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.e10(i,3));
    end
end
for i=1:length(region.e11)
    if abs(region.boundingHex.center(1)-region.e11(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e11(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e11(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e11(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e11(i,3))>ly
        ly = abs(region.boundingHex.center(3)-region.e11(i,3));
    end
end
for i=1:length(region.e12)
    if abs(region.boundingHex.center(1)-region.e12(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.e12(i,1));
    end
    if abs(region.boundingHex.center(2)-region.e12(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.e12(i,2));
    end
    if abs(region.boundingHex.center(3)-region.e12(i,3))>ly
        ly = abs(region.boundingHex.center(3)-region.e12(i,3));
    end
end

for i=1:length(region.f1)
    if abs(region.boundingHex.center(1)-region.f1(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.f1(i,1));
    end
    if abs(region.boundingHex.center(2)-region.f1(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.f1(i,2));
    end
    if abs(region.boundingHex.center(3)-region.f1(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.f1(i,3));
    end
end
for i=1:length(region.f2)
    if abs(region.boundingHex.center(1)-region.f2(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.f2(i,1));
    end
    if abs(region.boundingHex.center(2)-region.f2(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.f2(i,2));
    end
    if abs(region.boundingHex.center(3)-region.f2(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.f2(i,3));
    end
end
for i=1:length(region.f3)
    if abs(region.boundingHex.center(1)-region.f3(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.f3(i,1));
    end
    if abs(region.boundingHex.center(2)-region.f3(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.f3(i,2));
    end
    if abs(region.boundingHex.center(3)-region.f3(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.f3(i,3));
    end
end
for i=1:length(region.f4)
    if abs(region.boundingHex.center(1)-region.f4(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.f4(i,1));
    end
    if abs(region.boundingHex.center(2)-region.f4(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.f4(i,2));
    end
    if abs(region.boundingHex.center(3)-region.f4(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.f4(i,3));
    end
end
for i=1:length(region.f5)
    if abs(region.boundingHex.center(1)-region.f5(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.f5(i,1));
    end
    if abs(region.boundingHex.center(2)-region.f5(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.f5(i,2));
    end
    if abs(region.boundingHex.center(3)-region.f5(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.f5(i,3));
    end
end
for i=1:length(region.f6)
    if abs(region.boundingHex.center(1)-region.f6(i,1))>lx
        lx = abs(region.boundingHex.center(1)-region.f6(i,1));
    end
    if abs(region.boundingHex.center(2)-region.f6(i,2))>ly
        ly = abs(region.boundingHex.center(2)-region.f6(i,2));
    end
    if abs(region.boundingHex.center(3)-region.f6(i,3))>lz
        lz = abs(region.boundingHex.center(3)-region.f6(i,3));
    end
end

region.boundingHex.halfSides = [lx ly lz];

return