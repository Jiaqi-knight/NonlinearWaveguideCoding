function[structuralneighbours,shearneighbours,bendneighbours,firstdevneighbours]=build_neighbourhoods3D(N,Nx,Ny,Nz,periodicity,flagintbounds,indicesintbounds,typeintbounds,indicesbulk,indicesinternalbulk,indicesF1,indicesF2,indicesF3,indicesF4,indicesF5,indicesF6,indicesinternalF1,indicesinternalF2,indicesinternalF3,indicesinternalF4,indicesinternalF5,indicesinternalF6,indicesE1,indicesE2,indicesE3,indicesE4,indicesE5,indicesE6,indicesE7,indicesE8,indicesE9,indicesE10,indicesE11,indicesE12,indicesinternalE1,indicesinternalE2,indicesinternalE3,indicesinternalE4,indicesinternalE5,indicesinternalE6,indicesinternalE7,indicesinternalE8,indicesinternalE9,indicesinternalE10,indicesinternalE11,indicesinternalE12,indicesC1,indicesC2,indicesC3,indicesC4,indicesC5,indicesC6,indicesC7,indicesC8,indicesinternalC1,indicesinternalC2,indicesinternalC3,indicesinternalC4,indicesinternalC5,indicesinternalC6,indicesinternalC7,indicesinternalC8)

%%
%        Project: Fluid - structure interaction on deformable surfaces
%         Author: Luca Di Stasio
%    Institution: ETH Z眉rich
%                 Institute for Building Materials
% Research group: Computational Physics for Engineering Materials
%        Version: 0.1
%  Creation date: July 4th, 2014
%    Last update: July 25th, 2014
%
%    Description:
%          Input:
%         Output:
%
%                 Integer mask for periodicity:
%                 1  --> neighbour of F1 - F6
%                 2  --> neighbour of F2 - F4
%                 3  --> neighbour of F3 - F5
%                 4  --> neighbour of E1 - E11
%                 5  --> neighbour of E3 - E9
%                 6  --> neighbour of E2 - E12
%                 7  --> neighbour of E4 - E10
%                 8  --> neighbour of E6 - E8
%                 9  --> neighbour of E5 - E7
%                 10 --> neighbour of C1 - C7
%                 11 --> neighbour of C2 - C8
%                 12 --> neighbour of C3 - C5
%                 13 --> neighbour of C4 - C6
%
%                 Integer mask for internal boundary topology:
%                 1   --> hole C1
%                 2   --> hole C2
%                 3   --> hole C3
%                 4   --> hole C4
%                 5   --> hole C5
%                 6   --> hole C6
%                 7   --> hole C7
%                 8   --> hole C8
%                 9   --> hole E1
%                 10  --> hole E2
%                 11  --> hole E3
%                 12  --> hole E4
%                 13  --> hole E5
%                 14  --> hole E6
%                 15  --> hole E7
%                 16  --> hole E8
%                 17  --> hole E9
%                 18  --> hole E10
%                 19  --> hole E11
%                 20  --> hole E12
%                 21  --> F1
%                 22  --> F2
%                 23  --> F3
%                 24  --> F4
%                 25  --> F5
%                 26  --> F6
%                 27  --> C1
%                 28  --> C2
%                 29  --> C3
%                 30  --> C4
%                 31  --> C5
%                 32  --> C6
%                 33  --> C7
%                 34  --> C8
%                 35  --> E1
%                 36  --> E2
%                 37  --> E3
%                 38  --> E4
%                 39  --> E5
%                 40  --> E6
%                 41  --> E7
%                 42  --> E8
%                 43  --> E9
%                 44  --> E10
%                 45  --> E11
%                 46  --> E12

%%

structuralneighbours = [zeros(N,1) -1*ones(N,6)];

shearneighbours = [zeros(N,1) -1*ones(N,20)];

bendneighbours = [zeros(N,1) -1*ones(N,6)];

firstdevneighbours = zeros(N,12);

% ---> structural neighbours;上下左右前后

structuralneighbours(indicesbulk,:) = [6*ones(length(indicesbulk),1) indicesbulk-1          indicesbulk+1          indicesbulk-Nx         indicesbulk+Nx         indicesbulk-Nx*Ny       indicesbulk+Nx*Ny    ];

if any(periodicity==1) %TO_DO；other point
    structuralneighbours(indicesF1,:)   = [6*ones(length(indicesF1),1) indicesF1-1            indicesF1+1            indicesF1-Nx           indicesF1+Nx           indicesF1+(Nz-1)*Nx*Ny    indicesF1+Nx*Ny];
    structuralneighbours(indicesF6,:)   = [6*ones(length(indicesF6),1) indicesF6-1            indicesF6+1            indicesF6-Nx           indicesF6+Nx           indicesF6-Nx*Ny           indicesF6-(Nz-1)*Nx*Ny];
else
    structuralneighbours(indicesF1,:)   = [5*ones(length(indicesF1),1) indicesF1-1            indicesF1+1            indicesF1-Nx           indicesF1+Nx           -1*ones(length(indicesF1),1)    indicesF1+Nx*Ny      ];
    structuralneighbours(indicesF6,:)   = [5*ones(length(indicesF6),1) indicesF6-1            indicesF6+1            indicesF6-Nx           indicesF6+Nx           indicesF6-Nx*Ny         -1*ones(length(indicesF6),1) ];
end
if any(periodicity==2)
    structuralneighbours(indicesF2,:)   = [6*ones(length(indicesF2),1) indicesF2-1            indicesF2+1            indicesF2+(Ny-1)*Nx   indicesF2+Nx           indicesF2-Nx*Ny         indicesF2+Nx*Ny      ];
    structuralneighbours(indicesF4,:)   = [6*ones(length(indicesF4),1) indicesF4-1            indicesF4+1            indicesF4-Nx          indicesF4-(Ny-1)*Nx    indicesF4-Nx*Ny         indicesF4+Nx*Ny      ];

else
    structuralneighbours(indicesF2,:)   = [5*ones(length(indicesF2),1) indicesF2-1            indicesF2+1            -1*ones(length(indicesF2),1)   indicesF2+Nx           indicesF2-Nx*Ny         indicesF2+Nx*Ny      ];
    structuralneighbours(indicesF4,:)   = [5*ones(length(indicesF4),1) indicesF4-1            indicesF4+1            indicesF4-Nx           -1*ones(length(indicesF4),1)   indicesF4-Nx*Ny         indicesF4+Nx*Ny      ];
end
if any(periodicity==3)
    structuralneighbours(indicesF3,:)   = [6*ones(length(indicesF3),1) indicesF3-1            indicesF3-(Nx-1)   indicesF3-Nx           indicesF3+Nx           indicesF3-Nx*Ny         indicesF3+Nx*Ny      ];
    structuralneighbours(indicesF5,:)   = [6*ones(length(indicesF5),1) indicesF5+(Nx-1)       indicesF5+1            indicesF5-Nx           indicesF5+Nx           indicesF5-Nx*Ny         indicesF5+Nx*Ny      ];
else
    structuralneighbours(indicesF3,:)   = [5*ones(length(indicesF3),1) indicesF3-1            -1*ones(length(indicesF3),1)   indicesF3-Nx           indicesF3+Nx           indicesF3-Nx*Ny         indicesF3+Nx*Ny      ];
    structuralneighbours(indicesF5,:)   = [5*ones(length(indicesF5),1) -1*ones(length(indicesF5),1)   indicesF5+1            indicesF5-Nx           indicesF5+Nx           indicesF5-Nx*Ny         indicesF5+Nx*Ny      ];
end
%bug
if any(periodicity==1) && any(periodicity==2) && any(periodicity==3) %三面均对称（for ball）
    structuralneighbours(indicesE1,:)   = [6*ones(length(indicesE1),1)  indicesE1-1            indicesE1+1            indicesE1+(Ny-1)*Nx    indicesE1+Nx           indicesE1+(Nz-1)*Nx*Ny    indicesE1+Nx*Ny      ];
    structuralneighbours(indicesE2,:)   = [6*ones(length(indicesE2),1)  indicesE2-1            indicesE2-(Nx-1)       indicesE2-Nx           indicesE2+Nx           indicesE2+(Nz-1)*Nx*Ny    indicesE2+Nx*Ny      ];
    structuralneighbours(indicesE3,:)   = [6*ones(length(indicesE3),1)  indicesE3-1            indicesE3+1            indicesE3-Nx           indicesE3-(Ny-1)*Nx    indicesE3+(Nz-1)*Nx*Ny    indicesE3+Nx*Ny      ];
    structuralneighbours(indicesE4,:)   = [6*ones(length(indicesE4),1)  indicesE4+(Nx-1)       indicesE4+1            indicesE4-Nx           indicesE4+Nx           indicesE4+(Nz-1)*Nx*Ny    indicesE4+Nx*Ny      ];
    structuralneighbours(indicesE5,:)   = [6*ones(length(indicesE5),1)  indicesE5+(Nx-1)       indicesE5+1            indicesE5+(Ny-1)*Nx    indicesE5+Nx           indicesE5-Nx*Ny           indicesE5+Nx*Ny      ];
    structuralneighbours(indicesE6,:)   = [6*ones(length(indicesE6),1)  indicesE6-1            indicesE6-(Nx-1)       indicesE6+(Ny-1)*Nx    indicesE6+Nx           indicesE6-Nx*Ny           indicesE6+Nx*Ny      ];
    structuralneighbours(indicesE7,:)   = [6*ones(length(indicesE7),1)  indicesE7-1            indicesE7-(Nx-1)       indicesE7-Nx           indicesE7-(Ny-1)*Nx    indicesE7-Nx*Ny           indicesE7+Nx*Ny      ];
    structuralneighbours(indicesE8,:)   = [6*ones(length(indicesE8),1)  indicesE8+(Nx-1)       indicesE8+1            indicesE8-Nx           indicesE8-(Ny-1)*Nx    indicesE8-Nx*Ny           indicesE8+Nx*Ny      ];
    structuralneighbours(indicesE9,:)   = [6*ones(length(indicesE9),1)  indicesE9-1            indicesE9+1            indicesE9+(Ny-1)*Nx    indicesE9+Nx           indicesE9-Nx*Ny           indicesE9-(Nz-1)*Nx*Ny ];
    structuralneighbours(indicesE10,:)  = [6*ones(length(indicesE10),1) indicesE10-1           indicesE10-(Nx-1)      indicesE10-Nx          indicesE10+Nx          indicesE10-Nx*Ny          indicesE10-(Nz-1)*Nx*Ny];
    structuralneighbours(indicesE11,:)  = [6*ones(length(indicesE11),1) indicesE11-1           indicesE11+1           indicesE11-Nx          indicesE11-(Ny-1)*Nx   indicesE11-Nx*Ny          indicesE11-(Nz-1)*Nx*Ny];
    structuralneighbours(indicesE12,:)  = [6*ones(length(indicesE12),1) indicesE12+(Nx-1)      indicesE12+1           indicesE12-Nx          indicesE12+Nx          indicesE12-Nx*Ny          indicesE12-(Nz-1)*Nx*Ny];
    
    structuralneighbours(indicesC1,:)   = [6*ones(length(indicesC1),1)  indicesC1+(Nx-1)       indicesC1+1            indicesC1+(Ny-1)*Nx    indicesC1+Nx           indicesC1+(Nz-1)*Nx*Ny    indicesC1+Nx*Ny      ];
    structuralneighbours(indicesC2,:)   = [6*ones(length(indicesC2),1)  indicesC2-1            indicesC2-(Nx-1)       indicesC2+(Ny-1)*Nx    indicesC2+Nx           indicesC2+(Nz-1)*Nx*Ny    indicesC2+Nx*Ny      ];
    structuralneighbours(indicesC3,:)   = [6*ones(length(indicesC3),1)  indicesC3-1            indicesC3-(Nx-1)       indicesC3-Nx           indicesC3-(Ny-1)*Nx    indicesC3+(Nz-1)*Nx*Ny    indicesC3+Nx*Ny      ];
    structuralneighbours(indicesC4,:)   = [6*ones(length(indicesC4),1)  indicesC4+(Nx-1)       indicesC4+1            indicesC4-Nx           indicesC4-(Ny-1)*Nx    indicesC4+(Nz-1)*Nx*Ny    indicesC4+Nx*Ny      ];
    structuralneighbours(indicesC5,:)   = [6*ones(length(indicesC5),1)  indicesC5+(Nx-1)       indicesC5+1            indicesC5+(Ny-1)*Nx    indicesC5+Nx           indicesC5-Nx*Ny           indicesC5-(Nz-1)*Nx*Ny ];
    structuralneighbours(indicesC6,:)   = [6*ones(length(indicesC6),1)  indicesC6-1            indicesC6-(Nx-1)       indicesC6+(Ny-1)*Nx    indicesC6+Nx           indicesC6-Nx*Ny           indicesC6-(Nz-1)*Nx*Ny ];
    structuralneighbours(indicesC7,:)   = [6*ones(length(indicesC7),1)  indicesC7-1            indicesC7-(Nx-1)       indicesC7-Nx           indicesC7-(Ny-1)*Nx    indicesC7-Nx*Ny           indicesC7-(Nz-1)*Nx*Ny ];
    structuralneighbours(indicesC8,:)   = [6*ones(length(indicesC8),1)  indicesC8+(Nx-1)       indicesC8+1            indicesC8-Nx           indicesC8-(Ny-1)*Nx    indicesC8-Nx*Ny           indicesC8-(Nz-1)*Nx*Ny ];
else
    structuralneighbours(indicesE1,:)   = [4*ones(length(indicesE1),1) indicesE1-1            indicesE1+1            -1*ones(length(indicesE1),1)   indicesE1+Nx           -1*ones(length(indicesE1),1)    indicesE1+Nx*Ny      ];
    structuralneighbours(indicesE2,:)   = [4*ones(length(indicesE2),1) indicesE2-1            -1*ones(length(indicesE2),1)   indicesE2-Nx           indicesE2+Nx           -1*ones(length(indicesE2),1)    indicesE2+Nx*Ny      ];
    structuralneighbours(indicesE3,:)   = [4*ones(length(indicesE3),1) indicesE3-1            indicesE3+1            indicesE3-Nx           -1*ones(length(indicesE3),1)   -1*ones(length(indicesE3),1)    indicesE3+Nx*Ny      ];
    structuralneighbours(indicesE4,:)   = [4*ones(length(indicesE4),1) -1*ones(length(indicesE4),1)   indicesE4+1            indicesE4-Nx           indicesE4+Nx           -1*ones(length(indicesE4),1)    indicesE4+Nx*Ny      ];
    structuralneighbours(indicesE5,:)   = [4*ones(length(indicesE5),1) -1*ones(length(indicesE5),1)   indicesE5+1            -1*ones(length(indicesE5),1)   indicesE5+Nx           indicesE5-Nx*Ny         indicesE5+Nx*Ny      ];
    structuralneighbours(indicesE6,:)   = [4*ones(length(indicesE6),1) indicesE6-1            -1*ones(length(indicesE6),1)   -1*ones(length(indicesE6),1)   indicesE6+Nx           indicesE6-Nx*Ny         indicesE6+Nx*Ny      ];
    structuralneighbours(indicesE7,:)   = [4*ones(length(indicesE7),1) indicesE7-1            -1*ones(length(indicesE7),1)   indicesE7-Nx           -1*ones(length(indicesE7),1)   indicesE7-Nx*Ny         indicesE7+Nx*Ny      ];
    structuralneighbours(indicesE8,:)   = [4*ones(length(indicesE8),1) -1*ones(length(indicesE8),1)   indicesE8+1            indicesE8-Nx           -1*ones(length(indicesE8),1)   indicesE8-Nx*Ny         indicesE8+Nx*Ny      ];
    structuralneighbours(indicesE9,:)   = [4*ones(length(indicesE9),1) indicesE9-1            indicesE9+1            -1*ones(length(indicesE9),1)   indicesE9+Nx           indicesE9-Nx*Ny         -1*ones(length(indicesE9),1) ];
    structuralneighbours(indicesE10,:)  = [4*ones(length(indicesE10),1) indicesE10-1           -1*ones(length(indicesE10),1)  indicesE10-Nx          indicesE10+Nx          indicesE10-Nx*Ny        -1*ones(length(indicesE10),1)];
    structuralneighbours(indicesE11,:)  = [4*ones(length(indicesE11),1) indicesE11-1           indicesE11+1           indicesE11-Nx          -1*ones(length(indicesE11),1)  indicesE11-Nx*Ny        -1*ones(length(indicesE11),1)];
    structuralneighbours(indicesE12,:)  = [4*ones(length(indicesE12),1) -1*ones(length(indicesE12),1)  indicesE12+1           indicesE12-Nx          indicesE12+Nx          indicesE12-Nx*Ny        -1*ones(length(indicesE12),1)];
    
    structuralneighbours(indicesC1,:)   = [3*ones(length(indicesC1),1) -1*ones(length(indicesC1),1)   indicesC1+1            -1*ones(length(indicesC1),1)   indicesC1+Nx           -1*ones(length(indicesC1),1)    indicesC1+Nx*Ny      ];
    structuralneighbours(indicesC2,:)   = [3*ones(length(indicesC2),1) indicesC2-1            -1*ones(length(indicesC2),1)   -1*ones(length(indicesC2),1)   indicesC2+Nx           -1*ones(length(indicesC2),1)    indicesC2+Nx*Ny      ];
    structuralneighbours(indicesC3,:)   = [3*ones(length(indicesC3),1) indicesC3-1            -1*ones(length(indicesC3),1)   indicesC3-Nx           -1*ones(length(indicesC3),1)   -1*ones(length(indicesC3),1)    indicesC3+Nx*Ny      ];
    structuralneighbours(indicesC4,:)   = [3*ones(length(indicesC4),1) -1*ones(length(indicesC4),1)   indicesC4+1            indicesC4-Nx           -1*ones(length(indicesC4),1)   -1*ones(length(indicesC4),1)    indicesC4+Nx*Ny      ];
    structuralneighbours(indicesC5,:)   = [3*ones(length(indicesC5),1) -1*ones(length(indicesC5),1)   indicesC5+1            -1*ones(length(indicesC5),1)   indicesC5+Nx           indicesC5-Nx*Ny         -1*ones(length(indicesC5),1) ];
    structuralneighbours(indicesC6,:)   = [3*ones(length(indicesC6),1) indicesC6-1            -1*ones(length(indicesC6),1)   -1*ones(length(indicesC6),1)   indicesC6+Nx           indicesC6-Nx*Ny         -1*ones(length(indicesC6),1) ];
    structuralneighbours(indicesC7,:)   = [3*ones(length(indicesC7),1) indicesC7-1            -1*ones(length(indicesC7),1)   indicesC7-Nx           -1*ones(length(indicesC7),1)   indicesC7-Nx*Ny         -1*ones(length(indicesC7),1) ];
    structuralneighbours(indicesC8,:)   = [3*ones(length(indicesC8),1) -1*ones(length(indicesC8),1)   indicesC8+1            indicesC8-Nx           -1*ones(length(indicesC8),1)   indicesC8-Nx*Ny         -1*ones(length(indicesC8),1) ];

end



if flagintbounds
    for i=1:size(indicesintbounds,1)
        index = indicesintbounds(i,1);
        switch typeintbounds(i,1)
            case 1  % hole C1
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 2  % hole C2
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 3  % hole C3
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 4  % hole C4
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 5  % hole C5
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 6  % hole C6
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 7  % hole C7
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 8  % hole C8
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 9  % hole E1
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 10 % hole E2
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 11 % hole E3
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 12 % hole E4
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 13 % hole E5
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 14 % hole E6
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 15 % hole E7
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 16 % hole E8
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 17 % hole E9
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 18 % hole E10
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 19 % hole E11
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 20 % hole E12
                structuralneighbours(index,:)   = [6    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 21 % F1
                structuralneighbours(index,:)   = [5    index-1    index+1    index-Nx    index+Nx    -1             index+Nx*Ny];
            case 22 % F2
                structuralneighbours(index,:)   = [5    index-1    index+1    -1          index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 23 % F3
                structuralneighbours(index,:)   = [5    index-1    -1         index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 24 % F4
                structuralneighbours(index,:)   = [5    index-1    index+1    index-Nx    -1          index-Nx*Ny    index+Nx*Ny];
            case 25 % F5
                structuralneighbours(index,:)   = [5    -1         index+1    index-Nx    index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 26 % F6
                structuralneighbours(index,:)   = [5    index-1    index+1    index-Nx    index+Nx    index-Nx*Ny    -1         ];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 27 % C1
                structuralneighbours(index,:)   = [3    -1         index+1    -1          index+Nx    -1             index+Nx*Ny];
            case 28 % C2
                structuralneighbours(index,:)   = [3    index-1    -1         -1          index+Nx    -1             index+Nx*Ny];
            case 29 % C3
                structuralneighbours(index,:)   = [3    index-1    -1         index-Nx    -1          -1             index+Nx*Ny];
            case 30 % C4
                structuralneighbours(index,:)   = [3    -1         index+1    index-Nx    -1          -1             index+Nx*Ny];
            case 31 % C5
                structuralneighbours(index,:)   = [3    -1         index+1    -1          index+Nx    index-Nx*Ny    -1         ];
            case 32 % C6
                structuralneighbours(index,:)   = [3    index-1    -1         -1          index+Nx    index-Nx*Ny    -1         ];
            case 33 % C7
                structuralneighbours(index,:)   = [3    index-1    -1         index-Nx    -1          index-Nx*Ny    -1         ];
            case 34 % C8
                structuralneighbours(index,:)   = [3    -1         index+1    index-Nx    -1          index-Nx*Ny    -1         ];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 35 % E1
                structuralneighbours(index,:)   = [4    index-1    index+1    -1          index+Nx    -1             index+Nx*Ny];
            case 36 % E2
                structuralneighbours(index,:)   = [4    index-1    -1         index-Nx    index+Nx    -1             index+Nx*Ny];
            case 37 % E3
                structuralneighbours(index,:)   = [4    index-1    index+1    index-Nx    -1          -1             index+Nx*Ny];
            case 38 % E4
                structuralneighbours(index,:)   = [4    -1         index+1    index-Nx    index+Nx    -1             index+Nx*Ny];
            case 39 % E5
                structuralneighbours(index,:)   = [4    -1         index+1    -1          index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 40 % E6
                structuralneighbours(index,:)   = [4    index-1    -1         -1          index+Nx    index-Nx*Ny    index+Nx*Ny];
            case 41 % E7
                structuralneighbours(index,:)   = [4    index-1    -1         index-Nx    -1          index-Nx*Ny    index+Nx*Ny];
            case 42 % E8
                structuralneighbours(index,:)   = [4    -1         index+1    index-Nx    -1          index-Nx*Ny    index+Nx*Ny];
            case 43 % E9
                structuralneighbours(index,:)   = [4    index-1    index+1    -1          index+Nx    index-Nx*Ny    -1         ];
            case 44 % E10
                structuralneighbours(index,:)   = [4    index-1    -1         index-Nx    index+Nx    index-Nx*Ny    -1         ];
            case 45 % E11
                structuralneighbours(index,:)   = [4    index-1    index+1    index-Nx    -1          index-Nx*Ny    -1         ];
            case 46 % E12
                structuralneighbours(index,:)   = [4    -1         index+1    index-Nx    index+Nx    index-Nx*Ny    -1         ];
        end
    end
end

% ---> shear neighbours

A = -1 -Nx -Nx*Ny;
B =    -Nx -Nx*Ny;
C = +1 -Nx -Nx*Ny;
D = -1     -Nx*Ny;
E = +1     -Nx*Ny;
F = -1 +Nx -Nx*Ny;
G =    +Nx -Nx*Ny;
H = +1 +Nx -Nx*Ny;
I = -1 -Nx;
L = +1 -Nx;
M = -1 +Nx;
N = +1 +Nx;
O = -1 -Nx +Nx*Ny;
P =    -Nx +Nx*Ny;
Q = +1 -Nx +Nx*Ny;
R = -1     +Nx*Ny;
S = +1     +Nx*Ny;
T = -1 +Nx +Nx*Ny;
U =    +Nx +Nx*Ny;
V = +1 +Nx +Nx*Ny;

shearneighbours(indicesbulk,:) = [20*ones(length(indicesbulk),1) indicesbulk+A indicesbulk+B indicesbulk+C indicesbulk+D indicesbulk+E indicesbulk+F indicesbulk+G indicesbulk+H indicesbulk+I indicesbulk+L ...
    indicesbulk+M indicesbulk+N indicesbulk+O indicesbulk+P indicesbulk+Q indicesbulk+R indicesbulk+S indicesbulk+T indicesbulk+U indicesbulk+V ];

shearneighbours(indicesF1,:)   = [12*ones(length(indicesF1),1) -1*ones(length(indicesF1),1) -1*ones(length(indicesF1),1) -1*ones(length(indicesF1),1) -1*ones(length(indicesF1),1) -1*ones(length(indicesF1),1) -1*ones(length(indicesF1),1) -1*ones(length(indicesF1),1) -1*ones(length(indicesF1),1) indicesF1+I indicesF1+L ...
    indicesF1+M        indicesF1+N        indicesF1+O        indicesF1+P        indicesF1+Q        indicesF1+R        indicesF1+S        indicesF1+T        indicesF1+U indicesF1+V ];
shearneighbours(indicesF2,:)   = [12*ones(length(indicesF2),1) -1*ones(length(indicesF2),1) -1*ones(length(indicesF2),1) -1*ones(length(indicesF2),1) indicesF2+D        indicesF2+E        indicesF2+F indicesF2+G indicesF2+H -1*ones(length(indicesF2),1) -1*ones(length(indicesF2),1) ...
    indicesF2+M        indicesF2+N        -1*ones(length(indicesF2),1) -1*ones(length(indicesF2),1) -1*ones(length(indicesF2),1) indicesF2+R indicesF2+S indicesF2+T indicesF2+U        indicesF2+V ];
shearneighbours(indicesF3,:)   = [12*ones(length(indicesF3),1) indicesF3+A indicesF3+B        -1*ones(length(indicesF3),1) indicesF3+D -1*ones(length(indicesF3),1) indicesF3+F indicesF3+G        -1*ones(length(indicesF3),1) indicesF3+I -1*ones(length(indicesF3),1) ...
    indicesF3+M -1*ones(length(indicesF3),1) indicesF3+O        indicesF3+P -1*ones(length(indicesF3),1) indicesF3+R -1*ones(length(indicesF3),1) indicesF3+T        indicesF3+U -1*ones(length(indicesF3),1) ];
shearneighbours(indicesF4,:)   = [12*ones(length(indicesF4),1) indicesF4+A        indicesF4+B        indicesF4+C indicesF4+D indicesF4+E -1*ones(length(indicesF4),1) -1*ones(length(indicesF4),1) -1*ones(length(indicesF4),1) indicesF4+I        indicesF4+L ...
    -1*ones(length(indicesF4),1) -1*ones(length(indicesF4),1) indicesF4+O indicesF4+P indicesF4+Q indicesF4+R        indicesF4+S        -1*ones(length(indicesF4),1) -1*ones(length(indicesF4),1) -1*ones(length(indicesF4),1) ];
shearneighbours(indicesF5,:)   = [12*ones(length(indicesF5),1) -1*ones(length(indicesF5),1) indicesF5+B indicesF5+C        -1*ones(length(indicesF5),1) indicesF5+E -1*ones(length(indicesF5),1) indicesF5+G indicesF5+H        -1*ones(length(indicesF5),1) indicesF5+L ...
    -1*ones(length(indicesF5),1) indicesF5+N -1*ones(length(indicesF5),1) indicesF5+P        indicesF5+Q -1*ones(length(indicesF5),1) indicesF5+S -1*ones(length(indicesF5),1) indicesF5+U        indicesF5+V ];
shearneighbours(indicesF6,:)   = [12*ones(length(indicesF6),1) indicesF6+A indicesF6+B indicesF6+C        indicesF6+D        indicesF6+E        indicesF6+F        indicesF6+G        indicesF6+H        indicesF6+I        indicesF6+L ...
    indicesF6+M indicesF6+N -1*ones(length(indicesF6),1) -1*ones(length(indicesF6),1) -1*ones(length(indicesF6),1) -1*ones(length(indicesF6),1) -1*ones(length(indicesF6),1) -1*ones(length(indicesF6),1) -1*ones(length(indicesF6),1) -1*ones(length(indicesF6),1) ];

shearneighbours(indicesE1,:)   = [7*ones(length(indicesE1),1) -1*ones(length(indicesE1),10) ...
    indicesE1+M indicesE1+N -1*ones(length(indicesE1),3) indicesE1+R indicesE1+S indicesE1+T indicesE1+U indicesE1+V ];
shearneighbours(indicesE2,:)   = [7*ones(length(indicesE2),1) -1*ones(length(indicesE2),8) indicesE2+I -1*ones(length(indicesE2),1) ...
    indicesE2+M -1*ones(length(indicesE2),1) indicesE2+O indicesE2+P -1*ones(length(indicesE2),1) indicesE2+R -1*ones(length(indicesE2),1) indicesE2+T indicesE2+U -1*ones(length(indicesE2),1) ];
shearneighbours(indicesE3,:)   = [7*ones(length(indicesE3),1) -1*ones(length(indicesE3),8) indicesE3+I indicesE3+L ...
    -1*ones(length(indicesE3),2) indicesE3+O indicesE3+P indicesE3+Q indicesE3+R indicesE3+S -1*ones(length(indicesE3),3) ];
shearneighbours(indicesE4,:)   = [7*ones(length(indicesE4),1) -1*ones(length(indicesE4),9) indicesE4+L ...
    -1*ones(length(indicesE4),1) indicesE4+N -1*ones(length(indicesE4),1) indicesE4+P indicesE4+Q -1*ones(length(indicesE4),1) indicesE4+S -1*ones(length(indicesE4),1) indicesE4+U indicesE4+V ];
shearneighbours(indicesE5,:)   = [7*ones(length(indicesE5),1) -1*ones(length(indicesE5),4) indicesE5+E -1*ones(length(indicesE5),1) indicesE5+G indicesE5+H -1*ones(length(indicesE5),2) ...
    -1*ones(length(indicesE5),1) indicesE5+N -1*ones(length(indicesE5),4) indicesE5+S -1*ones(length(indicesE5),1) indicesE5+U indicesE5+V ];
shearneighbours(indicesE6,:)   = [7*ones(length(indicesE6),1) -1*ones(length(indicesE6),3) indicesE6+D -1*ones(length(indicesE6),1) indicesE6+F indicesE6+G -1*ones(length(indicesE6),3) ...
    indicesE6+M -1*ones(length(indicesE6),4) indicesE6+R -1*ones(length(indicesE6),1) indicesE6+T indicesE6+U -1*ones(length(indicesE6),1) ];
shearneighbours(indicesE7,:)   = [7*ones(length(indicesE7),1) indicesE7+A indicesE7+B -1*ones(length(indicesE7),1) indicesE7+D -1*ones(length(indicesE7),4) indicesE7+I -1*ones(length(indicesE7),1) ...
    -1*ones(length(indicesE7),2) indicesE7+O indicesE7+P -1*ones(length(indicesE7),1) indicesE7+R -1*ones(length(indicesE7),4) ];
shearneighbours(indicesE8,:)   = [7*ones(length(indicesE8),1) -1*ones(length(indicesE8),1) indicesE8+B indicesE8+C -1*ones(length(indicesE8),1) indicesE8+E -1*ones(length(indicesE8),4) indicesE8+L ...
    -1*ones(length(indicesE8),3) indicesE8+P indicesE8+Q -1*ones(length(indicesE8),1) indicesE8+S -1*ones(length(indicesE8),3) ];
shearneighbours(indicesE9,:)   = [7*ones(length(indicesE9),1) -1*ones(length(indicesE9),3) indicesE9+D indicesE9+E indicesE9+F indicesE9+G indicesE9+H -1*ones(length(indicesE9),2) ...
    indicesE9+M indicesE9+N -1*ones(length(indicesE9),8) ];
shearneighbours(indicesE10,:)  = [7*ones(length(indicesE10),1) indicesE10+A indicesE10+B -1*ones(length(indicesE10),1) indicesE10+D -1*ones(length(indicesE10),1) indicesE10+F indicesE10+G -1*ones(length(indicesE10),1) indicesE10+I -1*ones(length(indicesE10),1) ...
    indicesE10+M -1*ones(length(indicesE10),9) ];
shearneighbours(indicesE11,:)  = [7*ones(length(indicesE11),1) indicesE11+A indicesE11+B indicesE11+C indicesE11+D indicesE11+E -1*ones(length(indicesE11),3) indicesE11+I indicesE11+L ...
    -1*ones(length(indicesE11),10) ];
shearneighbours(indicesE12,:)  = [7*ones(length(indicesE12),1) -1*ones(length(indicesE12),1) indicesE12+B indicesE12+C -1*ones(length(indicesE12),1) indicesE12+E -1*ones(length(indicesE12),1) indicesE12+G indicesE12+H -1*ones(length(indicesE12),1) indicesE12+L ...
    -1*ones(length(indicesE12),1) indicesE12+N -1*ones(length(indicesE12),8) ];

shearneighbours(indicesC1,:)   = [4*ones(length(indicesC1),1) -1*ones(length(indicesC1),11) indicesC1+N -1*ones(length(indicesC1),4) indicesC1+S -1*ones(length(indicesC1),1) indicesC1+U indicesC1+V ];
shearneighbours(indicesC2,:)   = [4*ones(length(indicesC2),1) -1*ones(length(indicesC2),10) ...
    indicesC2+M -1*ones(length(indicesC2),4) indicesC2+R -1*ones(length(indicesC2),1) indicesC2+T indicesC2+U -1*ones(length(indicesC2),1) ];
shearneighbours(indicesC3,:)   = [4*ones(length(indicesC3),1) -1*ones(length(indicesC3),8) indicesC3+I -1*ones(length(indicesC3),3) indicesC3+O indicesC3+P -1*ones(length(indicesC3),1) indicesC3+R -1*ones(length(indicesC3),4) ];
shearneighbours(indicesC4,:)   = [4*ones(length(indicesC4),1) -1*ones(length(indicesC4),9) indicesC4+L ...
    -1*ones(length(indicesC4),3) indicesC4+P indicesC4+Q -1*ones(length(indicesC4),1) indicesC4+S -1*ones(length(indicesC4),3) ];
shearneighbours(indicesC5,:)   = [4*ones(length(indicesC5),1) -1*ones(length(indicesC5),4) indicesC5+E -1*ones(length(indicesC5),1) indicesC5+G indicesC5+H -1*ones(length(indicesC5),3) indicesC5+N -1*ones(length(indicesC5),8) ];
shearneighbours(indicesC6,:)   = [4*ones(length(indicesC6),1) -1*ones(length(indicesC6),3) indicesC6+D -1*ones(length(indicesC6),1) indicesC6+F indicesC6+G -1*ones(length(indicesC6),3) ...
    indicesC6+M -1*ones(length(indicesC6),9) ];
shearneighbours(indicesC7,:)   = [4*ones(length(indicesC7),1) indicesC7+A indicesC7+B -1*ones(length(indicesC7),1) indicesC7+D -1*ones(length(indicesC7),4) indicesC7+I -1*ones(length(indicesC7),11) ];
shearneighbours(indicesC8,:)   = [4*ones(length(indicesC8),1) -1*ones(length(indicesC8),1) indicesC8+B indicesC8+C -1*ones(length(indicesC8),1) indicesC8+E -1*ones(length(indicesC8),4) indicesC8+L ...
    -1*ones(length(indicesC8),10) ];

if flagintbounds
    for i=1:size(indicesintbounds,1)
        index = indicesintbounds(i,1);
        switch typeintbounds(i,1)
            case 1  % hole C1
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U -1];
            case 2  % hole C2
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N index+O index+P index+Q index+R index+S -1 index+U index+V];
            case 3  % hole C3
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N -1 index+P index+Q index+R index+S index+T index+U index+V];
            case 4  % hole C4
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N index+O index+P -1 index+R index+S index+T index+U index+V];
            case 5  % hole C5
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G -1 index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
            case 6  % hole C6
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E -1 index+G index+H index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
            case 7  % hole C7
                shearneighbours(index,:)   = [20 -1 index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
            case 8  % hole C8
                shearneighbours(index,:)   = [20 index+A index+B -1 index+D index+E index+F index+G index+H index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 9  % hole E1
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N index+O index+P index+Q index+R index+S -1 -1 -1];
            case 10 % hole E2
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N -1 index+P index+Q -1 index+S -1 index+U index+V];
            case 11 % hole E3
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N -1 -1 -1 index+R index+S index+T index+U index+V];
            case 12 % hole E4
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N index+O index+P -1 index+R -1 index+T index+U -1];
            case 13 % hole E5
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E index+F index+G -1 index+I index+L index+M -1 index+O index+P index+Q index+R index+S index+T index+U -1];
            case 14 % hole E6
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E -1 index+G index+H index+I index+L -1 index+N index+O index+P index+Q index+R index+S -1 index+U index+V];
            case 15 % hole E7
                shearneighbours(index,:)   = [20 -1 index+B index+C index+D index+E index+F index+G index+H -1 index+L index+M index+N -1 index+P index+Q index+R index+S index+T index+U index+V];
            case 16 % hole E8
                shearneighbours(index,:)   = [20 index+A index+B -1 index+D index+E index+F index+G index+H index+I -1 index+M index+N index+O index+P -1 index+R index+S index+T index+U index+V];
            case 17 % hole E9
                shearneighbours(index,:)   = [20 index+A index+B index+C index+D index+E -1 -1 -1 index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
            case 18 % hole E10
                shearneighbours(index,:)   = [20 -1 index+B index+C -1 index+E -1 index+G index+H index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
            case 19 % hole E11
                shearneighbours(index,:)   = [20 -1 -1 -1 index+D index+E index+F index+G index+H index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
            case 20 % hole E12
                shearneighbours(index,:)   = [20 index+A index+B -1 index+D -1 index+F index+G -1 index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 21 % F1
                shearneighbours(index,:)   = [12 -1      -1      -1      -1      -1      -1      -1      -1      index+I index+L index+M index+N index+O index+P index+Q index+R index+S index+T index+U index+V];
            case 22 % F2
                shearneighbours(index,:)   = [12 -1      -1      -1      index+D index+E index+F index+G index+H -1      -1      index+M index+N -1      -1      -1      index+R index+S index+T index+U index+V];
            case 23 % F3
                shearneighbours(index,:)   = [12 index+A index+B -1      index+D -1      index+F index+G -1      index+I -1      index+M -1      index+O index+P -1      index+R -1      index+T index+U -1     ];
            case 24 % F4
                shearneighbours(index,:)   = [12 index+A index+B index+C index+D index+E -1      -1      -1      index+I index+L -1      -1      index+O index+P index+Q index+R index+S -1      -1      -1     ];
            case 25 % F5
                shearneighbours(index,:)   = [12 -1      index+B index+C -1      index+E -1      index+G index+H -1      index+L -1      index+N -1      index+P index+Q -1      index+S -1      index+U index+V];
            case 26 % F6
                shearneighbours(index,:)   = [12 index+A index+B index+C index+D index+E index+F index+G index+H index+I index+L index+M index+N -1      -1      -1      -1      -1      -1      -1      -1     ];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 27 % C1
                shearneighbours(index,:)   = [4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 index+N -1 -1 -1 -1 index+S -1 index+U index+V];
            case 28 % C2
                shearneighbours(index,:)   = [4 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 index+M -1 -1 -1 -1 index+R -1 index+T index+U -1];
            case 29 % C3
                shearneighbours(index,:)   = [4 -1 -1 -1 -1 -1 -1 -1 -1 index+I -1 -1 -1 index+O index+P -1 index+R -1 -1 -1 -1];
            case 30 % C4
                shearneighbours(index,:)   = [4 -1 -1 -1 -1 -1 -1 -1 -1 -1 index+L -1 -1 -1 index+P index+Q -1 index+S -1 -1 -1];
            case 31 % C5
                shearneighbours(index,:)   = [4 -1 -1 -1 -1 index+E -1 index+G index+H -1 -1 -1 index+N -1 -1 -1 -1 -1 -1 -1 -1];
            case 32 % C6
                shearneighbours(index,:)   = [4 -1 -1 -1 index+D -1 index+F index+G -1 -1 -1 index+M -1 -1 -1 -1 -1 -1 -1 -1 -1];
            case 33 % C7
                shearneighbours(index,:)   = [4 index+A index+B -1 index+D -1 -1 -1 -1 index+I -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
            case 34 % C8
                shearneighbours(index,:)   = [4 -1 index+B index+C -1 index+E -1 -1 -1 -1 index+L -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 35 % E1
                shearneighbours(index,:)   = [7 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 index+M index+N -1 -1 -1 index+R index+S index+T index+U index+V];
            case 36 % E2
                shearneighbours(index,:)   = [7 -1 -1 -1 -1 -1 -1 -1 -1 index+I -1 index+M -1 index+O index+P -1 index+R -1 index+T index+U -1];
            case 37 % E3
                shearneighbours(index,:)   = [7 -1 -1 -1 -1 -1 -1 -1 -1 index+I index+L -1 -1 index+O index+P index+Q index+R index+S -1 -1 -1];
            case 38 % E4
                shearneighbours(index,:)   = [7 -1 -1 -1 -1 -1 -1 -1 -1 -1 index+L -1 index+N -1 index+P index+Q -1 index+S -1 index+U index+V];
            case 39 % E5
                shearneighbours(index,:)   = [7 -1 -1 -1 -1 index+E -1 index+G index+H -1 -1 -1 index+N -1 -1 -1 -1 index+S -1 index+U index+V];
            case 40 % E6
                shearneighbours(index,:)   = [7 -1 -1 -1 index+D -1 index+F index+G -1 -1 -1 index+M -1 -1 -1 -1 index+R -1 index+T index+U -1];
            case 41 % E7
                shearneighbours(index,:)   = [7 index+A index+B -1 index+D -1 -1 -1 -1 index+I -1 -1 -1 index+O index+P -1 index+R -1 -1 -1 -1];
            case 42 % E8
                shearneighbours(index,:)   = [7 -1 index+B index+C -1 index+E -1 -1 -1 -1 index+L -1 -1 -1 index+P index+Q -1 index+S -1 -1 -1];
            case 43 % E9
                shearneighbours(index,:)   = [7 -1 -1 -1 index+D index+E index+F index+G index+H -1 -1 index+M index+N -1 -1 -1 -1 -1 -1 -1 -1];
            case 44 % E10
                shearneighbours(index,:)   = [7 index+A index+B -1 index+D -1 index+F index+G -1 index+I -1 index+M -1 -1 -1 -1 -1 -1 -1 -1 -1];
            case 45 % E11
                shearneighbours(index,:)   = [7 index+A index+B index+C index+D index+E -1 -1 -1 index+I index+L -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
            case 46 % E12
                shearneighbours(index,:)   = [7 -1 index+B index+C -1 index+E -1 index+G index+H -1 index+L -1 index+N -1 -1 -1 -1 -1 -1 -1 -1];
        end
    end
end

% ---> bend neighbours

bendneighbours(indicesinternalbulk,:) = [6*ones(length(indicesinternalbulk),1) indicesinternalbulk-2 indicesinternalbulk+2 indicesinternalbulk-2*Nx indicesinternalbulk+2*Nx indicesinternalbulk-2*Nx*Ny indicesinternalbulk+2*Nx*Ny];

% bendneighbours(indicesF1,:)          = [5*ones(length(indicesF1),1) indicesF1-2                  indicesF1+2                  indicesF1-2*Nx               indicesF1+2*Nx               -1*ones(length(indicesF1),1) indicesF1+2*Nx*Ny           ];
% bendneighbours(indicesF2,:)          = [5*ones(length(indicesF2),1) indicesF2-2                  indicesF2+2                  -1*ones(length(indicesF2),1) indicesF2+2*Nx               indicesF2-2*Nx*Ny            indicesF2+2*Nx*Ny           ];
% bendneighbours(indicesF3,:)          = [5*ones(length(indicesF3),1) indicesF3-2                  -1*ones(length(indicesF3),1) indicesF3-2*Nx               indicesF3+2*Nx               indicesF3-2*Nx*Ny            indicesF3+2*Nx*Ny           ];
% bendneighbours(indicesF4,:)          = [5*ones(length(indicesF4),1) indicesF4-2                  indicesF4+2                  indicesF4-2*Nx               -1*ones(length(indicesF4),1) indicesF4-2*Nx*Ny            indicesF4+2*Nx*Ny           ];
% bendneighbours(indicesF5,:)          = [5*ones(length(indicesF5),1) -1*ones(length(indicesF5),1) indicesF5+2                  indicesF5-2*Nx               indicesF5+2*Nx               indicesF5-2*Nx*Ny            indicesF5+2*Nx*Ny           ];
% bendneighbours(indicesF6,:)          = [5*ones(length(indicesF6),1) indicesF6-2                  indicesF6+2                  indicesF6-2*Nx               indicesF6+2*Nx               indicesF6-2*Nx*Ny            -1*ones(length(indicesF6),1)];
%
% bendneighbours(indicesE1,:)          = [4*ones(length(indicesE1),1) indicesE1-2                  indicesE1+2                  -1*ones(length(indicesE1),1) indicesE1+2*Nx               -1*ones(length(indicesE1),1) indicesE1+2*Nx*Ny           ];
% bendneighbours(indicesE2,:)          = [4*ones(length(indicesE2),1) indicesE2-2                  -1*ones(length(indicesE2),1) indicesE2-2*Nx               indicesE2+2*Nx               -1*ones(length(indicesE2),1) indicesE2+2*Nx*Ny           ];
% bendneighbours(indicesE3,:)          = [4*ones(length(indicesE3),1) indicesE3-2                  indicesE3+2                  indicesE3-2*Nx               -1*ones(length(indicesE3),1) -1*ones(length(indicesE3),1) indicesE3+2*Nx*Ny           ];
% bendneighbours(indicesE4,:)          = [4*ones(length(indicesE4),1) -1*ones(length(indicesE4),1) indicesE4+2                  indicesE4-2*Nx               indicesE4+2*Nx               -1*ones(length(indicesE4),1) indicesE4+2*Nx*Ny           ];
% bendneighbours(indicesE5,:)          = [4*ones(length(indicesE5),1) -1*ones(length(indicesE5),1) indicesE5+2                  -1*ones(length(indicesE5),1) indicesE5+2*Nx               indicesE5-2*Nx*Ny            indicesE5+2*Nx*Ny           ];
% bendneighbours(indicesE6,:)          = [4*ones(length(indicesE6),1) indicesE6-2                  -1*ones(length(indicesE6),1) -1*ones(length(indicesE6),1) indicesE6+2*Nx               indicesE6-2*Nx*Ny            indicesE6+2*Nx*Ny           ];
% bendneighbours(indicesE7,:)          = [4*ones(length(indicesE7),1) indicesE7-2                  -1*ones(length(indicesE7),1) indicesE7-2*Nx               -1*ones(length(indicesE7),1) indicesE7-2*Nx*Ny            indicesE7+2*Nx*Ny           ];
% bendneighbours(indicesE8,:)          = [4*ones(length(indicesE8),1) -1*ones(length(indicesE8),1) indicesE8+2                  indicesE8-2*Nx               -1*ones(length(indicesE8),1) indicesE8-2*Nx*Ny            indicesE8+2*Nx*Ny           ];
% bendneighbours(indicesE9,:)          = [4*ones(length(indicesE9),1) indicesE9-2                  indicesE9+2                  -1*ones(length(indicesE9),1) indicesE9+2*Nx               indicesE9-2*Nx*Ny            -1*ones(length(indicesE9),1)];
% bendneighbours(indicesE10,:)         = [4*ones(length(indicesE10),1) indicesE10-2                 -1*ones(length(indicesE10),1) indicesE10-2*Nx              indicesE10+2*Nx              indicesE10-2*Nx*Ny           -1*ones(length(indicesE10),1)];
% bendneighbours(indicesE11,:)         = [4*ones(length(indicesE11),1) indicesE11-2                 indicesE11+2                 indicesE11-2*Nx              -1*ones(length(indicesE11),1) indicesE11-2*Nx*Ny           -1*ones(length(indicesE11),1)];
% bendneighbours(indicesE12,:)         = [4*ones(length(indicesE12),1) -1*ones(length(indicesE12),1) indicesE12+2                 indicesE12-2*Nx              indicesE12+2*Nx              indicesE12-2*Nx*Ny           -1*ones(length(indicesE12),1)];

bendneighbours(indicesC1,:)          = [3*ones(length(indicesC1),1) -1*ones(length(indicesC1),1) indicesC1+2                  -1*ones(length(indicesC1),1) indicesC1+2*Nx               -1*ones(length(indicesC1),1) indicesC1+2*Nx*Ny           ];
bendneighbours(indicesC2,:)          = [3*ones(length(indicesC2),1) indicesC2-2                  -1*ones(length(indicesC2),1) -1*ones(length(indicesC2),1) indicesC2+2*Nx               -1*ones(length(indicesC2),1) indicesC2+2*Nx*Ny           ];
bendneighbours(indicesC3,:)          = [3*ones(length(indicesC3),1) indicesC3-2                  -1*ones(length(indicesC3),1) indicesC3-2*Nx               -1*ones(length(indicesC3),1) -1*ones(length(indicesC3),1) indicesC3+2*Nx*Ny           ];
bendneighbours(indicesC4,:)          = [3*ones(length(indicesC4),1) -1*ones(length(indicesC4),1) indicesC4+2                  indicesC4-2*Nx               -1*ones(length(indicesC4),1) -1*ones(length(indicesC4),1) indicesC4+2*Nx*Ny           ];
bendneighbours(indicesC5,:)          = [3*ones(length(indicesC5),1) -1*ones(length(indicesC5),1) indicesC5+2                  -1*ones(length(indicesC5),1) indicesC5+2*Nx               indicesC5-2*Nx*Ny            -1*ones(length(indicesC5),1)];
bendneighbours(indicesC6,:)          = [3*ones(length(indicesC6),1) indicesC6-2                  -1*ones(length(indicesC6),1) -1*ones(length(indicesC6),1) indicesC6+2*Nx               indicesC6-2*Nx*Ny            -1*ones(length(indicesC6),1)];
bendneighbours(indicesC7,:)          = [3*ones(length(indicesC7),1) indicesC7-2                  -1*ones(length(indicesC7),1) indicesC7-2*Nx               -1*ones(length(indicesC7),1) indicesC7-2*Nx*Ny            -1*ones(length(indicesC7),1)];
bendneighbours(indicesC8,:)          = [3*ones(length(indicesC8),1) -1*ones(length(indicesC8),1) indicesC8+2                  indicesC8-2*Nx               -1*ones(length(indicesC8),1) indicesC8-2*Nx*Ny            -1*ones(length(indicesC8),1)];

bendneighbours(indicesinternalF1,:)  = [5*ones(length(indicesinternalF1),1) indicesinternalF1-2                  indicesinternalF1+2                  indicesinternalF1-2*Nx               indicesinternalF1+2*Nx               -1*ones(length(indicesinternalF1),1) indicesinternalF1+2*Nx*Ny           ];
bendneighbours(indicesinternalF2,:)  = [5*ones(length(indicesinternalF2),1) indicesinternalF2-2                  indicesinternalF2+2                  -1*ones(length(indicesinternalF2),1) indicesinternalF2+2*Nx               indicesinternalF2-2*Nx*Ny            indicesinternalF2+2*Nx*Ny           ];
bendneighbours(indicesinternalF3,:)  = [5*ones(length(indicesinternalF3),1) indicesinternalF3-2                  -1*ones(length(indicesinternalF3),1) indicesinternalF3-2*Nx               indicesinternalF3+2*Nx               indicesinternalF3-2*Nx*Ny            indicesinternalF3+2*Nx*Ny           ];
bendneighbours(indicesinternalF4,:)  = [5*ones(length(indicesinternalF4),1) indicesinternalF4-2                  indicesinternalF4+2                  indicesinternalF4-2*Nx               -1*ones(length(indicesinternalF4),1) indicesinternalF4-2*Nx*Ny            indicesinternalF4+2*Nx*Ny           ];
bendneighbours(indicesinternalF5,:)  = [5*ones(length(indicesinternalF5),1) -1*ones(length(indicesinternalF5),1) indicesinternalF5+2                  indicesinternalF5-2*Nx               indicesinternalF5+2*Nx               indicesinternalF5-2*Nx*Ny            indicesinternalF5+2*Nx*Ny           ];
bendneighbours(indicesinternalF6,:)  = [5*ones(length(indicesinternalF6),1) indicesinternalF6-2                  indicesinternalF6+2                  indicesinternalF6-2*Nx               indicesinternalF6+2*Nx               indicesinternalF6-2*Nx*Ny            -1*ones(length(indicesinternalF6),1)];

bendneighbours(indicesinternalE1,:)  = [4*ones(length(indicesinternalE1),1) indicesinternalE1-2                  indicesinternalE1+2                  -1*ones(length(indicesinternalE1),1) indicesinternalE1+2*Nx               -1*ones(length(indicesinternalE1),1) indicesinternalE1+2*Nx*Ny           ];
bendneighbours(indicesinternalE2,:)  = [4*ones(length(indicesinternalE2),1) indicesinternalE2-2                  -1*ones(length(indicesinternalE2),1) indicesinternalE2-2*Nx               indicesinternalE2+2*Nx               -1*ones(length(indicesinternalE2),1) indicesinternalE2+2*Nx*Ny           ];
bendneighbours(indicesinternalE3,:)  = [4*ones(length(indicesinternalE3),1) indicesinternalE3-2                  indicesinternalE3+2                  indicesinternalE3-2*Nx               -1*ones(length(indicesinternalE3),1) -1*ones(length(indicesinternalE3),1) indicesinternalE3+2*Nx*Ny           ];
bendneighbours(indicesinternalE4,:)  = [4*ones(length(indicesinternalE4),1) -1*ones(length(indicesinternalE4),1) indicesinternalE4+2                  indicesinternalE4-2*Nx               indicesinternalE4+2*Nx               -1*ones(length(indicesinternalE4),1) indicesinternalE4+2*Nx*Ny           ];
bendneighbours(indicesinternalE5,:)  = [4*ones(length(indicesinternalE5),1) -1*ones(length(indicesinternalE5),1) indicesinternalE5+2                  -1*ones(length(indicesinternalE5),1) indicesinternalE5+2*Nx               indicesinternalE5-2*Nx*Ny            indicesinternalE5+2*Nx*Ny           ];
bendneighbours(indicesinternalE6,:)  = [4*ones(length(indicesinternalE6),1) indicesinternalE6-2                  -1*ones(length(indicesinternalE6),1) -1*ones(length(indicesinternalE6),1) indicesinternalE6+2*Nx               indicesinternalE6-2*Nx*Ny            indicesinternalE6+2*Nx*Ny           ];
bendneighbours(indicesinternalE7,:)  = [4*ones(length(indicesinternalE7),1) indicesinternalE7-2                  -1*ones(length(indicesinternalE7),1) indicesinternalE7-2*Nx               -1*ones(length(indicesinternalE7),1) indicesinternalE7-2*Nx*Ny            indicesinternalE7+2*Nx*Ny           ];
bendneighbours(indicesinternalE8,:)  = [4*ones(length(indicesinternalE8),1) -1*ones(length(indicesinternalE8),1) indicesinternalE8+2                  indicesinternalE8-2*Nx               -1*ones(length(indicesinternalE8),1) indicesinternalE8-2*Nx*Ny            indicesinternalE8+2*Nx*Ny           ];
bendneighbours(indicesinternalE9,:)  = [4*ones(length(indicesinternalE9),1) indicesinternalE9-2                  indicesinternalE9+2                  -1*ones(length(indicesinternalE9),1) indicesinternalE9+2*Nx               indicesinternalE9-2*Nx*Ny            -1*ones(length(indicesinternalE9),1)];
bendneighbours(indicesinternalE10,:) = [4*ones(length(indicesinternalE10),1) indicesinternalE10-2                 -1*ones(length(indicesinternalE10),1) indicesinternalE10-2*Nx              indicesinternalE10+2*Nx              indicesinternalE10-2*Nx*Ny           -1*ones(length(indicesinternalE10),1)];
bendneighbours(indicesinternalE11,:) = [4*ones(length(indicesinternalE11),1) indicesinternalE11-2                 indicesinternalE11+2                 indicesinternalE11-2*Nx              -1*ones(length(indicesinternalE11),1) indicesinternalE11-2*Nx*Ny           -1*ones(length(indicesinternalE11),1)];
bendneighbours(indicesinternalE12,:) = [4*ones(length(indicesinternalE12),1) -1*ones(length(indicesinternalE12),1) indicesinternalE12+2                 indicesinternalE12-2*Nx              indicesinternalE12+2*Nx              indicesinternalE12-2*Nx*Ny           -1*ones(length(indicesinternalE12),1)];

bendneighbours(indicesinternalC1,:)  = [3*ones(length(indicesinternalC1),1) -1*ones(length(indicesinternalC1),1) indicesinternalC1+2                  -1*ones(length(indicesinternalC1),1) indicesinternalC1+2*Nx               -1*ones(length(indicesinternalC1),1) indicesinternalC1+2*Nx*Ny           ];
bendneighbours(indicesinternalC2,:)  = [3*ones(length(indicesinternalC2),1) indicesinternalC2-2                  -1*ones(length(indicesinternalC2),1) -1*ones(length(indicesinternalC2),1) indicesinternalC2+2*Nx               -1*ones(length(indicesinternalC2),1) indicesinternalC2+2*Nx*Ny           ];
bendneighbours(indicesinternalC3,:)  = [3*ones(length(indicesinternalC3),1) indicesinternalC3-2                  -1*ones(length(indicesinternalC3),1) indicesinternalC3-2*Nx               -1*ones(length(indicesinternalC3),1) -1*ones(length(indicesinternalC3),1) indicesinternalC3+2*Nx*Ny           ];
bendneighbours(indicesinternalC4,:)  = [3*ones(length(indicesinternalC4),1) -1*ones(length(indicesinternalC4),1) indicesinternalC4+2                  indicesinternalC4-2*Nx               -1*ones(length(indicesinternalC4),1) -1*ones(length(indicesinternalC4),1) indicesinternalC4+2*Nx*Ny           ];
bendneighbours(indicesinternalC5,:)  = [3*ones(length(indicesinternalC5),1) -1*ones(length(indicesinternalC5),1) indicesinternalC5+2                  -1*ones(length(indicesinternalC5),1) indicesinternalC5+2*Nx               indicesinternalC5-2*Nx*Ny            -1*ones(length(indicesinternalC5),1)];
bendneighbours(indicesinternalC6,:)  = [3*ones(length(indicesinternalC6),1) indicesinternalC6-2                  -1*ones(length(indicesinternalC6),1) -1*ones(length(indicesinternalC6),1) indicesinternalC6+2*Nx               indicesinternalC6-2*Nx*Ny            -1*ones(length(indicesinternalC6),1)];
bendneighbours(indicesinternalC7,:)  = [3*ones(length(indicesinternalC7),1) indicesinternalC7-2                  -1*ones(length(indicesinternalC7),1) indicesinternalC7-2*Nx               -1*ones(length(indicesinternalC7),1) indicesinternalC7-2*Nx*Ny            -1*ones(length(indicesinternalC7),1)];
bendneighbours(indicesinternalC8,:)  = [3*ones(length(indicesinternalC8),1) -1*ones(length(indicesinternalC8),1) indicesinternalC8+2                  indicesinternalC8-2*Nx               -1*ones(length(indicesinternalC8),1) indicesinternalC8-2*Nx*Ny            -1*ones(length(indicesinternalC8),1)];

if flagintbounds
    for i=1:size(indicesintbounds,1)
        index = indicesintbounds(i,1);
        switch typeintbounds(i,1)
            case 1  % hole C1
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 2  % hole C2
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 3  % hole C3
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 4  % hole C4
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 5  % hole C5
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 6  % hole C6
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 7  % hole C7
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 8  % hole C8
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 9  % hole E1
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 10 % hole E2
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 11 % hole E3
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 12 % hole E4
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 13 % hole E5
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 14 % hole E6
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 15 % hole E7
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 16 % hole E8
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 17 % hole E9
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 18 % hole E10
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 19 % hole E11
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
            case 20 % hole E12
                bendneighbours(index,:)   = [6    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 21 % F1
                bendneighbours(index,:)    = [5    index-2    index+2    index-2*Nx    index+2*Nx    -1               index+2*Nx*Ny];
                index2 = index+Nx*Ny;
                bendneighbours(index2,:)   = [5    index2-2   index2+2   index2-2*Nx   index2+2*Nx   -1               index2+2*Nx*Ny];
            case 22 % F2
                bendneighbours(index,:)    = [5    index-2    index+2    -1            index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
                index2 = index+Nx;
                bendneighbours(index2,:)   = [5    index2-2   index2+2   -1            index2+2*Nx   index2-2*Nx*Ny   index2+2*Nx*Ny];
            case 23 % F3
                bendneighbours(index,:)    = [5    index-2    -1         index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
                index2 = index-1;
                bendneighbours(index2,:)   = [5    index2-2   -1         index2-2*Nx   index2+2*Nx   index2-2*Nx*Ny   index2+2*Nx*Ny];
            case 24 % F4
                bendneighbours(index,:)    = [5    index-2    index+2    index-2*Nx    -1            index-2*Nx*Ny    index+2*Nx*Ny];
                index2 = index-Nx;
                bendneighbours(index2,:)   = [5    index2-2   index2+2   index2-2*Nx   -1            index2-2*Nx*Ny   index2+2*Nx*Ny];
            case 25 % F5
                bendneighbours(index,:)    = [5    -1         index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    index+2*Nx*Ny];
                index2 = index+1;
                bendneighbours(index2,:)   = [5    -1         index2+2   index2-2*Nx   index2+2*Nx   index2-2*Nx*Ny   index2+2*Nx*Ny];
            case 26 % F6
                bendneighbours(index,:)    = [5    index-2    index+2    index-2*Nx    index+2*Nx    index-2*Nx*Ny    -1           ];
                index2 = index-Nx*Ny;
                bendneighbours(index2,:)   = [5    index2-2   index2+2   index2-2*Nx   index2+2*Nx   index2-2*Nx*Ny   -1           ];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 27 % C1
                B =   -Nx-Nx*Ny;
                D = -1   -Nx*Ny;
                E =      -Nx*Ny;
                L = -1-Nx      ;
                M =   -Nx      ;
                O = -1         ;
                l = index;
                indices = [l;l+B;l+D;l+E;l+L;l+M;l+O];
                bendneighbours(indices,:)   = [3*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2    -1*ones(length(indices),1)    indices+2*Nx    -1*ones(length(indices),1)    indices+2*Nx*Ny];
            case 28 % C2
                B =   -Nx-Nx*Ny;
                E =      -Nx*Ny;
                F = +1   -Nx*Ny;
                M =   -Nx      ;
                N = +1-Nx      ;
                P = +1         ;
                l = index;
                indices = [l;l+B;l+E;l+F;l+M;l+N;l+P];
                bendneighbours(indices,:)   = [3*ones(length(indices),1)    indices-2    -1*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2*Nx    -1*ones(length(indices),1)    indices+2*Nx*Ny];
            case 29 % C3
                E =      -Nx*Ny;
                F = +1   -Nx*Ny;
                H =   +Nx-Nx*Ny;
                P = +1         ;
                R =   +Nx      ;
                S = +1+Nx      ;
                l = index;
                indices = [l;l+E;l+F;l+H;l+P;l+R;l+S];
                bendneighbours(indices,:)   = [3*ones(length(indices),1)    indices-2    -1*ones(length(indices),1)    indices-2*Nx    -1*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2*Nx*Ny];
            case 30 % C4
                D = -1   -Nx*Ny;
                E =      -Nx*Ny;
                H =   +Nx-Nx*Ny;
                O = -1         ;
                Q = -1+Nx      ;
                R =   +Nx      ;
                l = index;
                indices = [l;l+D;l+E;l+H;l+O;l+Q;l+R];
                bendneighbours(indices,:)   = [3*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2    indices-2*Nx    -1*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2*Nx*Ny];
            case 31 % C5
                L = -1-Nx      ;
                M =   -Nx      ;
                O = -1         ;
                U =   -Nx+Nx*Ny;
                Z = -1   +Nx*Ny;
                K =      +Nx*Ny;
                l = index;
                indices = [l;l+L;l+M;l+O;l+U;l+K;l+Z];
                bendneighbours(indices,:)   = [3*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2    -1*ones(length(indices),1)    indices+2*Nx    indices-2*Nx*Ny    -1*ones(length(indices),1)];
            case 32 % C6
                M =   -Nx      ;
                N = +1-Nx      ;
                P = +1         ;
                U =   -Nx+Nx*Ny;
                K =      +Nx*Ny;
                J = +1   +Nx*Ny;
                l = index;
                indices = [l;l+M;l+N;l+P;l+U;l+K;l+J];
                bendneighbours(indices,:)   = [3*ones(length(indices),1)    indices-2    -1*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2*Nx    indices-2*Nx*Ny    -1*ones(length(indices),1)];
            case 33 % C7
                P = +1         ;
                R =   +Nx      ;
                S = +1+Nx      ;
                K =      +Nx*Ny;
                J = +1   +Nx*Ny;
                X =   +Nx+Nx*Ny;
                l = index;
                indices = [l;l+P;l+R;l+S;l+K;l+J;l+X];
                bendneighbours(indices,:)   = [3*ones(length(indices),1)    indices-2    -1*ones(length(indices),1)    indices-2*Nx    -1*ones(length(indices),1)    indices-2*Nx*Ny    -1*ones(length(indices),1)];
            case 34 % C8
                O = -1         ;
                Q = -1+Nx      ;
                R =   +Nx      ;
                Z = -1   +Nx*Ny;
                K =      +Nx*Ny;
                X =   +Nx+Nx*Ny;
                l = index;
                indices = [l;l+O;l+Q;l+R;l+Z;l+K;l+X];
                bendneighbours(indices,:)   = [3*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2    indices-2*Nx    -1*ones(length(indices),1)    indices-2*Nx*Ny    -1*ones(length(indices),1)];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 35 % E1
                l = index;
                indices = [l;l-Nx*Ny;l-Nx-Nx*Ny;l-Nx];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    indices-2    indices+2    -1*ones(length(indices),1)    indices+2*Nx    -1*ones(length(indices),1)    indices+2*Nx*Ny];
            case 36 % E2
                l = index;
                indices = [l;l-Nx*Ny;l+1-Nx*Ny;l+1];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    indices-2    -1*ones(length(indices),1)    indices-2*Nx    indices+2*Nx    -1*ones(length(indices),1)    indices+2*Nx*Ny];
            case 37 % E3
                l = index;
                indices = [l;l-Nx*Ny;l-Nx*Ny+Nx;l+Nx];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    indices-2    indices+2    indices-2*Nx    -1*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2*Nx*Ny];
            case 38 % E4
                l = index;
                indices = [l;l-Nx*Ny;l-Nx*Ny-1;l-1];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2    indices-2*Nx    indices+2*Nx    -1*ones(length(indices),1)    indices+2*Nx*Ny];
            case 39 % E5
                l = index;
                indices = [l;l-1;l-Nx;l-1-Nx];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2    -1*ones(length(indices),1)    indices+2*Nx    indices-2*Nx*Ny    indices+2*Nx*Ny];
            case 40 % E6
                l = index;
                indices = [l;l+1;l-Nx;l+1-Nx];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    indices-2    -1*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2*Nx    indices-2*Nx*Ny    indices+2*Nx*Ny];
            case 41 % E7
                l = index;
                indices = [l;l+1;l+Nx;l+1+Nx];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    indices-2    -1*ones(length(indices),1)    indices-2*Nx    -1*ones(length(indices),1)    indices-2*Nx*Ny    indices+2*Nx*Ny];
            case 42 % E8
                l = index;
                indices = [l;l-1;l+Nx;l-1+Nx];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2    indices-2*Nx    -1*ones(length(indices),1)    indices-2*Nx*Ny    indices+2*Nx*Ny];
            case 43 % E9
                l = index;
                indices = [l;l+Nx*Ny;l-Nx;l-Nx+Nx*Ny];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    indices-2    indices+2    -1*ones(length(indices),1)    indices+2*Nx    indices-2*Nx*Ny    -1*ones(length(indices),1)];
            case 44 % E10
                l = index;
                indices = [l;l+Nx*Ny;l+1;l+1+Nx*Ny];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    indices-2    -1*ones(length(indices),1)    indices-2*Nx    indices+2*Nx    indices-2*Nx*Ny    -1*ones(length(indices),1)];
            case 45 % E11
                l = index;
                indices = [l;l+Nx*Ny;l+Nx;l+Nx+Nx*Ny];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    indices-2    indices+2    indices-2*Nx    -1*ones(length(indices),1)    indices-2*Nx*Ny    -1*ones(length(indices),1)];
            case 46 % E12
                l = index;
                indices = [l;l+Nx*Ny;l-1;l-1+Nx*Ny];
                bendneighbours(indices,:)   = [4*ones(length(indices),1)    -1*ones(length(indices),1)    indices+2    indices-2*Nx    indices+2*Nx    indices-2*Nx*Ny    -1*ones(length(indices),1)];
        end
    end
end

% ---> first derivatives neighbourhoods

% integer mask for type of finite difference approximation (at 2nd order) for 1st derivative
% 1 --> central finite difference
% 2 --> one-sided finite difference with nodes on + side (missing on - side)
% 3 --> one-sided finite difference with nodes on - side (missing on + side)

firstdevneighbours(indicesbulk,:) = [1*ones(length(indicesbulk),1) indicesbulk-1  indicesbulk+1 zeros(length(indicesbulk),1) 1*ones(length(indicesbulk),1) indicesbulk-Nx indicesbulk+Nx zeros(length(indicesbulk),1) 1*ones(length(indicesbulk),1) indicesbulk-Nx*Ny  indicesbulk+Nx*Ny  zeros(length(indicesbulk),1)];

firstdevneighbours(indicesF1,:)   = [1*ones(length(indicesF1),1)   indicesF1-1    indicesF1+1   zeros(length(indicesF1),1)   1*ones(length(indicesF1),1)   indicesF1-Nx   indicesF1+Nx   zeros(length(indicesF1),1)   2*ones(length(indicesF1),1)   indicesF1          indicesF1+Nx*Ny    indicesF1+2*Nx*Ny           ];
firstdevneighbours(indicesF2,:)   = [1*ones(length(indicesF2),1)   indicesF2-1    indicesF2+1   zeros(length(indicesF2),1)   2*ones(length(indicesF2),1)   indicesF2      indicesF2+Nx   indicesF2+2*Nx               1*ones(length(indicesF2),1)   indicesF2-Nx*Ny    indicesF2+Nx*Ny    zeros(length(indicesF2),1)  ];
firstdevneighbours(indicesF3,:)   = [3*ones(length(indicesF3),1)   indicesF3      indicesF3-1   indicesF3-2                  1*ones(length(indicesF3),1)   indicesF3-Nx   indicesF3+Nx   zeros(length(indicesF3),1)   1*ones(length(indicesF3),1)   indicesF3-Nx*Ny    indicesF3+Nx*Ny    zeros(length(indicesF3),1)  ];
firstdevneighbours(indicesF4,:)   = [1*ones(length(indicesF4),1)   indicesF4-1    indicesF4+1   zeros(length(indicesF4),1)   3*ones(length(indicesF4),1)   indicesF4      indicesF4-Nx   indicesF4-2*Nx               1*ones(length(indicesF4),1)   indicesF4-Nx*Ny    indicesF4+Nx*Ny    zeros(length(indicesF4),1)  ];
firstdevneighbours(indicesF5,:)   = [2*ones(length(indicesF5),1)   indicesF5      indicesF5+1   indicesF5+2                  1*ones(length(indicesF5),1)   indicesF5-Nx   indicesF5+Nx   zeros(length(indicesF5),1)   1*ones(length(indicesF5),1)   indicesF5-Nx*Ny    indicesF5+Nx*Ny    zeros(length(indicesF5),1)  ];
firstdevneighbours(indicesF6,:)   = [1*ones(length(indicesF6),1)   indicesF6-1    indicesF6+1   zeros(length(indicesF6),1)   1*ones(length(indicesF6),1)   indicesF6-Nx   indicesF6+Nx   zeros(length(indicesF6),1)   3*ones(length(indicesF6),1)   indicesF6          indicesF6-Nx*Ny    indicesF6-2*Nx*Ny           ];

firstdevneighbours(indicesE1,:)   = [1*ones(length(indicesE1),1)   indicesE1-1    indicesE1+1   zeros(length(indicesE1),1)   2*ones(length(indicesE1),1)   indicesE1      indicesE1+Nx   indicesE1+2*Nx               2*ones(length(indicesE1),1)   indicesE1          indicesE1+Nx*Ny    indicesE1+2*Nx*Ny           ];
firstdevneighbours(indicesE2,:)   = [3*ones(length(indicesE2),1)   indicesE2      indicesE2-1   indicesE2-2                  1*ones(length(indicesE2),1)   indicesE2-Nx   indicesE2+Nx   zeros(length(indicesE2),1)   2*ones(length(indicesE2),1)   indicesE2          indicesE2+Nx*Ny    indicesE2+2*Nx*Ny           ];
firstdevneighbours(indicesE3,:)   = [1*ones(length(indicesE3),1)   indicesE3-1    indicesE3+1   zeros(length(indicesE3),1)   3*ones(length(indicesE3),1)   indicesE3      indicesE3-Nx   indicesE3-2*Nx               2*ones(length(indicesE3),1)   indicesE3          indicesE3+Nx*Ny    indicesE3+2*Nx*Ny           ];
firstdevneighbours(indicesE4,:)   = [2*ones(length(indicesE4),1)   indicesE4      indicesE4+1   indicesE4+2                  1*ones(length(indicesE4),1)   indicesE4-Nx   indicesE4+Nx   zeros(length(indicesE4),1)   2*ones(length(indicesE4),1)   indicesE4          indicesE4+Nx*Ny    indicesE4+2*Nx*Ny           ];
firstdevneighbours(indicesE5,:)   = [2*ones(length(indicesE5),1)   indicesE5      indicesE5+1   indicesE5+2                  2*ones(length(indicesE5),1)   indicesE5      indicesE5+Nx   indicesE5+2*Nx               1*ones(length(indicesE5),1)   indicesE5-Nx*Ny    indicesE5+Nx*Ny    zeros(length(indicesE5),1)  ];
firstdevneighbours(indicesE6,:)   = [3*ones(length(indicesE6),1)   indicesE6      indicesE6-1   indicesE6-2                  2*ones(length(indicesE6),1)   indicesE6      indicesE6+Nx   indicesE6+2*Nx               1*ones(length(indicesE6),1)   indicesE6-Nx*Ny    indicesE6+Nx*Ny    zeros(length(indicesE6),1)  ];
firstdevneighbours(indicesE7,:)   = [3*ones(length(indicesE7),1)   indicesE7      indicesE7-1   indicesE7-2                  3*ones(length(indicesE7),1)   indicesE7      indicesE7-Nx   indicesE7-2*Nx               1*ones(length(indicesE7),1)   indicesE7-Nx*Ny    indicesE7+Nx*Ny    zeros(length(indicesE7),1)  ];
firstdevneighbours(indicesE8,:)   = [2*ones(length(indicesE8),1)   indicesE8      indicesE8+1   indicesE8+2                  3*ones(length(indicesE8),1)   indicesE8      indicesE8-Nx   indicesE8-2*Nx               1*ones(length(indicesE8),1)   indicesE8-Nx*Ny    indicesE8+Nx*Ny    zeros(length(indicesE8),1)  ];
firstdevneighbours(indicesE9,:)   = [1*ones(length(indicesE9),1)   indicesE9-1    indicesE9+1   zeros(length(indicesE9),1)   2*ones(length(indicesE9),1)   indicesE9      indicesE9+Nx   indicesE9+2*Nx               3*ones(length(indicesE9),1)   indicesE9          indicesE9-Nx*Ny    indicesE9-2*Nx*Ny           ];
firstdevneighbours(indicesE10,:)  = [3*ones(length(indicesE10),1)  indicesE10     indicesE10-1  indicesE10-2                 1*ones(length(indicesE10),1)  indicesE10-Nx  indicesE10+Nx  zeros(length(indicesE10),1)  3*ones(length(indicesE10),1)  indicesE10         indicesE10-Nx*Ny   indicesE10-2*Nx*Ny          ];
firstdevneighbours(indicesE11,:)  = [1*ones(length(indicesE11),1)  indicesE11-1   indicesE11+1  zeros(length(indicesE11),1)  3*ones(length(indicesE11),1)  indicesE11     indicesE11-Nx  indicesE11-2*Nx              3*ones(length(indicesE11),1)  indicesE11         indicesE11-Nx*Ny   indicesE11-2*Nx*Ny          ];
firstdevneighbours(indicesE12,:)  = [2*ones(length(indicesE12),1)  indicesE12     indicesE12+1  indicesE12+2                 1*ones(length(indicesE12),1)  indicesE12-Nx  indicesE12+Nx  zeros(length(indicesE12),1)  3*ones(length(indicesE12),1)  indicesE12         indicesE12-Nx*Ny   indicesE12-2*Nx*Ny          ];

firstdevneighbours(indicesC1,:)   = [2*ones(length(indicesC1),1)   indicesC1      indicesC1+1   indicesC1+2                  2*ones(length(indicesC1),1)   indicesC1      indicesC1+Nx   indicesC1+2*Nx               2*ones(length(indicesC1),1)   indicesC1          indicesC1+Nx*Ny    indicesC1+2*Nx*Ny           ];
firstdevneighbours(indicesC2,:)   = [3*ones(length(indicesC2),1)   indicesC2      indicesC2-1   indicesC2-2                  2*ones(length(indicesC2),1)   indicesC2      indicesC2+Nx   indicesC2+2*Nx               2*ones(length(indicesC1),1)   indicesC2          indicesC2+Nx*Ny    indicesC2+2*Nx*Ny           ];
firstdevneighbours(indicesC3,:)   = [3*ones(length(indicesC3),1)   indicesC3      indicesC3-1   indicesC3-2                  3*ones(length(indicesC3),1)   indicesC3      indicesC3-Nx   indicesC3-2*Nx               2*ones(length(indicesC1),1)   indicesC3          indicesC3+Nx*Ny    indicesC3+2*Nx*Ny           ];
firstdevneighbours(indicesC4,:)   = [2*ones(length(indicesC4),1)   indicesC4      indicesC4+1   indicesC4+2                  3*ones(length(indicesC4),1)   indicesC4      indicesC4-Nx   indicesC4-2*Nx               2*ones(length(indicesC1),1)   indicesC4          indicesC4+Nx*Ny    indicesC4+2*Nx*Ny           ];
firstdevneighbours(indicesC5,:)   = [2*ones(length(indicesC5),1)   indicesC5      indicesC5+1   indicesC5+2                  2*ones(length(indicesC5),1)   indicesC5      indicesC5+Nx   indicesC5+2*Nx               3*ones(length(indicesC1),1)   indicesC5          indicesC5-Nx*Ny    indicesC5-2*Nx*Ny           ];
firstdevneighbours(indicesC6,:)   = [3*ones(length(indicesC6),1)   indicesC6      indicesC6-1   indicesC6-2                  2*ones(length(indicesC6),1)   indicesC6      indicesC6+Nx   indicesC6+2*Nx               3*ones(length(indicesC1),1)   indicesC6          indicesC6-Nx*Ny    indicesC6-2*Nx*Ny           ];
firstdevneighbours(indicesC7,:)   = [3*ones(length(indicesC7),1)   indicesC7      indicesC7-1   indicesC7-2                  3*ones(length(indicesC7),1)   indicesC7      indicesC7-Nx   indicesC7-2*Nx               3*ones(length(indicesC1),1)   indicesC7          indicesC7-Nx*Ny    indicesC7-2*Nx*Ny           ];
firstdevneighbours(indicesC8,:)   = [2*ones(length(indicesC8),1)   indicesC8      indicesC8+1   indicesC8+2                  3*ones(length(indicesC8),1)   indicesC8      indicesC8-Nx   indicesC8-2*Nx               3*ones(length(indicesC1),1)   indicesC8          indicesC8-Nx*Ny    indicesC8-2*Nx*Ny           ];

%对于cylinder来说，需要保证前后的边界重合，我们对此需要进行“缝合”--》firstdevneighbours,包括外边
if any(periodicity==1)
firstdevneighbours(indicesF1,5:8)   = [1*ones(length(indicesF2),1)   indicesF1+(Nz-1)*Nx*Ny    indicesF1+Nx*Ny          zeros(length(indicesF1),1)];
firstdevneighbours(indicesE1,5:8)   = [1*ones(length(indicesE1),1)   indicesE1+(Nz-1)*Nx*Ny    indicesE1+Nx*Ny          zeros(length(indicesE1),1)];
firstdevneighbours(indicesE2,5:8)   = [1*ones(length(indicesE2),1)   indicesE2+(Nz-1)*Nx*Ny    indicesE2+Nx*Ny          zeros(length(indicesE2),1)];
firstdevneighbours(indicesE3,5:8)   = [1*ones(length(indicesE3),1)   indicesE3+(Nz-1)*Nx*Ny    indicesE3+Nx*Ny          zeros(length(indicesE3),1)];
firstdevneighbours(indicesE4,5:8)   = [1*ones(length(indicesE4),1)   indicesE4+(Nz-1)*Nx*Ny    indicesE4+Nx*Ny          zeros(length(indicesE4),1)];
firstdevneighbours(indicesC1,5:8)   = [1*ones(length(indicesC1),1)   indicesC1+(Nz-1)*Nx*Ny    indicesC1+Nx*Ny          zeros(length(indicesC1),1)];
firstdevneighbours(indicesC2,5:8)   = [1*ones(length(indicesC2),1)   indicesC2+(Nz-1)*Nx*Ny    indicesC2+Nx*Ny          zeros(length(indicesC2),1)];
firstdevneighbours(indicesC3,5:8)   = [1*ones(length(indicesC3),1)   indicesC3+(Nz-1)*Nx*Ny    indicesC3+Nx*Ny          zeros(length(indicesC3),1)];
firstdevneighbours(indicesC4,5:8)   = [1*ones(length(indicesC4),1)   indicesC4+(Nz-1)*Nx*Ny    indicesC4+Nx*Ny          zeros(length(indicesC4),1)];

firstdevneighbours(indicesF6,5:8)   = [1*ones(length(indicesF6),1)   indicesF6-Nx*Ny           indicesF6-(Nz-1)*Nx*Ny   zeros(length(indicesF6),1)];
firstdevneighbours(indicesE9,5:8)   = [1*ones(length(indicesE9),1)   indicesE9-Nx*Ny           indicesE9-(Nz-1)*Nx*Ny   zeros(length(indicesE9),1)];
firstdevneighbours(indicesE10,5:8)  = [1*ones(length(indicesE10),1)  indicesE10-Nx*Ny          indicesE10-(Nz-1)*Nx*Ny  zeros(length(indicesE10),1)];
firstdevneighbours(indicesE11,5:8)  = [1*ones(length(indicesE11),1)  indicesE11-Nx*Ny          indicesE11-(Nz-1)*Nx*Ny  zeros(length(indicesE11),1)];
firstdevneighbours(indicesE12,5:8)  = [1*ones(length(indicesE12),1)  indicesE12-Nx*Ny          indicesE12-(Nz-1)*Nx*Ny  zeros(length(indicesE12),1)];
firstdevneighbours(indicesC5,5:8)   = [1*ones(length(indicesC5),1)   indicesC5-Nx*Ny           indicesC5-(Nz-1)*Nx*Ny   zeros(length(indicesC5),1)];
firstdevneighbours(indicesC6,5:8)   = [1*ones(length(indicesC6),1)   indicesC6-Nx*Ny           indicesC6-(Nz-1)*Nx*Ny   zeros(length(indicesC6),1)];
firstdevneighbours(indicesC7,5:8)   = [1*ones(length(indicesC7),1)   indicesC7-Nx*Ny           indicesC7-(Nz-1)*Nx*Ny   zeros(length(indicesC7),1)];
firstdevneighbours(indicesC8,5:8)   = [1*ones(length(indicesC8),1)   indicesC8-Nx*Ny           indicesC8-(Nz-1)*Nx*Ny   zeros(length(indicesC8),1)];

end
if any(periodicity==2)
firstdevneighbours(indicesF2,5:8)   = [1*ones(length(indicesF2),1)   indicesF2+(Ny-1)*Nx     indicesF2+Nx          zeros(length(indicesF2),1)];
firstdevneighbours(indicesE1,5:8)   = [1*ones(length(indicesE1),1)   indicesE1+(Ny-1)*Nx     indicesE1+Nx          zeros(length(indicesE1),1)];
firstdevneighbours(indicesE5,5:8)   = [1*ones(length(indicesE5),1)   indicesE5+(Ny-1)*Nx     indicesE5+Nx          zeros(length(indicesE5),1)];
firstdevneighbours(indicesE6,5:8)   = [1*ones(length(indicesE6),1)   indicesE6+(Ny-1)*Nx     indicesE6+Nx          zeros(length(indicesE6),1)];
firstdevneighbours(indicesE9,5:8)   = [1*ones(length(indicesE9),1)   indicesE9+(Ny-1)*Nx     indicesE9+Nx          zeros(length(indicesE9),1)];
firstdevneighbours(indicesC1,5:8)   = [1*ones(length(indicesC1),1)   indicesC1+(Ny-1)*Nx     indicesC1+Nx          zeros(length(indicesC1),1)];
firstdevneighbours(indicesC2,5:8)   = [1*ones(length(indicesC2),1)   indicesC2+(Ny-1)*Nx     indicesC2+Nx          zeros(length(indicesC2),1)];
firstdevneighbours(indicesC6,5:8)   = [1*ones(length(indicesC6),1)   indicesC6+(Ny-1)*Nx     indicesC6+Nx          zeros(length(indicesC6),1)];
firstdevneighbours(indicesC5,5:8)   = [1*ones(length(indicesC5),1)   indicesC5+(Ny-1)*Nx     indicesC5+Nx          zeros(length(indicesC5),1)];

firstdevneighbours(indicesF4,5:8)   = [1*ones(length(indicesF4),1)   indicesF4-Nx          indicesF4-(Ny-1)*Nx   zeros(length(indicesF4),1)];
firstdevneighbours(indicesE3,5:8)   = [1*ones(length(indicesE3),1)   indicesE3-Nx          indicesE3-(Ny-1)*Nx   zeros(length(indicesE3),1)];
firstdevneighbours(indicesE8,5:8)   = [1*ones(length(indicesE8),1)   indicesE8-Nx          indicesE8-(Ny-1)*Nx   zeros(length(indicesE8),1)];
firstdevneighbours(indicesE7,5:8)   = [1*ones(length(indicesE7),1)   indicesE7-Nx          indicesE7-(Ny-1)*Nx   zeros(length(indicesE7),1)];
firstdevneighbours(indicesE11,5:8)  = [1*ones(length(indicesE11),1)  indicesE11-Nx         indicesE11-(Ny-1)*Nx  zeros(length(indicesE11),1)];
firstdevneighbours(indicesC4,5:8)   = [1*ones(length(indicesC4),1)   indicesC4-Nx          indicesC4-(Ny-1)*Nx   zeros(length(indicesC4),1)];
firstdevneighbours(indicesC3,5:8)   = [1*ones(length(indicesC3),1)   indicesC3-Nx          indicesC3-(Ny-1)*Nx   zeros(length(indicesC3),1)];
firstdevneighbours(indicesC7,5:8)   = [1*ones(length(indicesC7),1)   indicesC7-Nx          indicesC7-(Ny-1)*Nx   zeros(length(indicesC7),1)];
firstdevneighbours(indicesC8,5:8)   = [1*ones(length(indicesC8),1)   indicesC8-Nx          indicesC8-(Ny-1)*Nx   zeros(length(indicesC8),1)];

end
if any(periodicity==3)
firstdevneighbours(indicesF3,5:8)   = [1*ones(length(indicesF3),1)   indicesF3-1            indicesF3-(Nx-1)          zeros(length(indicesF3),1)];
firstdevneighbours(indicesE4,5:8)   = [1*ones(length(indicesE4),1)   indicesE4-1            indicesE4-(Nx-1)          zeros(length(indicesE4),1)];
firstdevneighbours(indicesE5,5:8)   = [1*ones(length(indicesE5),1)   indicesE5-1            indicesE5-(Nx-1)          zeros(length(indicesE5),1)];
firstdevneighbours(indicesE8,5:8)   = [1*ones(length(indicesE8),1)   indicesE8-1            indicesE8-(Nx-1)          zeros(length(indicesE8),1)];
firstdevneighbours(indicesE12,5:8)  = [1*ones(length(indicesE12),1)  indicesE12-1           indicesE12-(Nx-1)         zeros(length(indicesE12),1)];
firstdevneighbours(indicesC4,5:8)  =  [1*ones(length(indicesC4),1)   indicesC4-1            indicesC4-(Nx-1)          zeros(length(indicesC4),1)];
firstdevneighbours(indicesC1,5:8)  =  [1*ones(length(indicesC1),1)   indicesC1-1            indicesC1-(Nx-1)          zeros(length(indicesC1),1)];
firstdevneighbours(indicesC5,5:8)  =  [1*ones(length(indicesC5),1)   indicesC5-1            indicesC5-(Nx-1)          zeros(length(indicesC5),1)];
firstdevneighbours(indicesC8,5:8)  =  [1*ones(length(indicesC8),1)   indicesC8-1            indicesC8-(Nx-1)          zeros(length(indicesC8),1)];

firstdevneighbours(indicesF5,5:8)   = [1*ones(length(indicesF5),1)   indicesF5+(Nx-1)       indicesF5+1   zeros(length(indicesF5),1)];
firstdevneighbours(indicesE2,5:8)   = [1*ones(length(indicesE2),1)   indicesE2+(Nx-1)       indicesE2+1   zeros(length(indicesE2),1)];
firstdevneighbours(indicesE6,5:8)   = [1*ones(length(indicesE6),1)   indicesE6+(Nx-1)       indicesE6+1   zeros(length(indicesE6),1)];
firstdevneighbours(indicesE7,5:8)   = [1*ones(length(indicesE7),1)   indicesE7+(Nx-1)       indicesE7+1   zeros(length(indicesE7),1)];
firstdevneighbours(indicesE10,5:8)  = [1*ones(length(indicesE10),1)  indicesE10+(Nx-1)      indicesE10+1  zeros(length(indicesE10),1)];
firstdevneighbours(indicesC3,5:8)   = [1*ones(length(indicesC3),1)   indicesC3+(Nx-1)       indicesC3+1   zeros(length(indicesC3),1)];
firstdevneighbours(indicesC2,5:8)   = [1*ones(length(indicesC2),1)   indicesC2+(Nx-1)       indicesC2+1   zeros(length(indicesC2),1)];
firstdevneighbours(indicesC6,5:8)   = [1*ones(length(indicesC6),1)   indicesC6+(Nx-1)       indicesC6+1   zeros(length(indicesC6),1)];
firstdevneighbours(indicesC7,5:8)   = [1*ones(length(indicesC7),1)   indicesC7+(Nx-1)       indicesC7+1   zeros(length(indicesC7),1)];

end



if flagintbounds
    for i=1:size(indicesintbounds,1)
        index = indicesintbounds(i,1);
        switch typeintbounds(i,1)
            case 1  % hole C1
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 2  % hole C2
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 3  % hole C3
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 4  % hole C4
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 5  % hole C5
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 6  % hole C6
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 7  % hole C7
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 8  % hole C8
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 9  % hole E1
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 10 % hole E2
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 11 % hole E3
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 12 % hole E4
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 13 % hole E5
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 14 % hole E6
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 15 % hole E7
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 16 % hole E8
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 17 % hole E9
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 18 % hole E10
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 19 % hole E11
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
            case 20 % hole E12
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0      1   index-Nx*Ny      index+Nx*Ny   0 ];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 21 % F1
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0           2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 22 % F2
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         2   index         index+Nx   index+2*Nx  1   index-Nx*Ny      index+Nx*Ny   0             ];
            case 23 % F3
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   1   index-Nx      index+Nx   0           1   index-Nx*Ny      index+Nx*Ny   0             ];
            case 24 % F4
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         3   index         index-Nx   index-2*Nx  1   index-Nx*Ny      index+Nx*Ny   0             ];
            case 25 % F5
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   1   index-Nx      index+Nx   0           1   index-Nx*Ny      index+Nx*Ny   0             ];
            case 26 % F6
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         1   index-Nx      index+Nx   0           3   index            index-Nx*Ny   index-2*Nx*Ny ];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 27 % C1
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   2   index         index+Nx   index+2*Nx  2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 28 % C2
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   2   index         index+Nx   index+2*Nx  2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 29 % C3
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   3   index         index-Nx   index-2*Nx  2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 30 % C4
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   3   index         index-Nx   index-2*Nx  2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 31 % C5
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   2   index         index+Nx   index+2*Nx  3   index            index-Nx*Ny   index-2*Nx*Ny ];
            case 32 % C6
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   2   index         index+Nx   index+2*Nx  3   index            index-Nx*Ny   index-2*Nx*Ny ];
            case 33 % C7
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   3   index         index-Nx   index-2*Nx  3   index            index-Nx*Ny   index-2*Nx*Ny ];
            case 34 % C8
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   3   index         index-Nx   index-2*Nx  3   index            index-Nx*Ny   index-2*Nx*Ny ];
                %---------------------------------------------------------------------------------------------------------------------------------------------------
            case 35 % E1
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         2   index         index+Nx   index+2*Nx  2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 36 % E2
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   1   index-Nx      index+Nx   0           2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 37 % E3
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         3   index         index-Nx   index-2*Nx  2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 38 % E4
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   1   index-Nx      index+Nx   0           2   index            index+Nx*Ny   index+2*Nx*Ny ];
            case 39 % E5
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   2   index         index+Nx   index+2*Nx  1   index-Nx*Ny      index+Nx*Ny   0             ];
            case 40 % E6
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   2   index         index+Nx   index+2*Nx  1   index-Nx*Ny      index+Nx*Ny   0             ];
            case 41 % E7
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   3   index         index-Nx   index-2*Nx  1   index-Nx*Ny      index+Nx*Ny   0             ];
            case 42 % E8
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   3   index         index-Nx   index-2*Nx  1   index-Nx*Ny      index+Nx*Ny   0             ];
            case 43 % E9
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         2   index         index+Nx   index+2*Nx  3   index            index-Nx*Ny   index-2*Nx*Ny ];
            case 44 % E10
                firstdevneighbours(index,:)   = [3   index      index-1   index-2   1   index-Nx      index+Nx   0           3   index            index-Nx*Ny   index-2*Nx*Ny ];
            case 45 % E11
                firstdevneighbours(index,:)   = [1   index-1    index+1   0         3   index         index-Nx   index-2*Nx  3   index            index-Nx*Ny   index-2*Nx*Ny ];
            case 46 % E12
                firstdevneighbours(index,:)   = [2   index      index+1   index+2   1   index-Nx      index+Nx   0           3   index            index-Nx*Ny   index-2*Nx*Ny ];
        end
    end
end

return