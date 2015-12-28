function  [W,Z] = gaussquad(element,nint)

%Tetrahedra integration rules
if strcmpi(element,'tet')
    switch nint
        case 1
            W = 1;
            Z = [0.25,0.25,0.25,0.25];
        case 4
            W = 0.25*ones(4,1);
            r1=(5+3*sqrt(5))/20;
            s1=(5-sqrt(5))/20;
            t1=s1;
            u1=1-r1-s1-t1;
            
            Z = [r1,s1,t1,u1;...
                 s1,t1,u1,r1;...
                 t1,u1,r1,s1;...
                 u1,r1,s1,t1];
        case 5
            W=[-0.8;0.45*ones(4,1)];
            Z=[0.25,0.25,0.25,0.25;...
                0.5, 1/6, 1/6, 1/6;...
                1/6, 0.5, 1/6, 1/6;...
                1/6, 1/6, 0.5, 1/6;
                1/6, 1/6, 1/6, 0.5];
        case 11
            fprintf('Integration order not implemented yet');
        case 15
            fprintf('Integration order not implemented yet');
    end
end
end