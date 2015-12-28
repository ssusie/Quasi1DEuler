function [Q1,R,Q2]=givens_rotations1(R)
    Q1=eye(length(R));
    Q2=eye(length(R));
    
    %% Main iterations
    for i=1:length(R)-3
        % Left Multiplication
        [G1,y1] = planerot(flipud(R(i:i+1,length(R))));
        G1=G1';
        Q1_1=eye(length(R));
        Q1_1(i:i+1,i:i+1)=G1;
        Q1=Q1_1*Q1;
        R(i:i+1,:)=G1*R(i:i+1,:);
%         R=Q1_1*R;
        
        % Right Multiplication
        [G2,y2]=planerot(fliplr(R(i+1,i:i+1))');
        Q2_1=eye(length(R));
        Q2_1(i:i+1,i:i+1)=G2;
        Q2=Q2*Q2_1;
        R(:,i:i+1)=R(:,i:i+1)*G2;
%         R=R*Q2_1;
        
        % Correction: To eleminate 4 elements
        [Q1,R,Q2]= correction(Q1,R,Q2,i);
    end
    %% Correction: To eleminate 4 elements
    function [Q1,R,Q2]= correction(Q1,R,Q2,i)
        % Determine number of correction iterations
        if (mod(i,2)==0)
            j=i/2-1;
        else
            j=(i-1)/2;
        end
        
        % Carry out corrections
        for k=1:j
            [G1,y1] = planerot(flipud(R(i-2:i-1,i+1)));
            G1=G1';
            Q1_1=eye(length(R));
            Q1_1(i-2:i-1,i-2:i-1)=G1;
%             R=Q1_1*R;
            R(i-2:i-1,:)=G1*R(i-2:i-1,:);
            Q1=Q1_1*Q1;
            
            [G2,y2]=planerot(fliplr(R(i-1,i-2:i-1))');
            Q2_1=eye(length(R));
            Q2_1(i-2:i-1,i-2:i-1)=G2;
%             R=R*Q2_1;
            R(:,i-2:i-1)=R(:,i-2:i-1)*G2;
            Q2=Q2*Q2_1;
            
            i=i-2; % Move 2 rows and columns up
        end
    end
end