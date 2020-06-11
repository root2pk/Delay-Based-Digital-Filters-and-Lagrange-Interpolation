%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    FUNCTION : GENERAL LAGRANGE INTERPOLATION
%    
%    THIS FUNCTION TAKES IN ARGUMENTS N,Q AND FMODE.
%    
%    FMODE = 1 : COMPUTES A LOOK UP TABLE FOR VALUES OF ALPHA[-0.5 0.5) AND 
%    THE CORRESPONDING COEFFICIENT OF THE LAGRANGE POLYNOMIAL FOR N NEIGHBOURS.
%    
%    FMODE = 2 : COMPUTES A LOOK UP TABLE FOR VALUES OF 
%    ALPHA(-(N-1)/2,(N-1)/2) AND THE CORRESPONDING COEFFICIENT OF THE 
%    LAGRANGE POLYNOMIAL FOR N NEIGHBOURS AND PLOTS IT.
%    
%    AUTHOR : RUTHU PREM KUMAR
%    DATE :06/12/2019
%    
%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  



function [P] = MA2_s1760559_Premkumar_Linterp(N,Q,fmode)

    %If N isnt divisible by 2, then exit function and return void aray.
    if mod(N,2)~=0 || N<1
        disp ('Error, N must be even and greater than 1');
        P=[];

    %%%%%%%%%%%%%%%%%%%%%%%   FMODE 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif fmode==1

        q = [1:Q];
        aq = (-Q/2+q-1)/Q;            %Vector of alpha values [-0.5,0.5)
        list1 = [-(N-1)/2:(N-1)/2];   %Vector of values from -(N-1)/2 to (N-1)/2

        %For loop to calculate the array P

        for n = 1:N
            P(:,n) = ones(Q,1);
            for m =1:Q
                k = list1(list1~=list1(n));
                P(m,n) = prod((aq(m) - k)./(list1(n)-k)); %Coefficient
            end
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%   FMODE = 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif fmode==2

        aq = linspace(-(N-1)/2,(N-1)/2,Q);        %Vector of alpha values
        list1 = [-(N-1)/2:(N-1)/2];               %Vector of values from -(N-1)/2 to (N-1)/2

        %For loop to calculate the array P

        for n = 1:N
           P(:,n) = ones(Q,1);
           for m =1:Q
                k = list1(list1~=list1(n));
                P(m,n) = prod((aq(m) - k)./(list1(n)-k));
           end
        end

        %%%%%%%%%       POLYNOMIAL PLOT    %%%%%%%%%%%%%%%%%%%%%%%%%
        figure;
        plot(aq,P,'Linewidth',1);  %Plot graph

        %legend and axes limits
        legend(num2str((-(N-1)/2:(N-1)/2)'));
        ylim([-1 2]); xlim([-(N-1)/2 (N-1)/2]);

    end

end