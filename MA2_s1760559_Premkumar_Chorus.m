% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%      PROGRAM : CHORUS EFFECT
%      
%  THIS PROGRAM TAKES IN AN IPUT SIGNAL,
%  CONVERTS IT TO MONO AND THEN USING TWO LFOS PER CHANNEL,
%  PERFORMS A CHORUS EFFECT ON IT AND PRODUCES AN OUTPUT SIGNAL.
%  THE LFOS ARE RANDOMLY GENERATED SUM OF SINUSOIDS WITH USER DEFINED DEPTHS
%  AND MAXIMUM VALUES OF DELAY. 
%  THE NUMBER OF SINUSOIDS TO SUM OVER IS DEFINED BY THE USER IN THE 
%  BEGINNING(PARAMETER N_rand) 
%  THE PROGRAM USES A SPECIFIC ALGORITHM TO TO INTERPOLATE FOR N=1,2 AND 4,
%  AND A GENERAL ALGORITHM FOR N>4
%  
%  THE VALUE OF N CAN  BE CHANGED FROM WITHIN THE SCRIPT, AND THE VALUE OF 
%  Q IS OBTAINED FROM THE USER WHILE THE PROGRAM IS RUNNING
% 
% 
%     AUTHOR : RUTHU PREM KUMAR
%     DATE : 06/12/2019
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 






clear all; close all;

%read in our WAV file, and store sample rate in Fs
[x,Fs]=audioread('Guitar_dry.wav');

%in case of stereo, to mono
x = 0.5*sum(x,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Setting the value of number of nearest neighbours
N = 4;

%Setting the value of fmode
fmode = 1;

%Number of sinusoids to sum over to create random LFOs(User defined)
N_rand = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LFO Parameters for both the LFOs

%--------------LEFT CHANNEL------------------------
%LFO 1
MO1 = 0.010;              %MO1 = 10 ms
D1 = 0.003;               %D1 = 2 ms
f1 = 0.5;                 %f1 = 0.5 Hz
g1 = 0.1;                 %g1 = 0.1

%LFO 2
MO2 = 0.015;              %MO2 = 15 ms
D2 = 0.005;               %D2 = 5 ms
f2 = 1;                   %f2 = 1 Hz
g2 = 0.25;                %g2 = 0.25

%----------------------------------------------------
%---------------RIGHT CHANNEL------------------------

%LFO 3
MO3 = 0.020;              %MO3 = 20 ms
D3 = 0.007;               %D3 = 7 ms
f3 = 1.5;                 %f3 = 1.5 Hz
g3 = -0.25;               %g3 = -0.25

%LFO 4
MO4 = 0.035;              %MO4 = 35 ms
D4 = 0.010;               %D4 = 10 ms
f4 = 2;                   %f4 = 2 Hz
g4 = -0.5;                %g4 = -0.5

%---------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Converting MO and D values in terms of samples
MO1 = round(MO1*Fs); MO2 = round(MO2*Fs);  MO3 = round(MO3*Fs); 
MO4 = round(MO4*Fs);
D1 = round(D1*Fs); D2 = round(D2*Fs); D3 = round(D3*Fs); D4 = round(D4*Fs);

%Checking which LFO has highest delay and zero padding accordingly
start_pos = max([MO1, MO2, MO3, MO4]);
%Zero padding the input signal to avoid negative indices
x = vertcat(zeros(start_pos,1),x);

%Calculating new length of input signal
L_new = length(x);

%Creating output signals
yR = zeros(L_new,1);   %Right
yL = yR;               %Left      
y = zeros(L_new,2);    %Total


%---------------Creating LFOs-----------------------------------------

%Sinusoids with the frequencies f1,f2,f3,f4
n = [0:L_new];
s(1,:) = sin(2*pi*f1*n/Fs);
s(2,:) = sin(2*pi*f2*n/Fs);
s(3,:) = sin(2*pi*f3*n/Fs);
s(4,:) = sin(2*pi*f4*n/Fs);


%Modifying the sinusoids by adding N_rand sinusoids of random frequencies
for t =1:4
    s(t+4,:) = zeros(1,L_new+1);
    for k =1:N_rand
        %Vector with random values of frequency in range (0.2 Hz,2 HZ)
        f_rand(k) = 0.2 + rand(1)*1.8;       
        %Adding random sinusoids to s
        s(t+4,:) = s(t+4,:)+sin(2*pi*f_rand(k)*n/Fs); 
    end
    %Normalizing the sinusoids obtained
    s(t+4,:) = (s(t+4,:)/max(abs(s(t+4,:))))*max(abs(s(t,:)));
                
end

%Using the sinusoid combinations to create the LFOs

%Creating LFO M1
M1(n+1) = MO1 + D1*(s(5,:)-1);
M1=M1.';

%Creating LFO M2

M2(n+1) = MO2 + D2*(s(6,:)-1);
M2=M2.';
 
%Creating LFO M3
M3(n+1) = MO3 + D3*(s(7,:)-1);
M3=M3.';

%Creating LFO M4
M4(n+1) = MO4 + D4*(s(8,:)-1);
M4=M4.';


%---------------------PLOTTING LFOS-------------------------------
%Plotting LFO M1

f1 = figure; 

subplot(2,3,1);

plot([1:length(M1)]/Fs, M1,'r');
ylim([0 MO1+100]);
ylabel('LFO delay time(samples)'); xlabel('Time(s)'); title('M1(Left)');

%Plotting LFO M2
subplot(2,3,2);

plot([1:length(M2)]/Fs,M2,'r');
ylim([0 MO2+100]);
ylabel('LFO delay time(samples)'); xlabel('Time(s)'); title('M2(Left)');


%Plotting LFO M3
subplot(2,3,4);

plot([1:length(M3)]/Fs, M3,'g');
ylim([0 MO3+100]);
ylabel('LFO delay time(samples)'); xlabel('Time(s)'); title('M3(Right)');


%Plotting LFO M4
subplot(2,3,5);

plot([1:length(M4)]/Fs, M4,'g');
ylim([0 MO4+100]);
ylabel('LFO delay time(samples)'); xlabel('Time(s)'); title('M4(Left)');


%--------------------------INTERPOLATION-----------------------------
%Zeroth order interpolation
if N==1

    %Loop to calculate output vector y
    for n=start_pos:L_new
        %Left channel M1 and M2
        yR(n) = x(n) + g1*x(n - round(M1(n))) + g2*x(n - round(M2(n)));
        %Right channel M3 and M4
        yL(n) = x(n) + g3*x(n - round(M3(n))) + g4*x(n - round(M4(n)));
        
        
    end

 
%Linear Interpolation    
elseif N==2
    
    %Loop to calculate output vector y
    for n=start_pos:L_new
        
        %-----------------M1----------------------------------
        M_int = floor(M1(n));
        alpha = M1(n) - M_int - 0.5;
        x_minus = x(n-floor(M1(n)));
        x_plus = x(n-(floor(M1(n))+1));
        
        Ptotal = (alpha+0.5)*x_plus-(alpha - 0.5)*x_minus;
        
        yL(n) = x(n) + g1*Ptotal;
        
        %-----------------M2----------------------------------
        M_int = floor(M2(n));
        alpha = M2(n) - M_int - 0.5;
        x_minus = x(n-floor(M2(n)));
        x_plus = x(n-(floor(M2(n))+1));
        
        Ptotal = (alpha+0.5)*x_plus-(alpha - 0.5)*x_minus;
        
        yL(n) = yL(n) + g2*Ptotal;
        
        %-----------------M3----------------------------------
        M_int = floor(M3(n));
        alpha = M3(n) - M_int - 0.5;
        x_minus = x(n-floor(M3(n)));
        x_plus = x(n-(floor(M3(n))+1));
        
        Ptotal = (alpha+0.5)*x_plus-(alpha - 0.5)*x_minus;
        
        yR(n) = x(n) + g3*Ptotal;
        
        %-----------------M4----------------------------------
        M_int = floor(M4(n));
        alpha = M4(n) - M_int - 0.5;
        x_minus = x(n-floor(M4(n)));
        x_plus = x(n-(floor(M4(n))+1));
        
        Ptotal = (alpha+0.5)*x_plus-(alpha - 0.5)*x_minus;
        
        yR(n) = yR(n) + g4*Ptotal;
        
    end

%Cubic Interpolation    
elseif N==4
    
    %Loop to calculate output vector y
    for n=start_pos:L_new
        
        %---------------------------M1-------------------------------
        M_int = floor(M1(n));
        alpha = M1(n) - M_int - 0.5;
        
        M_intarray = [-(N-2)/2:N/2];               %Vector storing values of x positions for N neighbours
        l = [1:N];
        X(l) = x(n-(M_int+M_intarray(l)));         %Vector storing values of input signal at x positions of N neighbours         
        
        P(1) = (alpha+0.5)*(alpha-0.5)*(alpha-1.5)/-6;
        P(2) = (alpha+1.5)*(alpha-0.5)*(alpha-1.5)/2;
        P(3) = (alpha+0.5)*(alpha-1.5)*(alpha+1.5)/-2;
        P(4) = (alpha+0.5)*(alpha-0.5)*(alpha+1.5)/6;
        
        %Calculating the polynomial Ptotal
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(l);
        end
        
        yL(n) = x(n) + g1*Ptotal;
        
        
        %----------------------------M2-----------------------------
        M_int = floor(M2(n));
        alpha = M2(n) - M_int - 0.5;
        
      
        X(l) = x(n-(M_int+M_intarray(l)));         %Vector storing values of input signal at x positions of N neighbours         
        
        P(1) = (alpha+0.5)*(alpha-0.5)*(alpha-1.5)/-6;
        P(2) = (alpha+1.5)*(alpha-0.5)*(alpha-1.5)/2;
        P(3) = (alpha+0.5)*(alpha-1.5)*(alpha+1.5)/-2;
        P(4) = (alpha+0.5)*(alpha-0.5)*(alpha+1.5)/6;
        
        %Calculating the polynomial Ptotal
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(l);
        end
              
        yL(n) = yL(n) + g2*Ptotal;
        
         
        %-----------------------------M3--------------------------------
        M_int = floor(M3(n));
        alpha = M3(n) - M_int - 0.5;
        
      
        X(l) = x(n-(M_int+M_intarray(l)));         %Vector storing values of input signal at x positions of N neighbours         
        
        P(1) = (alpha+0.5)*(alpha-0.5)*(alpha-1.5)/-6;
        P(2) = (alpha+1.5)*(alpha-0.5)*(alpha-1.5)/2;
        P(3) = (alpha+0.5)*(alpha-1.5)*(alpha+1.5)/-2;
        P(4) = (alpha+0.5)*(alpha-0.5)*(alpha+1.5)/6;
        
        %Calculating the polynomial Ptotal
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(l);
        end
              
        yL(n) = x(n) + g3*Ptotal;
        
         
        %------------------------------M4---------------------------------
        M_int = floor(M4(n));
        alpha = M4(n) - M_int - 0.5;
        
      
        X(l) = x(n-(M_int+M_intarray(l)));         %Vector storing values of input signal at x positions of N neighbours         
        
        P(1) = (alpha+0.5)*(alpha-0.5)*(alpha-1.5)/-6;
        P(2) = (alpha+1.5)*(alpha-0.5)*(alpha-1.5)/2;
        P(3) = (alpha+0.5)*(alpha-1.5)*(alpha+1.5)/-2;
        P(4) = (alpha+0.5)*(alpha-0.5)*(alpha+1.5)/6;
        
        %Calculating the polynomial Ptotal
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(l);
        end
              
        yR(n) = yR(n) + g4*Ptotal;
        
    end
    

   
elseif N>4
    
    %Accept value of Q from user
    Q = input ('Enter value of Q');
    
    %Call function to calculate look up table
    P = MA2_s1760559_Premkumar_Linterp(N,Q,fmode);
    
    for n=start_pos:L_new
        
        q = [1:Q];
        aq = (-Q/2+q-1)/Q;
        M_intarray = [-(N-2)/2:N/2];               %Vector storing values of x positions for N neighbours
        l = [1:N];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%   M1    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculating M_int
        M_int = floor(M1(n));    
        
        %calculating alpha
        alpha = M1(n) - M_int - 0.5;                  
        [~,idx]=min(abs(aq-alpha));   
        
        %Vector storing values of input signal at x positions of N neighbours                    
        X(l) = x(n-(M_int+M_intarray(l)));      
        
        %Calculating Polynomial Ptotal
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(idx,l);
        end   
        
        %Assigning to output vector
        yL(n) = x(n) + g1*Ptotal;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%    M2    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        %Calculating M_int
        M_int = floor(M2(n));    
        
        %calculating alpha
        alpha = M2(n) - M_int - 0.5;                  
        [val,idx]=min(abs(aq-alpha));   
        
        %Vector storing values of input signal at x positions of N neighbours                    
        X(l) = x(n-(M_int+M_intarray(l)));      
        
        %Calculating Polynomial Ptotal
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(idx,l);
        end   
        
        %Assigning to output vector
        yL(n) = yL(n) + g2*Ptotal;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%   M3    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculating M_int
        M_int = floor(M3(n));    
        
        %calculating alpha
        alpha = M3(n) - M_int - 0.5;                  
        [~,idx]=min(abs(aq-alpha));   
        
        %Vector storing values of input signal at x positions of N neighbours                    
        X(l) = x(n-(M_int+M_intarray(l)));      
        
        %Calculating Polynomial Ptotal
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(idx,l);
        end   
        
        %Assigning to output vector
        yR(n) = x(n) + g3*Ptotal;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%   M4    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculating M_int
        M_int = floor(M4(n));    
        
        %calculating alpha
        alpha = M4(n) - M_int - 0.5;                  
        [~,idx]=min(abs(aq-alpha));   
        
        %Vector storing values of input signal at x positions of N neighbours                    
        X(l) = x(n-(M_int+M_intarray(l)));      
        
        %Calculating Polynomial Ptotal
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(idx,l);
        end   
        
        %Assigning to output vector
        yR(n) = yR(n) + g4*Ptotal;
        
    end
    
    
    
%else display Error    
else
    disp('Error');
end

%%Display input and output signal only if fmode=1, else only display
%%polynomials

if fmode == 1
    
   
    figure(f1);
    
    %Assigning channels to output array y
    y(:,1) = yL; y(:,2) = yR;
    
    %Graph plot for input signal
           
    subplot(2,3,3);
    plot([1:length(x)]/Fs,x,'r');
    xlabel('Time(s)'); ylabel('Amplitude'); title('Input Signal');
    
    %plot for Output signal
    subplot(2,3,6);
    plot([1:length(y)]/Fs,y,'b');
    xlabel('Time(s)'); ylabel('Amplitude'); title('Output Signal');
    
    
    
    %Playing the output audio file
    soundsc(y,Fs);
    
end


