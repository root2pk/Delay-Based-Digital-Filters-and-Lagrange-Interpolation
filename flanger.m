%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                     PROGRAM : FLANGER
%
%   PROGRAM TO TAKE IN AN INPUT SIGNAL AND PRODUCE A FLANGING EFFECT ON IT
%   USING AN LFO WHO'S PARAMETERS ARE USER DEFINED
%
%   THE PROGRAM USES A SPECIFIC ALGORITHM TO TO INTERPOLATE FOR N=1,2 AND 4,
%   AND A GENERAL ALGORITHM FOR N>4
%
%   THE VALUE OF N CAN BE CHANGED FROM WITHIN THE SCRIPT BY THE USER
%
%
%   AUTHOR : RUTHU PREM KUMAR
%   DATE : 06/12/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all;

%read in our WAV file, and store sample rate in Fs
[x,Fs]=audioread('Guitar_dry.wav');

%in case of stereo, to mono
x = 0.5*sum(x,2);
L = length(x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Number of nearest neighbours
N = 2;

fmode = 1;     %Mode of operation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LFO Parameters


MO = 0.002;              %MO = 2 ms
g = 0.5;                 %g = 0.5 Effect strength
fo = 2;                  %fo = 2 Hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Converting MO in terms of sample rate

MO = round(MO*Fs);


%Starting position for reading x values to avoid reading negative indices
start_pos = 2*MO;

%Zero padding the input signal to prevent reading negative indices
x = vertcat(zeros(start_pos,1),x);

%Recalculating the length of signal x%Zero padding the x signal to prevent reading negative values
L_new = length(x);

%Creating the output vector y
y = zeros(L_new,1);

%Creating LFO M
n = [0:L_new];
M(n+1) = MO*(1+sin(2*pi*fo*n/Fs));

%Converting M to column vector
M=M.';

%Plotting LFO

f1 = figure;

subplot (3,1,1)
plot ([1:length(M)]/Fs,M);
xlabel('Time(s)'); ylabel ('LFO Delay time(samples)'); title('LFO Plot'); ylim([-MO 3*MO]);

 % ----------------------INTERPOLATION--------------------------------
%Zeroth order interpolation
if N==1

    %Loop to calculate output vector y
    for n=start_pos:L_new

        y(n) = x(n)+g*x(n-round(M(n)));

    end

%Linear Interpolation
elseif N==2

    %Loop to calculate output vector y
    for n=start_pos:L_new

        M_int = floor(M(n));
        alpha = M(n) - M_int - 0.5;
        x_minus = x(n-floor(M(n)));
        x_plus = x(n-(floor(M(n))+1));

        P = (alpha+0.5)*x_plus-(alpha - 0.5)*x_minus;

        y(n) = x(n) + g*P;
    end

%Cubic Interpolation
elseif N==4

    %Loop to calculate output vector y
    for n=start_pos:L_new

        M_int = floor(M(n));                       %Calculating M_int
        alpha = M(n) - M_int - 0.5;                %calculating alpha
        M_intarray = [-(N-2)/2:N/2];               %Vector storing values of x positions for N neighbours
        l = [1:N];
        X(l) = x(n-(M_int+M_intarray(l)));         %Vector storing values of input signal at x positions of N neighbours


        %Calculating the coefficents
        P(1) = (alpha+0.5)*(alpha-0.5)*(alpha-1.5)/-6;
        P(2) = (alpha+1.5)*(alpha-0.5)*(alpha-1.5)/2;
        P(3) = (alpha+0.5)*(alpha-1.5)*(alpha+1.5)/-2;
        P(4) = (alpha+0.5)*(alpha-0.5)*(alpha+1.5)/6;

        %Calculating the polynomial P
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(l);
        end

        %Assigning to output vector
        y(n) = x(n) + g*Ptotal;
    end

%Interpolation above N = 4
elseif N>4

    Q = input ('Enter value of Q');
    P = linear_interp(N,Q,fmode);

    for n=start_pos:L_new

        %Calculating M_int
        M_int = floor(M(n));

        %calculating alpha
        alpha = M(n) - M_int - 0.5;
        q = [1:Q];
        aq = (-Q/2+q-1)/Q;
        [val,idx]=min(abs(aq-alpha));
        M_intarray = [-(N-2)/2:N/2];               %Vector storing values of x positions for N neighbours
        l = [1:N];
        X(l) = x(n-(M_int+M_intarray(l)));         %Vector storing values of input signal at x positions of N neighbours
        Ptotal = 0;
        for l =1:N
            Ptotal = Ptotal + X(l)*P(idx,l);
        end
        %Assigning to output vector
        y(n) = x(n) + g*Ptotal;
    end
else
    disp('Error');
end


%%Plot the input and output vectors only if fmode = 1
%%Else plot only the polynomials


if fmode == 1
    figure(f1);
    %Graph plot for input signal
    subplot(3,1,2);
    plot([1:length(x)]/Fs,x,'r');
    xlabel('Time(s)'); ylabel('Amplitude'); title({'';'';'Input Signal'});

    %Output signal
    subplot(3,1,3);
    plot([1:length(y)]/Fs,y,'b');
    xlabel('Time(s)'); ylabel('Amplitude'); title({'Output Signal'});

    %Playing the output audio file
    soundsc(y,Fs);
end

