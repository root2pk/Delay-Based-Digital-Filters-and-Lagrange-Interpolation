%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  COMB FILTERING
% 
%    PROGRAM TO PERFORM FEEDBACK AND FEEDFORWARD
%     ON A MONO SIGNAL AND LISTEN TO THE OUTPUT
%
%    AUTHOR : RUTHU PREM KUMAR
%    DATE : 06/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



clear all; close all;


%read in our WAV file, and store sample rate in Fs
[x,Fs]=audioread('Guitar_dry.wav');

%in case of stereo, to mono
x = 0.5*sum(x,2);

%Set up comb filtering parameters
feed_time = 0.1;            %Feed time in seconds
M = Fs*feed_time;           %Feed time converted in terms of samples
g = 0.5;                    %Strength of effect

%Zero padding to prevent negaive indices
x=vertcat(zeros(M,1),x);

%creating empty vectors to store output signals
L = length(x);
y_ff = zeros(L,1);
y_fb = zeros(L,1);

%Loop to creat eoutput vectors
for n=M+1:L
    y_ff(n) = x(n)+g*x(n-M);
    y_fb(n) = x(n)-g*y_fb(n-M);
end

%Plotting the feedback and feedfoward signals
figure;
%------------INPUT SIGNAL--------------------
subplot(3,1,1);
plot([1:length(x)]/Fs,x,'r');
xlabel('Time(s)'); ylabel('Amplitude'); title('Input Signal'); 
ylim([min(x) max(x)]);

%---------------FEEDFORWARD SIGNAL----------------
subplot(3,1,2);
plot([1:length(y_ff)]/Fs,y_ff,'b');
xlabel('Time(s)'); ylabel('Amplitude'); title('Feedforward Signal');
ylim([min(y_ff) max(y_fb)]);

%---------------FEEDBACK SIGNAL---------------------
subplot(3,1,3);
plot([1:length(y_fb)]/Fs,y_fb,'g');
xlabel('Time(s)'); ylabel('Amplitude'); title('Feedback Signal');
ylim([min(y_fb) max(y_fb)]);
%Playing output files (Comment out unwanted signal)

soundsc(y_fb,Fs);
% soundsc(y_ff,Fs);