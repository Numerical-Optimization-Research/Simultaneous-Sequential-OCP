%% EE6013_A5_P1.m


%% Clean Up

clear all;
close all;
clc 

warning('off','all')
warning
% id = 'MATLAB:nearlySingularMatrix';

%% Variables
eq = 9;

initialX = [-4; 1];
soln = [-pi/2; pi];

maxIteration = 35E3;
iterations = 1:maxIteration;

gradient = [2; 3];
hessian = [3, 0; 0, 4];

dt = zeros(1,eq);
val = zeros(1,eq);
error = zeros(length(initialX),maxIteration);
errorTot = zeros(eq,maxIteration);
outputX1 = zeros(eq,maxIteration);
outputX2 = zeros(eq,maxIteration);


%% Gradient Descent Fixed - 1
n = 1;
aNext = initialX;
outputX1(n,1) = aNext(1);
outputX2(n,1) = aNext(2);
lambdaFix = 0.002;

setting = "No";
timeStart = datetime('now');
for i = 1:maxIteration
    
    error(:,i) = soln(:) - aNext(:);
    aNext = gradientDescent_P2(aNext,gradient,setting,lambdaFix);
    outputX1(n,i+1) = aNext(1);
    outputX2(n,i+1) = aNext(2);
    
end
timeEnd = datetime('now');
dt_s = timeEnd - timeStart;
dt(n) = seconds(dt_s);
val(n) = dt(n)/maxIteration;
errorTot(n,:) = sqrt(error(1,:).^2 + error(2,:).^2);


%% Gradient Descent Line Search - 2
n = 2;
aNext = initialX;
outputX1(n,1) = aNext(1);
outputX2(n,1) = aNext(2);
lambdaFix = 1.5E-4;

setting = "Line";
timeStart = datetime('now');
for i = 1:maxIteration
    
    error(:,i) = soln(:) - aNext(:);
    aNext = gradientDescent_P2(aNext,gradient,setting,lambdaFix);
    outputX1(n,i+1) = aNext(1);
    outputX2(n,i+1) = aNext(2);
    
end
timeEnd = datetime('now');
dt_s = timeEnd - timeStart;
dt(n) = seconds(dt_s);
val(n) = dt(n)/maxIteration;
errorTot(n,:) = sqrt(error(1,:).^2 + error(2,:).^2);


%% Gradient Descent Barzilai-Borwein - 3




%% True Hessian - 4
n = 4;
aNext = initialX;
outputX1(n,1) = aNext(1);
outputX2(n,1) = aNext(2);

setting = "True";
timeStart = datetime('now');
for i = 1:maxIteration
    
    error(:,i) = soln(:) - aNext(:);
    aNext = hessianFunc_P2(aNext,gradient,setting,0,hessian);
    outputX1(n,i+1) = aNext(1);
    outputX2(n,i+1) = aNext(2);
        
end
timeEnd = datetime('now');
dt_s = timeEnd - timeStart;
dt(n) = seconds(dt_s);
val(n) = dt(n)/maxIteration;
errorTot(n,:) = sqrt(error(1,:).^2 + error(2,:).^2);


%% Gauss Newton - 5
n = 5;
aNext = initialX;
outputX1(n,1) = aNext(1);
outputX2(n,1) = aNext(2);

setting = "GN";
timeStart = datetime('now');
for i = 1:maxIteration
    
    error(:,i) = soln(:) - aNext(:);
    aNext = gaussNewtonFunc_P2(aNext,setting,gradient,0,0);
    outputX1(n,i+1) = aNext(1);
    outputX2(n,i+1) = aNext(2);
    
end
timeEnd = datetime('now');
dt_s = timeEnd - timeStart;
dt(n) = seconds(dt_s);
val(n) = dt(n)/maxIteration;
errorTot(n,:) = sqrt(error(1,:).^2 + error(2,:).^2);

%% Levenberg-Marquardt - 6
n = 6;
aNext = initialX;
outputX1(n,1) = aNext(1);
outputX2(n,1) = aNext(2);


lambdaLM = 200;

setting = "LM";
timeStart = datetime('now');
for i = 1:maxIteration
    
    error(:,i) = soln(:) - aNext(:);
    aNext = gaussNewtonFunc_P2(aNext,setting,gradient,lambdaLM,hessian);
    outputX1(n,i+1) = aNext(1);
    outputX2(n,i+1) = aNext(2);
    
end
timeEnd = datetime('now');
dt_s = timeEnd - timeStart;
dt(n) = seconds(dt_s);
val(n) = dt(n)/maxIteration;
errorTot(n,:) = sqrt(error(1,:).^2 + error(2,:).^2);

%% Constant Hessian - 7
n = 7;
lambdaHess = 80000;
constHess = lambdaHess*eye(length(initialX));
aNext = initialX;
outputX1(n,1) = aNext(1);
outputX2(n,1) = aNext(2);

setting = "Fixed";
timeStart = datetime('now');
for i = 1:maxIteration
    
    error(:,i) = soln(:) - aNext(:);
    aNext = hessianFunc_P2(aNext,gradient,setting,constHess,0);
    outputX1(n,i+1) = aNext(1);
    outputX2(n,i+1) = aNext(2);
    
end
timeEnd = datetime('now');
dt_s = timeEnd - timeStart;
dt(n) = seconds(dt_s);
val(n) = dt(n)/maxIteration;
errorTot(n,:) = sqrt(error(1,:).^2 + error(2,:).^2);


%% BFGS - 8
n = 8;
step = -pi/2;
Bstart = [4;1000];
B_i = Bstart.*eye(length(initialX));

aNext = initialX + step;
aCurr = initialX;
outputX1(n,1) = aCurr(1);
outputX2(n,1) = aCurr(2);
outputX1(n,2) = aNext(1);
outputX2(n,2) = aNext(2);

timeStart = datetime('now');
for i = 1:maxIteration
    
    aPrev = aCurr;
    aCurr = aNext;
    
    error(:,i) = soln(:) - aNext(:);
    aNext = BFGS_P2(aCurr,aPrev,gradient,B_i);
    outputX1(n,i+1) = aNext(1);
    outputX2(n,i+1) = aNext(2);
    
end
timeEnd = datetime('now');
dt_s = timeEnd - timeStart;
dt(n) = seconds(dt_s);
val(n) = dt(n)/maxIteration;
errorTot(n,:) = sqrt(error(1,:).^2 + error(2,:).^2);


%% Forced Hessian - 9
n = 9;
lambdaHess = 0.8;
constHess = lambdaHess*eye(length(initialX));
aNext = initialX;
outputX1(n,1) = aNext(1);
outputX2(n,1) = aNext(2);

setting = "Forced";
timeStart = datetime('now');
for i = 1:maxIteration
    
    error(:,i) = soln(:) - aNext(:);
    aNext = hessianFunc_P2(aNext,gradient,setting,constHess,hessian);
    outputX1(n,i+1) = aNext(1);
    outputX2(n,i+1) = aNext(2);
        
end
timeEnd = datetime('now');
dt_s = timeEnd - timeStart;
dt(n) = seconds(dt_s);
val(n) = dt(n)/maxIteration;
errorTot(n,:) = sqrt(error(1,:).^2 + error(2,:).^2);


%% 3D Plot 

[X1,X2] = meshgrid(-4*pi:0.1:4*pi,-4*pi:0.1:4*pi);

f = sin(X1) + cos(X2) + pi;

figure('name','Mesh Plot - P1 - Zoomed Out')
surf(X1,X2,f)
xlabel('X1');
ylabel('X2');
zlabel('F');

%% Contour Plot 
figure('name','Contour Plot - P1 - Zoomed Out')

subplot(3,1,1);
contour(X1,X2,f);
hold on 
plot(outputX1(1,:),outputX2(1,:));
for i = 2:eq
    if val(i) ~= 0
        plot(outputX1(i,:),outputX2(i,:));
    end
end
legend('Function Contour','1. Gradient Descent','2. Gradient Descent- Line Search','4. True Hessian','5. Gauss Newton','6. Levenberg-Marquardt','7. Constant Hessian','8. BFGS','9. Forced Hessian');
txt = 'Tuning Variables:';
% text(txt,lambda,)
xlabel('X1');
ylabel('X2');
xlim([-4*pi 4*pi]);
ylim([-pi/2 2*pi]);
hold off

subplot(3,1,2);
semilogx(errorTot(1,:),iterations);
hold on
for i = 2:eq
    if val(i) ~= 0
        semilogx(errorTot(i,:),iterations);
    end
end
xlabel('Iteration (semilogx)');
ylabel('Error');
xlim([1 maxIteration/20]);
% ylim([0 8000]);
hold off


subplot(3,1,3); 
hold on
for i = 1:eq
    if val(i) ~= 0
        timeError = 0:val(i):(dt(i)-val(i));
        plot(timeError,errorTot(i,:));
    end
end
xlabel('Time (s)');
ylabel('Error');
xlim([0 0.06]);
hold off


%% Display

disp('**********************************');
disp('Tuning Variables:');
disp('Gradient Descent fixed Lambda:  ');
disp(lambdaFix);
disp('Levenberg-Marquardt start Lambda:  ');
disp(lambdaLM);
disp('Constant Hessian start:  ');
disp(constHess);
disp('BFGS B_i start:  ');
disp(B_i);
disp('BFGS deltaX step:  ');
disp(step);

%% Edits
%
%
%
%
%