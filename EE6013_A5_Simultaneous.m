%% EE6013_A5.m

%% Clean Up

clear all;
close all;
clc

%% Variables

jacobian = [1,1,0];
% iterations = 50;
max = 1000;
multipleShootingIt = 20;
t = 1:1:max;

%%

a = zeros(max,3);
lambda = zeros(max,length(jacobian));
x1 = zeros(max,1);
x2 = zeros(max,1);
u = zeros(max,1);
u(1) = pi;
initialGuess = [1/4*pi;1/8*pi];
x1(1) = initialGuess(1);
x2(1) = initialGuess(2);

%% Sequential Shooting
count = 0;




%% Plotting

figure('name','X1')
plot(t,x1(:));

figure('name','X2')
plot(t,x2(:));


%% Function

function val = isPosDef(aPrev)

    hess = find_h(aPrev);
    h_eig = eig(hess);
    posEig = true;
    if length(hess) == 3
        for i = 1:length(h_eig)
            if (h_eig(i) <= 0)
                posEig = false;
            end
        end
    end
    
    if posEig
        val = true;
    else
        val = false;
    end
    
end


function a_i = find_g(aPrev,lambda)

    a_i = [cos(aPrev(1));-sin(aPrev(2));1] + lambda';

end


function hess = find_h(aPrev)

    hess = [-sin(aPrev(1)), 0,0; 0, -cos(aPrev(2)),0;0,0,0];

end
