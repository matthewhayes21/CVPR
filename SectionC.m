% clear the work space and close all open plots
clear;
close all;

% set up
slash = '/'; %<------ use this to change all "/" to "\" or visa versa

colours = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#A2142F"];
names = ["Acrylic", "Black foam","Car sponge", "Flour sack", "Kitchen sponge","Steel vase"];

% define the colours ... because I'm a bit extra
orange = "#ff7f0e";
grey = "#4f4f4f";
black = "#000000";


% load the data
load('F0_PVT.mat')
load('F0_Electrodes.mat')
%% part a

% pull out the data for the 
Pressure = pressure(:,2:3);
Vibration = vibration(:,2:3);
Temprature = temprature(:,2:3);

%%
Ap = mean(Pressure, 'all');
Av = mean(Vibration, 'all');
At = mean(Temprature, 'all');

sigmap = std(Pressure, 0,'all');
sigmav = std(Vibration, 0,'all');
sigmat = std(Temprature, 0,'all');


Pressure = (Pressure - Ap)./sigmap;
Vibration = (Vibration - Av)./sigmav;
Temp = (Temprature - At)./sigmat;

%%
% calculate the group means
means = zeros(3,2);

for i = 1:2
    means(:,i) = [mean(Pressure(:, i));
                  mean(Vibration(:, i));
                  mean(Temp(:, i))];
end
%% Pressure vs Vibration

PV = [transpose(makeVector(Pressure)); transpose(makeVector(Vibration))];
mean = means(1:2,:);


% calculate Sw
Sw = real(zeros(2));

diff = PV - [kron(ones(1,10), mean(:, i)),kron(ones(1,10), mean(:, 2))];
Sw = Sw + diff*transpose(diff);


% calculate Sb
meanDiff = (mean(:, 1) - mean(:, 2));
Sb = meanDiff * transpose(meanDiff);
disp("---------------------------")
disp("For Pressure vs Vibration")
[eigenVectors, eigenValues] = eig(Sw\Sb)

% get the index of the colum with the max variance
[~,vecIndex] = max(max(eigenValues,[],1));

n = 0;
figure(1);
    hold on 
for i = 1:2
    scatter(Pressure(:,i),Vibration(:,i),'filled', "DisplayName", names(i+1),"MarkerFaceColor", colours(i+1));
    n = n+10;
    plot(mean(1,i),mean(2,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(i+1), "DisplayName",names(i+1)+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')
axis([floor(min(Pressure, [], 'all')) ceil(max(Pressure, [], 'all')) floor(min(Vibration, [], 'all')) ceil(max(Vibration, [], 'all'))])
xlabel('Pressure')
ylabel('Vibreaiton')
title("Pressure Vs Vibration for the 'Black foam' and 'Car sponge' with the LDA function")

%%
% Pressure vs Temp

PT = [transpose(makeVector(Pressure)); transpose(makeVector(Temp))];
mean = means(1:2:3,:);

% calculate Sw
Sw = real(zeros(2));

diff = PT - [kron(ones(1,10), mean(:, i)),kron(ones(1,10), mean(:, 2))];
Sw = Sw + diff*transpose(diff);

% calculate Sb
meanDiff = (mean(:, 1) - mean(:, 2));
Sb = meanDiff * transpose(meanDiff);
disp("---------------------------")
disp("For Pressure vs Temprature")
[eigenVectors, eigenValues] = eig(Sw\Sb)

% get the index of the colum with the max variance
[~,vecIndex] = max(max(eigenValues,[],1));

n = 0;
figure(2);
    hold on 
for i = 1:2
    scatter(Pressure(:,i),Temp(:,i),'filled', "DisplayName", names(i+1),"MarkerFaceColor", colours(i+1));
    n = n+10;
    plot(mean(1,i),mean(2,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(i+1), "DisplayName",names(i+1)+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')
axis([floor(min(Pressure, [], 'all')) ceil(max(Pressure, [], 'all')) ...
    floor(min(Temp, [], 'all')) ceil(max(Temp, [], 'all'))])
xlabel('Pressure')
ylabel('Temprature')
title("Pressure Vs Temprature for the 'Black foam' and 'Car sponge' with the LDA function")


%%
% Temp vs Vibration
TV = [transpose(makeVector(Temp)); transpose(makeVector(Vibration))];
mean = means(3:-1:2,:);

% calculate Sw
Sw = real(zeros(2));

diff = TV - [kron(ones(1,10), mean(:, i)),kron(ones(1,10), mean(:, 2))];
Sw = Sw + diff*transpose(diff);

% calculate Sb
meanDiff = (mean(:, 1) - mean(:, 2));
Sb = meanDiff * transpose(meanDiff);
disp("---------------------------")
disp("For Temprature vs Vibrations ")
[eigenVectors, eigenValues] = eig(Sw\Sb)

% get the index of the colum with the max variance
[~,vecIndex] = max(max(eigenValues,[],1));

n = 0;
figure(3);
    hold on 
for i = 1:2
    scatter(Temp(:,i),Vibration(:,i),'filled', "DisplayName", names(i+1),"MarkerFaceColor", colours(i+1));
    n = n+10;
    plot(mean(1,i),mean(2,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(i+1), "DisplayName",names(i+1)+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')
axis([floor(min(Temp, [], 'all')) ceil(max(Temp, [], 'all'))...
    floor(min(Vibration, [], 'all')) ceil(max(Vibration, [], 'all'))])
xlabel('Temprature')
ylabel('Vibrations')
title("Temprature Vs Vibrations for the 'Black foam' and 'Car sponge' with the LDA function")
