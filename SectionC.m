% clear the work space and close all open plots
clear;
close all;

% set up
slash = '/'; %<------ use this to change all "/" to "\" or visa versa

colours = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#A2142F"];
names = ["Acrylic", "Black foam","Car sponge", "Flour sack", "Kitchen sponge", "Steel vase"];

% define the colours ... because I'm a bit extra
orange = "#ff7f0e";
grey = "#4f4f4f";
black = "#000000";

% load the data
load('F0_PVT.mat')
load('F0_Electrodes.mat')

clc;
%% part a
disp("----------------------------")
disp("Part A")
disp("----------------------------")
% pull out the data for the balck foam and the car sponge
Pressure = pressure(:,2:3);
Vibration = vibration(:,2:3);
Temprature = temprature(:,2:3);


% standadise the data
Ap = mean(Pressure, 'all');
Av = mean(Vibration, 'all');
At = mean(Temprature, 'all');

sigmap = std(Pressure, 0,'all');
sigmav = std(Vibration, 0,'all');
sigmat = std(Temprature, 0,'all');


Pressure = (Pressure - Ap)./sigmap;
Vibration = (Vibration - Av)./sigmav;
Temp = (Temprature - At)./sigmat;

% calculate the group means
means = zeros(3,2);

for i = 1:2
    means(:,i) = [mean(Pressure(:, i));
                  mean(Vibration(:, i));
                  mean(Temp(:, i))];
end

%--------------------------------------------------------------------------
% Pressure vs Vibration
disp("---------------------------")
disp("For Pressure vs Vibration")
[eigenVectors, eigenValues] = LDA(Pressure, Vibration)

% get the index of the colum with the max variance
[~,vecIndex] = find(eigenValues == max(eigenValues,[],'all'));

n = 0;
figure(1);
    hold on 
for i = 1:2
    scatter(Pressure(:,i),Vibration(:,i),'filled', "DisplayName", names(i+1),"MarkerFaceColor", colours(i+1));
    n = n+10;
    plot(means(1,i),means(2,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(i+1), "DisplayName",names(i+1)+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')%, "FontSize",25)
axis([floor(min(Pressure, [], 'all')) ceil(max(Pressure, [], 'all')) floor(min(Vibration, [], 'all')) ceil(max(Vibration, [], 'all'))])
xlabel('Pressure')%, "FontSize",25)
ylabel('Vibreaiton')%, "FontSize")%,25)
title("Pressure Vs Vibration for the 'Black foam' and 'Car sponge' with the LDA function")
%--------------------------------------------------------------------------
% Pressure vs Temp
disp("---------------------------")
disp("For Pressure vs Temperature")
[eigenVectors, eigenValues] = LDA(Pressure, Temp)

% get the index of the colum with the max variance
[~,vecIndex] = find(eigenValues == max(eigenValues,[],'all'));

n = 0;
figure(2);
    hold on 
for i = 1:2
    scatter(Pressure(:,i),Temp(:,i),'filled', "DisplayName", names(i+1),"MarkerFaceColor", colours(i+1));
    n = n+10;
    plot(means(1,i),means(3,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(i+1), "DisplayName",names(i+1)+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')%, "FontSize",25)
axis([floor(min(Pressure, [], 'all')) ceil(max(Pressure, [], 'all')) ...
    floor(min(Temp, [], 'all')) ceil(max(Temp, [], 'all'))])
xlabel('Pressure')%, "FontSize",25)
ylabel('Temperature')%, "FontSize",25)
title("Pressure Vs Temperature for the 'Black foam' and 'Car sponge' with the LDA function")
%--------------------------------------------------------------------------
% Temp vs Vibration
disp("---------------------------")
disp("For Temperature vs Vibrations ")
[eigenVectors, eigenValues] = LDA(Temp, Vibration)

% get the index of the colum with the max variance
[~,vecIndex] = find(eigenValues == max(eigenValues,[],'all'));

n = 0;
figure(3);
    hold on 
for i = 1:2
    scatter(Temp(:,i),Vibration(:,i),'filled', "DisplayName", names(i+1),"MarkerFaceColor", colours(i+1));
    n = n+10;
    plot(means(3,i),means(2,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(i+1), "DisplayName",names(i+1)+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')%, "FontSize",25)
axis([floor(min(Temp, [], 'all')) ceil(max(Temp, [], 'all'))...
    floor(min(Vibration, [], 'all')) ceil(max(Vibration, [], 'all'))])
xlabel('Temperature')%, "FontSize",25)
ylabel('Vibrations')%, "FontSize",25)
title("Temprature Vs Vibrations for the 'Black foam' and 'Car sponge' with the LDA function")

%% Part B
disp("----------------------------")
disp("Part B")
disp("----------------------------")
% Apply LDA to the 3D Data
disp("---------------------------")
disp("For PVT Data")
[eigenVectors, eigenValues] = LDA3(Pressure, Vibration, Temp)

% get the index of the eigen vector with the largest Variance
[~,vec1Index] = find(eigenValues == max(eigenValues,[],'all'));
% get the index of the eigen vector with the second largest Variance

reduced = eigenValues;
reduced(:,vec1Index)= [];
[~, vec2Index] = find(eigenValues == max(reduced, [], 'all'));

% plot the data
figure(4);
hold on
for i = 1:2
    scatter3(Pressure(:, i),Vibration(:,i),Temp(:, i),'filled', "DisplayName", names(i+1),"MarkerFaceColor", colours(i+1));
    plot3(means(1,i),means(2,i),means(3,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(i+1), "DisplayName",names(i+1)+" mean");
end
plot3(0,0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
% plot the first LDA function
plot3([-eigenVectors(1,vec1Index),eigenVectors(1,vec1Index)].*10,...
      [-eigenVectors(2,vec1Index),eigenVectors(2,vec1Index)].*10,...
      [-eigenVectors(3,vec1Index),eigenVectors(3,vec1Index)].*10,...
      'k--', "LineWidth",1, 'DisplayName','First LDA function');
% plot the first LDA function
plot3([-eigenVectors(1,vec2Index),eigenVectors(1,vec2Index)].*10,...
      [-eigenVectors(2,vec2Index),eigenVectors(2,vec2Index)].*10,...
      [-eigenVectors(3,vec2Index),eigenVectors(3,vec2Index)].*10,..., 
      '--', "Color",grey, "LineWidth",1, 'DisplayName','Second LDA function');

patch([-eigenVectors(1,vec1Index),-eigenVectors(1,vec2Index),eigenVectors(1,vec1Index),eigenVectors(1,vec2Index)].*100,...
      [-eigenVectors(2,vec1Index),-eigenVectors(2,vec2Index),eigenVectors(2,vec1Index),eigenVectors(2,vec2Index)].*100,...
      [-eigenVectors(3,vec1Index),-eigenVectors(3,vec2Index),eigenVectors(3,vec1Index),eigenVectors(3,vec2Index)].*100,...
      'k', 'FaceColor','#0072BD','FaceAlpha', 0.5, "DisplayName", '', 'HandleVisibility', 'off');
hold off

legend('Location','bestoutside')%, "FontSize",15)
grid()
view(3)

axis([floor(min(Pressure, [], 'all')) ceil(max(Pressure, [], 'all'))...
      floor(min(Vibration, [], 'all')) ceil(max(Vibration, [], 'all'))...
      floor(min(Temp, [], 'all')) ceil(max(Temp, [], 'all'))])

xlabel("Pressure")%, "FontSize",15);
ylabel("Vibration")%, "FontSize",15);
zlabel("Temperature")%, "FontSize",15);

title("PVT data for the 'Black foam' and 'Car sponge' with the LDA function")

%% Part C
disp("----------------------------")
disp("Part C")
disp("----------------------------")

% pull out the data for the Acrylic and the Steel vase
Pressure = [pressure(:,1),pressure(:,6)];
Vibration = [vibration(:,1),vibration(:,6)];
Temprature = [temprature(:,1),temprature(:,6)];

% standadise the data
Ap = mean(Pressure, 'all');
Av = mean(Vibration, 'all');
At = mean(Temprature, 'all');

sigmap = std(Pressure, 0,'all');
sigmav = std(Vibration, 0,'all');
sigmat = std(Temprature, 0,'all');

Pressure = (Pressure - Ap)./sigmap;
Vibration = (Vibration - Av)./sigmav;
Temp = (Temprature - At)./sigmat;

% calculate the group means
means = zeros(3,2);
for i = 1:2
    means(:,i) = [mean(Pressure(:, i));
                  mean(Vibration(:, i));
                  mean(Temp(:, i))];
end
%--------------------------------------------------------------------------
% Pressure vs Vibration
disp("---------------------------")
disp("For Pressure vs Vibration")
[eigenVectors, eigenValues] = LDA(Pressure, Vibration)

% get the index of the colum with the max variance
[~,vecIndex] = find(eigenValues == max(eigenValues,[],'all'));
Indexes = [1, 6];
n = 0;
figure(5);
    hold on 
for i = 1:2
    scatter(Pressure(:,i),Vibration(:,i),'filled', "DisplayName", names(Indexes(i)),"MarkerFaceColor", colours(Indexes(i)));
    n = n+10;
    plot(means(1,i),means(2,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(Indexes(i)), "DisplayName",names(Indexes(i))+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')%, "FontSize",25)
axis([floor(min(Pressure, [], 'all')) ceil(max(Pressure, [], 'all')) floor(min(Vibration, [], 'all')) ceil(max(Vibration, [], 'all'))])
xlabel('Pressure')%, "FontSize",25)
ylabel('Vibreaiton')%, "FontSize",25)
title("Pressure Vs Vibration for the 'Acrylic' and 'Steel vase' with the LDA function")

%--------------------------------------------------------------------------
% Pressure vs Temp
disp("---------------------------")
disp("For Pressure vs Temperature")
[eigenVectors, eigenValues] = LDA(Pressure, Temp)

% get the index of the colum with the max variance
[~,vecIndex] = find(eigenValues == max(eigenValues,[],'all'));

n = 0;
figure(6);
    hold on 
for i = 1:2
    scatter(Pressure(:,i),Temp(:,i),'filled', "DisplayName", names(Indexes(i)),"MarkerFaceColor", colours(Indexes(i)));
    n = n+10;
    plot(means(1,i),means(3,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(Indexes(i)), "DisplayName",names(Indexes(i))+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')%, "FontSize",25)
axis([floor(min(Pressure, [], 'all')) ceil(max(Pressure, [], 'all')) ...
    floor(min(Temp, [], 'all')) ceil(max(Temp, [], 'all'))])
xlabel('Pressure')%, "FontSize",25)
ylabel('Temperature')%, "FontSize",25)
title("Pressure Vs Temprature for the 'Acrylic' and 'Steel vase' with the LDA function")
%--------------------------------------------------------------------------
% Temp vs Vibration
disp("---------------------------")
disp("For Temperature vs Vibrations ")
[eigenVectors, eigenValues] = LDA(Temp, Vibration)

% get the index of the colum with the max variance
[~,vecIndex] = find(eigenValues == max(eigenValues,[],'all'));

n = 0;
figure(7);
    hold on 
for i = 1:2
    scatter(Temp(:,i),Vibration(:,i),'filled', "DisplayName", names(Indexes(i)),"MarkerFaceColor", colours(Indexes(i)));
    n = n+10;
    plot(means(3,i),means(2,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(Indexes(i)), "DisplayName",names(Indexes(i))+" mean");
end
plot(0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
plot([-eigenVectors(1,vecIndex),eigenVectors(1,vecIndex)].*10,[-eigenVectors(2,vecIndex),eigenVectors(2,vecIndex)].*10, 'k', "LineWidth",1, 'DisplayName','LDA function')
hold off 
grid()
legend('Location','bestoutside')%, "FontSize",25)
axis([floor(min(Temp, [], 'all')) ceil(max(Temp, [], 'all'))...
    floor(min(Vibration, [], 'all')) ceil(max(Vibration, [], 'all'))])
xlabel('Temperature')%, "FontSize",25)
ylabel('Vibrations')%, "FontSize",25)
title("Temprature Vs Vibrations for the 'Acrylic' and 'Steel vase' with the LDA function")

%--------------------------------------------------------------------------
% PVT

% Apply LDA to the 3D Data
disp("---------------------------")
disp("For PVT Data")
[eigenVectors, eigenValues] = LDA3(Pressure, Vibration, Temp)

% get the index of the eigen vector with the largest Variance
[~,vec1Index] = find(eigenValues == max(eigenValues,[],'all'));
% get the index of the eigen vector with the second largest Variance

reduced = eigenValues;
reduced(:,vec1Index)= [];
[~, vec2Index] = find(eigenValues == max(reduced, [], 'all'));

% plot the data
figure(8);
hold on
for i = 1:2
    scatter3(Pressure(:, i),Vibration(:,i),Temp(:, i),'filled', "DisplayName", names(Indexes(i)),"MarkerFaceColor", colours(Indexes(i)));
    plot3(means(1,i),means(2,i),means(3,i),'o', 'Color','k', "MarkerSize",10, "MarkerFaceColor",colours(Indexes(i)), "DisplayName",names(Indexes(i))+" mean");
end
plot3(0,0,0,'ko',"MarkerSize",10, 'DisplayName','Group Mean');
% plot the first LDA function
plot3([-eigenVectors(1,vec1Index),eigenVectors(1,vec1Index)].*10,...
      [-eigenVectors(2,vec1Index),eigenVectors(2,vec1Index)].*10,...
      [-eigenVectors(3,vec1Index),eigenVectors(3,vec1Index)].*10,...
      'k--', "LineWidth",1, 'DisplayName','First LDA function');
% plot the first LDA function
plot3([-eigenVectors(1,vec2Index),eigenVectors(1,vec2Index)].*10,...
      [-eigenVectors(2,vec2Index),eigenVectors(2,vec2Index)].*10,...
      [-eigenVectors(3,vec2Index),eigenVectors(3,vec2Index)].*10,..., 
      '--', "Color",grey, "LineWidth",1, 'DisplayName','Second LDA function');

patch([-eigenVectors(1,vec1Index),-eigenVectors(1,vec2Index),eigenVectors(1,vec1Index),eigenVectors(1,vec2Index)].*100,...
      [-eigenVectors(2,vec1Index),-eigenVectors(2,vec2Index),eigenVectors(2,vec1Index),eigenVectors(2,vec2Index)].*100,...
      [-eigenVectors(3,vec1Index),-eigenVectors(3,vec2Index),eigenVectors(3,vec1Index),eigenVectors(3,vec2Index)].*100,...
      'k', 'FaceColor','#77AC30','FaceAlpha', 0.5, "DisplayName", '', 'HandleVisibility', 'off');
hold off

legend('Location','bestoutside')%, "FontSize",25)
grid()
view(3)

axis([floor(min(Pressure, [], 'all')) ceil(max(Pressure, [], 'all'))...
      floor(min(Vibration, [], 'all')) ceil(max(Vibration, [], 'all'))...
      floor(min(Temp, [], 'all')) ceil(max(Temp, [], 'all'))])

xlabel("Pressure")%, "FontSize",25);
ylabel("Vibration")%, "FontSize",25);
zlabel("Temperature")%, "FontSize",25);

title("PVT data for the 'Black foam' and 'Car sponge' with the LDA function")





