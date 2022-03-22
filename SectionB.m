slash = '/'; %<------ use this to change all "/" to "\" or visa versa

colours = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#A2142F"];
names = ["Acrylic", "Black foam","Car sponge", "Flour sack", "Kitchen sponge","Steel vase"];

% define the colours ... because I'm a bit extra
orange = "#ff7f0e";
grey = "#4f4f4f";
black = "#000000";

%% part 1

%  standardise the data

Ap = mean(pressure);
At = mean(temprature2);
Av = mean(vibration);

Pressure = pressure - Ap;
Vibration = vibration - Av;
Temp = temprature2 - At;

% check the new means
fprintf("The new mean is (%1.2f,%1.2f,%1.2f)\n", mean(Pressure), mean(Temp), mean(Vibration));

%%
% compute the covariance matrcies
% S = cov([Pressure, Vibration, Temp])
% % calculate the eigne vectors of S
% [V, D] = eig(S)

[Vec,Projected,Val] = pca([Pressure, Vibration, Temp]);

disp("The eigen vecotrs");
Vec 
disp("The eigen values");
Val
%% visualise the data
figure(1);
% plot of the origional data with the eigen vecotrs
view(3)
hold on
n = 0;
for i = 1:6
    scatter3(Pressure(n+1:n+10),Vibration(n+1:n+10),Temp(n+1:n+10),'filled', "DisplayName", names(i),"MarkerFaceColor", colours(i));
    n = n+10;
end

% plot the axies
% plot3([0,0],[floor(min(Vibration)-1) ceil(max(Vibration)+1) ],[0,0],'k--');
% plot3([floor(min(Pressure)-1) ceil(max(Pressure)+1) ],[0,0],[0,0],'k--');
% plot3([0,0],[0,0],[floor(min(Temp)-1),ceil(max(Temp)+1)],'k--');
axis([floor(min(Pressure)-1) ceil(max(Pressure)+1) floor(min(Vibration)-1) ceil(max(Vibration)+1) floor(min(Temp)-1) ceil(max(Temp)+1)])
axis('equal')

% plot the eigen vectors
c = [orange; grey; black];
for i = 1:3
    name = "Principal Componant " +  string(i);
    plot3([0, Vec(1,i)*20],[0, Vec(2,i)*20],[0, Vec(3,i)*20], "Color",c(i), 'LineWidth', 2, 'DisplayName', name);
end
hold off


% other  stuff
xlabel("Pressure");
ylabel("Vibration");
zlabel("Temprature");
title("Scatter plot of the PVT data for the different objects.");

grid()
legend("Location","eastoutside")
%%
figure(2);
% plot of data reduced to two dimentions
hold on
n = 0;
for i = 1:6
    scatter(Projected(n+1:n+10, 1),Projected(n+1:n+10, 2) ,'filled', "DisplayName", names(i),"MarkerFaceColor", colours(i));
    n = n+10;
end
hold off
grid()
legend()

xlabel("First Principal Component");
ylabel("Second Principal Component");
title("Scatter plot of the PVT data projected on to the first two principal components");
%%
figure(3);
% plot the data on a single axis
n = 0;
for i = 1:6
    for j = 1:3
        subplot(3,1,j)
        hold on
        scatter(Projected(n+1:n+10, j),zeros(length(Projected(n+1:n+10, j))),'filled', "DisplayName", names(i),"MarkerFaceColor", colours(i));
        hold off
    end
    n = n+10;
end

for j = 1:3
    subplot(3,1,j)
    grid();
    title("Data projected on to proncipal companent "+string(j))
end

%% Part 2

% standadise the data
Ae = mean(electordes, 2);
Electordes = zeros(size(electordes));

for i = 1:width(Electordes)
    Electordes(:,i) = electordes(:,i)-Ae;
end

[eVec,eProjected, eVal] = pca(transpose(Electordes));

figure(4);
plot(1:1:length(eVal), eVal, 'o-')
axis([1  18 0 max(eVal)+10])
grid()
title("Scree Plot")
ylabel("Eignnvalue")
xlabel("Principal Component")

%%
figure(5);
% scatter polt of the electrod data projected on to the first three
% pricipal components
view(3);
grid()
n = 0;
hold on
for i = 1:6
    scatter3(eProjected(n+1:n+10, 1),eProjected(n+1:n+10, 2),eProjected(n+1:n+10, 3),'filled', "DisplayName", names(i),"MarkerFaceColor", colours(i));
    n = n+10;
end
hold off
legend()



