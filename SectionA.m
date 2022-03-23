%% load all the data 
slash = '/'; %<------ use this to change all "/" to "\" or visa versa
dataDir = 'PR_CW_DATA_2021';
Data = dir(fullfile(dataDir,'*.mat'));

fprintf("There are %d sets of data.\n", length(Data))
%% loop through and plot the data

% define the colours ... because I'm a bit extra
orange = '#ff7f0e';
grey = '#4f4f4f';

row = 1;
figure(1);
for i = 2:10:60
    name = [Data(i).folder, slash ,Data(i).name];
    load(name);

    % do some string formatting to add to the plots
    endName = strfind(Data(i).name, '_');
    endName = endName(1) - 1;
    objectName = Data(i).name(1:endName);

    % numper of points
    n = 1:length(F0tac);
    %The pressure data

    subplot(6,4,row);
    hold on;
    plot(n, F0pdc,'Color',grey);
    plot(n, F1pdc,'Color',orange);
    hold off;
    grid();
    ylabel(objectName)
    title("Pressure");

    legend("F0", "F1");
   
    %The vibration data
    figure(1);
    subplot(6,4,row+1)
    hold on;
    plot(n, F0pac(2,:),'Color',grey);
    plot(n, F1pac(2,:),'Color',orange);
    hold off;
    grid();
    title("Vibration");
    legend("F0", "F1");

    %The temprature data
    figure(1);
    subplot(6,4,row+2)
    hold on;
    plot(n, F0tac,'Color',grey);
    plot(n, F1tac,'Color',orange);
    hold off;
    grid();
    title("Temprature");
    legend("F0", "F1");

    %The Electrode data
    figure(1);
    subplot(6,4,row+3)
    hold on;
    plot(n, F0Electrodes,'Color',grey);
    plot(n, F1Electrodes,'Color',orange);
    hold off;
    grid();
    title("Electrodes");
    legend("F0", "F1");
    row = row + 4;
end
%%
% I think n = 400 will be a good time interval 
N = 700;

% extract all the data and then split in objects later
p = [];
v = [];
t = [];
e = [];

% extract the data
for i = 1:60
   name = [Data(i).folder, slash ,Data(i).name];
   load(name);
   p = [p; F0pdc(N)];
   v = [v; F0pac(2,N)];
   t = [t; F0tac(N)];
   e = [e, F0Electrodes(:,N)];
end

% split the data in to the different objects. 
% The data will be represented as matrcies were the colums relate to an object

pressure = [];
vibration = [];
temprature = [];

% the electrode data is already split
electordes = e;

for n = 0:10:50
    pressure = [pressure, p(n+1:n+10)];
    vibration = [vibration, v(n+1:n+10)];
    temprature = [temprature, t(n+1:n+10)];
end


% NOTE: for the electrode we know which object data is associated with 
% baced off their coloum i.e colums 1 to 10 are for one object and 21 to 30
% are for another 


PVT = [transpose(makeVector(Pressure));transpose(makeVector(Vibration));transpose(makeVector(Temp))];

%% Make the sater plot

% colours
colours = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#A2142F"];
names = ["Acrylic", "Black foam","Car sponge", "Flour sack", "Kitchen sponge","Steel vase"];
figure(2);
view(3)
hold on
for i = 1:6
    scatter3(pressure(i,:),vibration(i,:),temprature(i,:),'filled', "DisplayName", names(i),"MarkerFaceColor", colours(i));
end
hold off

xlabel("Pressure");
ylabel("Vibration");
zlabel("Temprature");
title("Scatter plot of the PVT data for the different objects.");

grid()
legend()
