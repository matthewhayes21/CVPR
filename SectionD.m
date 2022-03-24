load F0_PVT.mat
PVT = [makeVector(pressure), makeVector(vibration), makeVector(temprature)];
[nrows ncol] = size(PVT);
objects = [];
rep_no = nrows/length(names);
for i = 1:length(names)
    objects = [objects  repelem(names(i),rep_no)];
end

objects = reshape(objects, [length(objects) 1]); % convert to column vector

% define "superset" of objects i.e. the materials
materials = objects;
materials(materials == "Black foam" | materials == "Car sponge" | materials == "Kitchen sponge" ) = "Foam";
materials(materials == "Flour sack") = "Fabric";
materials(materials == "Acrylic") = "Plastic";
materials(materials == "Steel vase") = "Metal";

material_names = ["Foam", "Fabric", "Plastic", "Metal"];

%% L2 Norm
dist_metric = 'euclidean';
X = PVT;
Y = pdist(X,dist_metric);
squareform(Y);
Z = linkage(Y);
dendrogram(Z,0)
figure(1);
H = dendrogram(Z,0,'ColorThreshold',29,'Labels',objects,"Orientation","right");
set(gca,'FontSize',6.5)
% notes: for L2 norm cutoff of 29 will produce 6 object classes
% and 36 will produce 4 object classes
figure(2);
H = dendrogram(Z,0,'ColorThreshold',36,'Labels',materials,"Orientation","right");
set(gca,'FontSize',6.5)

L2_6clusters = cluster(Z,'maxclust',6);
L2_4cluster = cluster(Z,'maxclust',4);
% Clusters for L2 (objects)
cluster_labels = [L2_6clusters objects];
[~,idx] = sort(cluster_labels(:,1)); % sort just the first column
cluster_labels_sorted = cluster_labels(idx,:);   % sort the whole matrix using the sort indice

clusters = unique(cluster_labels(:,1));
y = [];

% get counts of the objects for each cluster
for i = 1:length(clusters)
    clust = clusters(i);
    % get sub matrix corresponding to the matrix class
    sub_cluster_indx = find(cluster_labels(:,1) == clust); 
    sub_cluster_labels = cluster_labels(sub_cluster_indx,:);
    histo_objects = [];
    for j = 1:length(names)
        object = names(j);
        object_count = sum(sub_cluster_labels(:,2) == object);
        histo_objects = [histo_objects object_count];
    end
    
    y = [y; histo_objects];
end

figure(3);
colours = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#A2142F"];
bh = bar(y,'stacked');
set(bh, 'FaceColor', 'Flat')
set(gca,'FontSize',20)
for i = 1:length(colours)
    bh(i).FaceColor = colours(i);
end
% bh(1).FaceColor = colours(i);
% bh(2).CData = colours;
% bh(3).CData = colours;
% bh(4).CData = colours;
% bh(5).CData = colours;
% bh(6).CData = colours;
% set(bh, {'CData'}, colours)
xlabel('Class');
ylabel('Object count')
legend(names)

% Clusters for L2 (materials)
cluster_labels = [L2_4cluster materials];
[~,idx] = sort(cluster_labels(:,1)); % sort just the first column
cluster_labels_sorted = cluster_labels(idx,:);   % sort the whole matrix using the sort indice

clusters = unique(cluster_labels(:,1));
y = [];

% get counts of the objects for each cluster
for i = 1:length(clusters)
    clust = clusters(i);
    % get sub matrix corresponding to the matrix class
    sub_cluster_indx = find(cluster_labels(:,1) == clust); 
    sub_cluster_labels = cluster_labels(sub_cluster_indx,:);
    histo_objects = [];
    for j = 1:length(material_names)
        object = material_names(j);
        object_count = sum(sub_cluster_labels(:,2) == object);
        histo_objects = [histo_objects object_count];
    end
    
    y = [y; histo_objects];
end

figure(4);
colours = ["#000000","#FFFF00","#FF00FF","#00FF00"];
bh = bar(y,'stacked');
set(bh, 'FaceColor', 'Flat')
set(gca,'FontSize',20)
for i = 1:length(colours)
    bh(i).FaceColor = colours(i);
end
xlabel('Class');
ylabel('Material count')
legend(material_names)
%% L1 Norm
dist_metric = 'cityblock';
X = PVT;
Y = pdist(X,dist_metric);
squareform(Y);
Z = linkage(Y);
dendrogram(Z,0)
figure(5);
H = dendrogram(Z,0,'ColorThreshold',40,'Labels',objects,"Orientation","right");
set(gca,'FontSize',6.5)
% notes: for L2 norm cutoff of 29 will produce 6 object classes
% and 50 will produce 4 object classes
figure(6);
H = dendrogram(Z,0,'ColorThreshold',50,'Labels',materials,"Orientation","right");
set(gca,'FontSize',6.5)

L2_6clusters = cluster(Z,'maxclust',6);
L2_4cluster = cluster(Z,'maxclust',4);
% Clusters for L2 (objects)
cluster_labels = [L2_6clusters objects];
[~,idx] = sort(cluster_labels(:,1)); % sort just the first column
cluster_labels_sorted = cluster_labels(idx,:);   % sort the whole matrix using the sort indice

clusters = unique(cluster_labels(:,1));
y = [];

% get counts of the objects for each cluster
for i = 1:length(clusters)
    clust = clusters(i);
    % get sub matrix corresponding to the matrix class
    sub_cluster_indx = find(cluster_labels(:,1) == clust); 
    sub_cluster_labels = cluster_labels(sub_cluster_indx,:);
    histo_objects = [];
    for j = 1:length(names)
        object = names(j);
        object_count = sum(sub_cluster_labels(:,2) == object);
        histo_objects = [histo_objects object_count];
    end
    
    y = [y; histo_objects];
end

figure(7);
colours = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#A2142F"];
bh = bar(y,'stacked');
set(bh, 'FaceColor', 'Flat')
set(gca,'FontSize',20)
for i = 1:length(colours)
    bh(i).FaceColor = colours(i);
end
% bh(1).FaceColor = colours(i);
% bh(2).CData = colours;
% bh(3).CData = colours;
% bh(4).CData = colours;
% bh(5).CData = colours;
% bh(6).CData = colours;
% set(bh, {'CData'}, colours)
xlabel('Class');
ylabel('Object count')
legend(names)

% Clusters for L2 (materials)
cluster_labels = [L2_4cluster materials];
[~,idx] = sort(cluster_labels(:,1)); % sort just the first column
cluster_labels_sorted = cluster_labels(idx,:);   % sort the whole matrix using the sort indice

clusters = unique(cluster_labels(:,1));
y = [];

% get counts of the objects for each cluster
for i = 1:length(clusters)
    clust = clusters(i);
    % get sub matrix corresponding to the matrix class
    sub_cluster_indx = find(cluster_labels(:,1) == clust); 
    sub_cluster_labels = cluster_labels(sub_cluster_indx,:);
    histo_objects = [];
    for j = 1:length(material_names)
        object = material_names(j);
        object_count = sum(sub_cluster_labels(:,2) == object);
        histo_objects = [histo_objects object_count];
    end
    
    y = [y; histo_objects];
end

figure(8);
colours = ["#000000","#FFFF00","#FF00FF","#00FF00"];
bh = bar(y,'stacked');
set(bh, 'FaceColor', 'Flat')
set(gca,'FontSize',20)
for i = 1:length(colours)
    bh(i).FaceColor = colours(i);
end
xlabel('Class');
ylabel('Material count')
legend(material_names)

%% Part 2 Bagging
% get the PCA components
pca3 = eProjected(:,1:3);
pca3_data = [pca3 objects];

% split data into training and testing
rng(7); % original number is 7
[nrows ncol] = size(pca3);
train_indx = randsample(nrows,0.6*nrows);
test_indx = setdiff(1:nrows, train_indx);
pca3_train = pca3(train_indx,:);
objects_train = objects(train_indx);
pca3_test = pca3(test_indx,:);
objects_test = objects(test_indx);

%%

% electrodes_data = [electrodes_data objects];

X = pca3_train;
Y = categorical(objects_train);

rng(11); % For reproducibility 
Mdl = TreeBagger(50,X,Y,'OOBPrediction','On',...
    'Method','classification');

% view(Mdl.Trees{1},'Mode','graph');

figure(9);
oobErrorBaggedEnsemble = oobError(Mdl);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';

%% Retrain with "optimal" number of trees and test
no_trees = 13;
Mdl = TreeBagger(no_trees,X,Y,'OOBPrediction','On',...
    'Method','classification');
Yfit = predict(Mdl,pca3_test);
figure(10);
cm = confusionchart(categorical(objects_test),categorical(Yfit))
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
set(gca,'FontSize',15)
%%
%  get the average accuracy
mean(categorical(objects_test) == categorical(Yfit))

%% Visualize 2 of the trees
figure(11);
view(Mdl.Trees{3},'Mode','graph');
figure(12);
view(Mdl.Trees{7},'Mode','graph');



%% Repeat the process but with all the electrode data

% get electrodes data to pca performance
electrodes_data = electordes.'; % transpose to get a vector of 60x19

X = electrodes_data(train_indx,:);
Y = categorical(objects_train);

rng(11); % For reproducibility
Mdl = TreeBagger(50,X,Y,'OOBPrediction','On',...
    'Method','classification');

% view(Mdl.Trees{1},'Mode','graph');

figure(13);
oobErrorBaggedEnsemble = oobError(Mdl);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';
%% Retrain with "optimal" number of trees and test
no_trees = 13;
Mdl = TreeBagger(no_trees,X,Y,'OOBPrediction','On',...
    'Method','classification');
Yfit = predict(Mdl,electrodes_data(test_indx,:));
figure(14);
cm = confusionchart(categorical(objects_test),categorical(Yfit))
cm.RowSummary = 'row-normalized';
cm.ColumnSummary = 'column-normalized';
set(gca,'FontSize',15)

%%
%  get the average accuracy
mean(categorical(objects_test) == categorical(Yfit))

%% Repeated tests to determine average performance using PCA components and All components
N = 100;
no_trees = 13;
average_accuracy_pca = zeros(N,1);
average_accuracy_all = zeros(N,1);

for i = 1:N
    % sample a different set for test and train
    [nrows ncol] = size(pca3);
    train_indx = randsample(nrows,0.6*nrows);
    test_indx = setdiff(1:nrows, train_indx);
    pca3_train = pca3(train_indx,:);
    objects_train = objects(train_indx);
    pca3_test = pca3(test_indx,:);
    objects_test = objects(test_indx);
    
    % run on pca3 data
    X_pca = pca3_train;
    Y = categorical(objects_train);
    Mdl_pca = TreeBagger(no_trees,X_pca,Y,'OOBPrediction','On',...
    'Method','classification');
    Yfit_pca = predict(Mdl_pca,pca3_test);

    average_accuracy_pca(i) = mean(categorical(objects_test) == categorical(Yfit_pca));

    % run on all the data
    X_all = electrodes_data(train_indx,:);
%     Y = categorical(objects_train);
    Mdl_all = TreeBagger(no_trees,X_all,Y,'OOBPrediction','On',...
    'Method','classification');
    Yfit_all = predict(Mdl_all,electrodes_data(test_indx,:));
    average_accuracy_all(i) = mean(categorical(objects_test) == categorical(Yfit_all));
end

fprintf("The average accuracy for the principal component trained bagged model is %f \n",mean(average_accuracy_pca))
fprintf("The average accuracy for the complete trained bagged model is %f \n",mean(average_accuracy_all))
%%



