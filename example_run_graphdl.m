%http://www.cad.zju.edu.cn/home/dengcai/Data/ReproduceExp.html#GraphSC
%Head to head comparison of GraphSC and GraphDL.

clear;	

%Different faces under lighting conditions. Viz w/o normalization.

load('Yale_32x32.mat');
nClass = length(unique(gnd));

%Normalize each data vector to have L2-norm equal to 1  
fea = NormalizeFea(fea);


%Visualize a feature


%Clustering in the original space
rng('default');
label = litekmeans(fea,nClass,'Replicates',10);
MIhat = MutualInfo(gnd,label);
disp(['kmeans use all the features. MIhat: ',num2str(MIhat)]);
%kmeans in the original space. MIhat: 0.66064

%PCA reduction
eigvector = pca(fea, 'NumComponents', 64);
newfea = fea*eigvector;
newfea = NormalizeFea(newfea);
%Clustering in 64-dim PCA subspace
rng('default');
label = litekmeans(newfea,nClass,'Replicates',10);
MIhat = MutualInfo(gnd,label);
disp(['kmeans in the 64-dim PCA subspace. MIhat: ',num2str(MIhat)]);
%kmeans in the 64-dim PCA subspace. MIhat: 0.65946


%Construct graphs, and graph laplacians
rng('default');
options = [];
options.WeightMode = 'Cosine';

W = constructW(fea, options);
W_2 = constructW(newfea', options);
W_1 = constructW(fea', options);


L_manifold = graph_laplacian(W);
L_full_data = graph_laplacian(W_1);
L_data = graph_laplacian(W_2);

%Run GraphSC
nBasis = 30;
alpha = 1;
beta = 0.1;
nIters = 100;
rng('default');
warning('off', 'all');
%[B, S, stat] = GraphSC(fea', W, nBasis, alpha, beta, nIters); %'


%Set up GraphDL
%Set up dictionary learning params
params = struct();
params.Y = fea';
params.T = 3;
params.K = nBasis;
params.alpha = 0.0;
params.beta = 0.6;
params.mu = 0.08;
params.iternum = 20;


params.L = L_full_data;
params.Lc = L_manifold;
params.sdebug = 1;

[D,X] = DGRDL(params);


%Clustering 
rng('default');
label = litekmeans(S',nClass,'Replicates',10); %'
MIhat = MutualInfo(gnd,label);
disp(['Clustering in the ',num2str(nBasis),'-dim GraphSC subspace with ',num2str(nIters),' iterations. MIhat: ',num2str(MIhat)]);
%Example results:
%Clustering in the 30-dim GraphSC subspace with 15 iterations. MIhat: 0.76741
%Clustering in the 30-dim GraphSC subspace with 100 iterations. MIhat: 0.83143


%Clustering 
rng('default');
label = litekmeans(X',nClass,'Replicates',10); %'
MIhat = MutualInfo(gnd,label);
disp(['Clustering in the ',num2str(nBasis),'-dim GraphDL subspace with ',num2str(nIters),' iterations. MIhat: ',num2str(MIhat)]);
%Example results:
%Clustering in the 30-dim GraphDL subspace with 15 iterations. MIhat: 0.95706


