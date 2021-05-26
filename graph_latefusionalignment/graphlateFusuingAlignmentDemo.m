clear
clc
warning off;

DataName{1} = 'caltech101_nTrain5_48';
DataName{2} = 'caltech101_nTrain10_48';
DataName{3} = 'caltech101_nTrain15_48';
DataName{4} = 'caltech101_nTrain20_48';
DataName{5} = 'caltech101_nTrain25_48';
DataName{6} = 'caltech101_nTrain30_48';
DataName{7} = 'CCV';
DataName{8} = 'ORL';
DataName{9} = 'flower17';
DataName{10} = 'YALE';
DataName{11} = 'plant';
DataName{12} = 'mfeat';
DataName{13} = 'AR10P';
DataName{14} = 'bbcsport2view';
DataName{15} = 'flower102';
DataName{16} = 'proteinfold';
DataName{17} = 'UCI_DIGIT';
DataName{18} = 'nonpl';

pathdata = 'E:\work2015';
datapath = '/home/zoey/#abababa/';
path = './';
addpath(genpath(path));

%%addpath('C:\Program Files\Mosek\8\toolbox\r2014a');

for i=11
    dataName = DataName{i}; %%% flower17; flower102; proteinFold,caltech101_mit,UCI_DIGIT,ccv
    %% %% washington; wisconsin; texas; cornell
    %% caltech101_nTrain5_48
    %% proteinFold
    load([datapath,'datasets/',dataName,'_Kmatrix'],'KH','Y');
    % load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numclass = length(unique(Y));
    numker = size(KH,3);
    num = size(KH,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    KH = kcenter(KH);
    KH = knorm(KH);
    K0 = zeros(num,num);

    qnorm = 2;

    opt.disp = 0;


    %%---The Proposed1---%%
    tic
    lambdaset = 2.^[-5:1:5];
    betaset = 10.^[-5:2:5];

    for it =1:length(lambdaset)     
        for ij = 1:length(betaset)
            [H_normalized9,iter] = graphlatefusionalignmentclustering(KH,numclass,lambdaset(it),betaset(ij),Y);
            res9 = myNMIACCwithmean(H_normalized9,Y,numclass);
            accval9(it,ij) = res9(1,1);
            accstd9(it,ij) = res9(1,2);
            nmival9(it,ij)= res9(2,1);
            nmistd9(it,ij) = res9(2,2);
            purval9(it,ij) = res9(3,1);
            purstd9(it,ij) = res9(3,2);
            iter9(it,ij) = iter;
        end
    end
    res(:,9) = [max(max(accval9)); max(max(nmival9));max(max(purval9))];

    save(['./results/',dataName,'.mat'],'accval9','nmival9','purval9','iter9');
end


