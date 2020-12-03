function [Hstar,iter] = graphlatefusionalignmentclustering(KH,k,lambda,beta,Y)

num = size(KH, 1); %the number of samples
numker = size(KH, 3); %m represents the number of kernels
maxIter = 10; %the number of iterations
H = zeros(num,k,numker);
G = eye(num,k,numker);
Z= eye(num);
opt.disp = 0;

flag = 1;
iter = 0;
while flag
    iter = iter +1;
    %the first step-- optimize H_i

   
    
    for p=1:numker
        A = KH(:,:,p)+lambda*(Z+Z'-eye(num)-Z*Z');
        [AP,~] = eigs(A, k, 'la', opt);
        H(:,:,p) = AP;
    end
    
    %%the second step-- optimize S
    %the second step-- optimize S_i
    K = zeros(num,num);
    for i=1:numker
        K = K + H(:,:,i)*H(:,:,i)';
    end
    
    %tmp = K/(K+(beta/lambda)*eye(num));
    tmp = (beta*Z + lambda*K)/(lamda*K + beta*eye(num));
    for ii=1:num
        idx = 1:num;
        idx(ii) = [];
        G(ii,idx) = EProjSimplex_new(tmp(ii,idx));
    end
    %%Z = UpdateS(K,0,beta/lambda,0);
    
    %the third step-- optimize S
    K = zeros(num,num);
    for i=1:numker
        K = K + H(:,:,i)*H(:,:,i)';
    end
    
    tmp = (beta\(beta+gamma))*G;
     for ii=1:num
        idx = 1:num;
        idx(ii) = [];
        Z(ii,idx) = EProjSimplex_new(tmp(ii,idx));
    end
    
      
    term1 =0;
    term2 =0;
    term3 =0;
    for j =1:numker
        term1 = term1+ trace(KH(:,:,j)-KH(:,:,j)*H(:,:,j)*H(:,:,j)'); 
        term2 = term2+ lambda*norm((H(:,:,j)-G*H(:,:,j)),'fro')^2;
        term3 = term3+ beta*norm(G-Z,'for')^2;
    end

   % term3 = beta*norm(Z,'fro')^2;
    term4 = gamma*norm(Z,'fro')^2;
   % terms(iter,:)=[term1,term2,term3, term1+term2,term2+term3];
   terms(iter,:)=[term1,term2,term3, term4,term1+term2,term2+term3,term3+term4];
    obj(iter) = term1+term2+term3+term4;
    
    if (iter>2) && (abs((obj(iter-1)-obj(iter))/(obj(iter-1)))<1e-6 || iter>maxIter)
%     if iter==maxIter
       Z= (Z+Z')/2;
       Z=Z-diag(Z);
       [Hstar,~] = eigs(Z, k, 'la', opt);
       flag =0;
    end
end
