%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            misc functions            %
%from ksvdbox13, ksvd.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function err2 = reperror2(X,D,Gamma)

% compute in blocks to conserve memory
err2 = zeros(1,size(X,2));
blocksize = 2000;
for i = 1:blocksize:size(X,2)
  blockids = i : min(i+blocksize-1,size(X,2));
  err2(blockids) = sum((X(:,blockids) - D*Gamma(:,blockids)).^2);
end

end
