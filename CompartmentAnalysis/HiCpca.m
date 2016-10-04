function pcatop5 = HiCpca(correlationmat)

pcaresults = pca(correlationmat);

pcatop5 = pcaresults(:,1:5);