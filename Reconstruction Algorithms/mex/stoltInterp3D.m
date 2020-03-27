function SkYkXkZ = stoltInterp3D(SkYkXk,k,KU)
SkYkXkZ = complex(zeros(size(KU)));
sizeKU2 = size(KU,2);
parfor ii = 1:size(KU,1)
    for jj = 1:sizeKU2
        tempS = squeeze(SkYkXk(ii,jj,:));
        tempKU = squeeze(KU(ii,jj,:));
        SkYkXkZ(ii,jj,:) = interp1(k(:),tempS(:),tempKU(:),'v5cubic');
    end
end