function [ind]=maxvol2(a)
    [n,r]=size(a);

    for i=1:n
        p(i)=i;
    end

    %Initialize
    %[l,u,p]=lu(a,'vector');
    ind = p(1:r);
    b = a(p,:);
    sbm = a(ind,:);
    z = b(r+1:n,:)*inv(sbm);
    %Start iterations
    niters = 20;
    gm = 2;
    eps = 1e-06; %Stopping accuracy
    iter = 0;

    while (iter <= niters);
        [mx,pv] = max((abs(z))');
        [mx0,pv0] = max(mx);
        if ( mx0 <= 1 + eps )
            break;
        end
        %Swap i0-th row with j0-th row
        i0 = pv0;
        j0 = pv(i0);
        i = p(i0+r);
        p(i0+r) = p(j0);
        p(j0) = i;
        bb = z(:,j0);
        bb(i0) = bb(i0)+1;
        cc = z(i0,:);
        cc(j0) = cc(j0)-1;
        z = z-bb*cc./z(i0,j0);
        iter = iter+1;
        ind = p(1:r);
        %fprintf('%e, %e \n',abs(det(b(ind,:))),abs(det(b(ind,:))*mx0))
        %fprintf('log(det): %f\n',log(abs(det(b(ind,:)))));
        %please check consistency
    end
    return
end

