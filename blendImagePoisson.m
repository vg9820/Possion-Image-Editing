function imret = blendImagePoisson(im1, im2, roi, targetPosition)

% input: im1 (background), im2 (foreground), roi (in im2), targetPosition (in im1)
[m,n,d] = size(im1);
% exchange the column of roi and target,to make sure roi and im1 have same
% cordinate system!!!
% as we use inpolygon,we have to make vertex clockwise or counterclockwise
roi(:,[1,2]) = roi(:,[2,1]);
roi = [roi;roi(1,:)];
targetPosition(:,[1,2]) = targetPosition(:,[2,1]);
targetPosition = [targetPosition;targetPosition(1,:)];
M = [0,0,m,m,0];
N = [0,n,n,0,0];
%dx dy denote the displacement from roi to target
dx = targetPosition(1,1)-roi(1,1);
dy = targetPosition(1,2)-roi(1,2);
%% Construct Sparse Matrix
%here im2 is 3 unit8,so we need 3 matrixs to hold values
[A1,A2,A3] = deal(sparse(m*n,m*n));
[b1,b2,b3] = deal(zeros(m*n,1));

for i=1:m
    for j=1:n
        % check if [i,j] is inside targetPosition
        in = inpolygon(i,j,targetPosition(:,1),targetPosition(:,2));
        if in==0
            % in=0 means outside the Polygon
            [A1((i-1)*n+j,(i-1)*n+j),A2((i-1)*n+j,(i-1)*n+j),A3((i-1)*n+j,(i-1)*n+j)]=deal(1,1,1);
            [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(im1(i,j,1),im1(i,j,2),im1(i,j,3));
            continue
        end
        % for [i,j] inside targetPosition
        % Np denotes the cardinate of Np
        Np = 0;
        %[i0,j0] indicates the translation map from im2 to im1
        [i0,j0]=deal(round(i-dx),round(j-dy));
        
        %[i-1,j]
        in1 = inpolygon(i-1,j,M,N);
        if in1 == 1
            % [i-1.j] is interior the im1
            Np = Np+1;
            in2 = inpolygon(i-1,j,targetPosition(:,1),targetPosition(:,2));
            if in2 == 1
                % [i-1,j] is interior the targetPosition
                [A1((i-1)*n+j,(i-2)*n+j),A2((i-1)*n+j,(i-2)*n+j),A3((i-1)*n+j,(i-2)*n+j)]=deal(-1,-1,-1);
            else
                % [i-1,j] is exterior the targetPosition
                % here because [i,j] is interior to omega, [i-1.j] must be
                % boundary point
                [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(b1((i-1)*n+j)+double(im1(i-1,j,1)),b2((i-1)*n+j)+double(im1(i-1,j,2)),b3((i-1)*n+j)+double(im1(i-1,j,3)));
            end
            %here use double to prevent unit8 type from overflow.
            [v1,v2,v3] = deal(double(im2(i0,j0,1))-double(im2(i0-1,j0,1)),double(im2(i0,j0,2))-double(im2(i0-1,j0,2)),double(im2(i0,j0,3))-double(im2(i0-1,j0,3)));
            [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(b1((i-1)*n+j)+v1,b2((i-1)*n+j)+v2,b3((i-1)*n+j)+v3);
        end
        
         %[i+1,j]
        in1 = inpolygon(i+1,j,M,N);
        if in1 == 1
            Np = Np+1;
            in2 = inpolygon(i+1,j,targetPosition(:,1),targetPosition(:,2));
            if in2 == 1
                [A1((i-1)*n+j,i*n+j),A2((i-1)*n+j,i*n+j),A3((i-1)*n+j,i*n+j)]=deal(-1,-1,-1);
            else
                [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(b1((i-1)*n+j)+im1(i+1,j,1),b2((i-1)*n+j)+im1(i+1,j,2),b3((i-1)*n+j)+im1(i+1,j,3));
            end
            [v1,v2,v3] = deal(double(im2(i0,j0,1))-double(im2(i0+1,j0,1)),double(im2(i0,j0,2))-double(im2(i0+1,j0,2)),double(im2(i0,j0,3))-double(im2(i0+1,j0,3)));
            [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(b1((i-1)*n+j)+v1,b2((i-1)*n+j)+v2,b3((i-1)*n+j)+v3);
        end
        
        %[i,j-1]
        in1 = inpolygon(i,j-1,M,N);
        if in1 == 1
            Np = Np+1;
            in2 = inpolygon(i,j-1,targetPosition(:,1),targetPosition(:,2));
            if in2 == 1
                [A1((i-1)*n+j,(i-1)*n+j-1),A2((i-1)*n+j,(i-1)*n+j-1),A3((i-1)*n+j,(i-1)*n+j-1)]=deal(-1,-1,-1);
            else
                [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(b1((i-1)*n+j)+double(im1(i,j-1,1)),b2((i-1)*n+j)+double(im1(i,j-1,2)),b3((i-1)*n+j)+double(im1(i,j-1,3)));
            end
            [v1,v2,v3] = deal(double(im2(i0,j0,1))-double(im2(i0,j0-1,1)),double(im2(i0,j0,2))-double(im2(i0,j0-1,2)),double(im2(i0,j0,3))-double(im2(i0,j0-1,3)));
            [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(b1((i-1)*n+j)+v1,b2((i-1)*n+j)+v2,b3((i-1)*n+j)+v3);
        end
        
        %[i,j+1]
        in1 = inpolygon(i,j+1,M,N);
        if in1 == 1
            Np = Np+1;
            in2 = inpolygon(i,j+1,targetPosition(:,1),targetPosition(:,2));
            if in2 == 1
                [A1((i-1)*n+j,(i-1)*n+j+1),A2((i-1)*n+j,(i-1)*n+j+1),A3((i-1)*n+j,(i-1)*n+j+1)]=deal(-1,-1,-1);
            else
                [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(b1((i-1)*n+j)+double(im1(i,j+1,1)),b2((i-1)*n+j)+double(im1(i,j+1,2)),b3((i-1)*n+j)+double(im1(i,j+1,3)));
            end
            [v1,v2,v3] = deal(double(im2(i0,j0,1))-double(im2(i0,j0+1,1)),double(im2(i0,j0,2))-double(im2(i0,j0+1,2)),double(im2(i0,j0,3))-double(im2(i0,j0+1,3)));
            [b1((i-1)*n+j),b2((i-1)*n+j),b3((i-1)*n+j)]=deal(b1((i-1)*n+j)+v1,b2((i-1)*n+j)+v2,b3((i-1)*n+j)+v3);
        end
        %we get Np
        [A1((i-1)*n+j,(i-1)*n+j),A2((i-1)*n+j,(i-1)*n+j),A3((i-1)*n+j,(i-1)*n+j)]=deal(Np,Np,Np);
    end
end
%% LU,QR decomposition using linsolve:don't use this because full(A) require too much space(156445x156445 (182.4GB)).
% %As A is not Symmetric positive definite matrix,we cannot use chol here.
% [L1,U1] = lu(A1);
% [Q2,R2] = qr(A2);
% [L3,U3] = lu(A3);
% % linsolve get X1
% tic
% opts.LT = true;
% Y1 = linsolve(full(L1),b1,opts);
% opts.UT = true;
% X1 = linsolve(full(U1),Y1,opts);
% toc
% % linsolve get X2
% tic
% %Q is not  opts in linsolve, so we use \.
% Y2 = Q2\b2;
% opts.UT = true;
% X2 = linsolve(full(R2),Y2,opts);
% toc
% %X3 = Q3\(R3\b3); 
% tic
% opts.LT = true;
% Y3 = linsolve(full(L3),b3,opts);
% opts.UT = true;
% X3 = linsolve(full(U3),Y3,opts);
% toc
% % Compared to LU,computing X3 using QR really takes a long time here  
%% LU QR decompostion using '\'
% To be honest, 'A\b' is best method for most cases,it takes some steps to
% verify matrix A (require some time) and choose best method to compute(such as QR,LU,Chol)
% So here I do LU,QR decompostion on A is meaningless since  'A\b'will
% still judge A,even A is LT/UT. Do pre-decompostion is meaningful when we
% use linsolve and choose opts so that we can save the time for judge A.
% But at least.'A\b' accept sparse matrix,which is pleasing.
 [L1,U1] = lu(A1);
 [Q2,R2] = qr(A2);
 [L3,U3] = lu(A3);
 tic
 X1 = U1\(L1\b1);
 toc
 tic
 X2 = R2\(Q2\b2);
 toc
 tic
 X3 = U3\(L3\b3);
 toc
 %Comparing time comsuming, LU perform better than QR(4s to 23s)
%% TODO: compute blended image
imret = im1*0;
for i=1:m
    for j=1:n
        [imret(i,j,1),imret(i,j,2),imret(i,j,3)]=deal(X1((i-1)*n+j),X2((i-1)*n+j),X3((i-1)*n+j));
    end
end
