function imret = blendImagePoisson1(im1, im2, roi, targetPosition)
% input: im1 (background), im2 (foreground), roi (in im2), targetPosition (in im1)
[m,n,d] = size(im1);
clc
% exchange the column of roi and target,to make sure roi and im1 have same cordinate system!!!
% as we use inpolygon,we have to make vertex clockwise or counterclockwise
roi(:,[1,2]) = roi(:,[2,1]);
roi = [roi;roi(1,:)];

targetPosition(:,[1,2]) = targetPosition(:,[2,1]);
targetPosition = [targetPosition;targetPosition(1,:)];

%dx dy denote the displacement from roi to target
dx = targetPosition(1,1)-roi(1,1);
dy = targetPosition(1,2)-roi(1,2);

%% Construct Sparse Matrix
%Find the point in targetPosition
A = repmat([1:m]',[1,n]);
B = repmat([1:n],[m,1]);
in = inpolygon(A,B,targetPosition(:,1),targetPosition(:,2));
%get these point position
[x,y] = find(in == 1);
% construct matrix
[A1,A2,A3] = deal(sparse(length(x),length(x)));
[b1,b2,b3] = deal(zeros(length(x),1));
%judge matrix: for (x(i),y(i)), whether its neighbor is inside the targetPosition
judge = ismember([x-1,y;x,y-1;x+1,y;x,y+1],[x,y],'row');
judge = reshape(judge,[length(x),4]);

for i=1:length(x)
    [i0,j0]=deal(round(x(i)-dx),round(y(i)-dy));
    [A1(i,i),A2(i,i),A3(i,i)] = deal(4);
    s = sum(judge(i,:));
    % for inner point
    if s == 4
        [t1,t2,t3,t4] = deal(find(x==x(i)-1 & y==y(i)),find(x==x(i) & y==y(i)-1),find(x==x(i)+1 & y==y(i)),find(x==x(i) & y==y(i)+1));
        [A1(i,[t1,t2,t3,t4]),A2(i,[t1,t2,t3,t4]),A3(i,[t1,t2,t3,t4])] = deal([-1,-1,-1,-1]);
        b1(i) = 4*double(im2(i0,j0,1))-double(im2(i0-1,j0,1))-double(im2(i0+1,j0,1))-double(im2(i0,j0-1,1))-double(im2(i0,j0+1,1));
        b2(i) = 4*double(im2(i0,j0,2))-double(im2(i0-1,j0,2))-double(im2(i0+1,j0,2))-double(im2(i0,j0-1,2))-double(im2(i0,j0+1,2));
        b3(i) = 4*double(im2(i0,j0,3))-double(im2(i0-1,j0,3))-double(im2(i0+1,j0,3))-double(im2(i0,j0-1,3))-double(im2(i0,j0+1,3));
        continue;
    end
    % for boundary point
    % (x-1,y) out
    if(judge(i,1)==0)
        b1(i)=b1(i)+double(im1(x(i)-1,y(i),1));
        b2(i)=b2(i)+double(im1(x(i)-1,y(i),2));
        b3(i)=b3(i)+double(im1(x(i)-1,y(i),3));
    else
        t1 = find(x==x(i)-1 & y==y(i));
        [A1(i,t1),A2(i,t1),A3(i,t1)] = deal(-1);
    end
    
    % (x,y-1) out
    if(judge(i,2)==0)
        b1(i)=b1(i)+double(im1(x(i),y(i)-1,1));
        b2(i)=b2(i)+double(im1(x(i),y(i)-1,2));
        b3(i)=b3(i)+double(im1(x(i),y(i)-1,3));
    else
        t2 = find(x==x(i) & y==y(i)-1);
        [A1(i,t2),A2(i,t2),A3(i,t2)] = deal(-1);
    end 
    
    % (x+1,y) out
    if(judge(i,3)==0)
        b1(i)=b1(i)+double(im1(x(i)+1,y(i),1));
        b2(i)=b2(i)+double(im1(x(i)+1,y(i),2));
        b3(i)=b3(i)+double(im1(x(i)+1,y(i),3));
    else
        t3 = find(x==x(i)+1 & y==y(i));
        [A1(i,t3),A2(i,t3),A3(i,t3)] = deal(-1);
    end
    
    % (x,y+1) out
    if(judge(i,4)==0)
        b1(i)=b1(i)+double(im1(x(i),y(i)+1,1));
        b2(i)=b2(i)+double(im1(x(i),y(i)+1,2));
        b3(i)=b3(i)+double(im1(x(i),y(i)+1,3));
    else
        t4 = find(x==x(i) & y==y(i)+1);
        [A1(i,t4),A2(i,t4),A3(i,t4)] = deal(-1);
    end
    %vpq
    b1(i) = b1(i)+4*double(im2(i0,j0,1))-double(im2(i0-1,j0,1))-double(im2(i0+1,j0,1))-double(im2(i0,j0-1,1))-double(im2(i0,j0+1,1));
    b2(i) = b2(i)+4*double(im2(i0,j0,2))-double(im2(i0-1,j0,2))-double(im2(i0+1,j0,2))-double(im2(i0,j0-1,2))-double(im2(i0,j0+1,2));
    b3(i) = b3(i)+4*double(im2(i0,j0,3))-double(im2(i0-1,j0,3))-double(im2(i0+1,j0,3))-double(im2(i0,j0-1,3))-double(im2(i0,j0+1,3));
end
        
%% LU QR decompostion using '\'
tic
[L1,U1] = lu(A1);
 X1 = U1\(L1\b1);
toc

[Q2,R2] = qr(A2);
X2 = R2\(Q2\b2);


X3 = A3\b3;


 %Comparing time comsuming, LU perform better than QR(4s to 23s)
%% TODO: compute blended image
imret = im1;
for i =1: length(x)
    imret(x(i),y(i),1) = X1(i);
    imret(x(i),y(i),2) = X2(i);
    imret(x(i),y(i),3) = X3(i);
end

