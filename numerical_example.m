% NUMERICAL SCRIPT FOR: A mathematical model for elasticity using calculus on discrete manifolds

%% KEY
%A=Incidence matrix
%a_{i,j}=the element of A. If i->j=1, if j->i=-1, else 0.
%D=Initial positions of nodes
%b=Initial edges
%X=Final positions of nodes
%y=Final edges
%B=External forces at nodes
%F=Forces in primary bonds (spring energy due to change of edge length)
%K=Material coefficient vector, given for each bond
%a0i=Lower bound on the change of primary edge length
%a1i=Upper bound on the change of primary edge length
%a2i=Lower bound on the change of primary edge length
%a3i=Upper bound on the change of primary edge length

%% SETTING UP D1, D2, A1, A2 
D1=[0,0,1;
    1,0,0;
    0,-1,0;
    -1,0,0;
    0,1,0;
    0,0,-1];
A1= [1,0,0,-1,0,0;
     1,0,0,0,-1,0;
     1,-1,0,0,0,0;
     1,0,-1,0,0,0;
     0,0,1,-1,0,0;
     0,0,0,1,-1,0;
     0,-1,0,0,1,0;
     0,1,-1,0,0,0;
     0,0,0,-1,0,1;
     0,0,0,0,-1,1;
     0,-1,0,0,0,1;
     0,0,-1,0,0,1];
A2 =[1,-1,0,0,0,0,0,0,0,0,0,0;
     0,1,-1,0,0,0,0,0,0,0,0,0;
     0,0,1,-1,0,0,0,0,0,0,0,0;
     -1,0,0,1,0,0,0,0,0,0,0,0;
     1,0,0,0,-1,0,0,0,0,0,0,0;
     1,0,0,0,0,-1,0,0,0,0,0,0;
     0,1,0,0,0,-1,0,0,0,0,0,0;
     0,1,0,0,0,0,-1,0,0,0,0,0;
     0,0,1,0,0,0,-1,0,0,0,0,0;
     0,0,1,0,0,0,0,-1,0,0,0,0;
     0,0,0,1,0,0,0,-1,0,0,0,0;
     0,0,0,1,-1,0,0,0,0,0,0,0;
     0,0,0,0,1,-1,0,0,0,0,0,0;
     0,0,0,0,0,1,-1,0,0,0,0,0;
     0,0,0,0,0,0,1,-1,0,0,0,0;
     0,0,0,0,-1,0,0,1,0,0,0,0;
     0,0,0,0,1,0,0,0,0,0,0,-1;
     0,0,0,0,1,0,0,0,-1,0,0,0;
     0,0,0,0,0,1,0,0,-1,0,0,0;
     0,0,0,0,0,1,0,0,0,-1,0,0;
     0,0,0,0,0,0,1,0,0,-1,0,0;
     0,0,0,0,0,0,1,0,0,0,-1,0;
     0,0,0,0,0,0,0,1,0,0,-1,0;
     0,0,0,0,0,0,0,1,0,0,0,-1;
     0,0,0,0,0,0,0,0,1,-1,0,0;
     0,0,0,0,0,0,0,0,0,1,-1,0;
     0,0,0,0,0,0,0,0,0,0,1,-1;
     0,0,0,0,0,0,0,0,-1,0,0,1;
     1,0,-1,0,0,0,0,0,0,0,0,0;
     0,1,0,-1,0,0,0,0,0,0,0,0;
     0,0,0,0,0,0,0,0,1,0,-1,0;
     0,0,0,0,0,0,0,0,0,1,0,-1;
     1,0,0,0,0,0,0,0,-1,0,0,0;
     0,1,0,0,0,0,0,0,0,-1,0,0;
     0,0,1,0,0,0,0,0,0,0,-1,0;
     0,0,0,1,0,0,0,0,0,0,0,-1];
mid=0.5*abs(A1);      
D2=[mid*D1(:,1),mid*D1(:,2),mid*D1(:,3)];

%% CALCULATING b1 & b2 from A1, D1 & A2,D2
b1=[A1*D1(:,1),A1*D1(:,2),A1*D1(:,3)];
b2=[A2*D2(:,1),A2*D2(:,2),A2*D2(:,3)];

%% CALCULATING LENGTH OF EDGES
lengthb1=sqrt(b1(:,1).^2+b1(:,2).^2+b1(:,3).^2);
lengthb2=sqrt(b2(:,1).^2+b2(:,2).^2+b2(:,3).^2);

%% ASSIGNING SOME VALUE FOR K1 and K2 (trial and error to find this value)
K1=.0001*lengthb1;
K2=.0001*lengthb2;
%% ASSIGNING SOME VALUE FOR a0-a3 proportional to edge length
a0=.15*lengthb1;
a1=5*lengthb1;
a2=.15*lengthb2;
a3=5*lengthb2;
%% Defining m1 and m2 the number of inner & outer edges
m1=length(b1);m2=length(b2);
%% SETTING UP Ktilda and Atilda
KtildaVector=0*K1;
%% Assuming a0i>0 and a2i>0
for i=1:m1
    for j=1:m2
        if A2(j,i)^2>0
            KtildaVector(i)=KtildaVector(i)+0.5*((K2(j)*a2(j)^2)/(2*a1(i)*(a1(i)+abs(b1(i))))+(K2(j)*a3(j)^2)/(2*a0(i)*(a0(i)+abs(b1(i)))));
        end
    end
    KtildaVector(i)=K1(i)+KtildaVector(i)-0.5*(abs(b1(i))*K1(i)/(a0(i)+abs(b1(i)))+abs(b1(i))*K1(i)/(a1(i)+abs(b1(i))));
end
Ktilda=gallery('tridiag',zeros(1,m1-1),KtildaVector,zeros(1,m1-1));

Atilda=A1'*Ktilda*A1;

%% LINEAR ALGEBRA SECTION
% SETTING UP THE 4 ATILDA MATRICES
P=5;%Here we say 5 nodes move, one node remains fixed.
Q=length(D1)-P;

Atilda11=Atilda(1:P,1:P);
Atilda12=Atilda(1:P,P+1:P+Q);
Atilda21=Atilda(P+1:P+Q,1:P);
Atilda22=Atilda(P+1:P+Q,P+1:P+Q);


XQ=D1(P+1:P+Q,:);
% DEFINING THE FORCES ON ALL THE NODES THAT MOVE
BP=zeros(P,3);
for i=1:P
    BP(i,:)=D1(i,:)/7*(1.11)^i; %Forces on nodes away from origin so a0&a2>0
end

% LINEAR ALGEBRA CALCULATIONS
XP=Atilda11\(BP-Atilda12*XQ);
BQ=Atilda21*XP+Atilda22*XQ;
X1=[XP;XQ];
X2=[mid*X1(:,1),mid*X1(:,2),mid*X1(:,3)];
y1=[A1*X1(:,1),A1*X1(:,2),A1*X1(:,3)];
lengthy1=sqrt(y1(:,1).^2+y1(:,2).^2+y1(:,3).^2);

y2=[A2*mid*X1(:,1),A2*mid*X1(:,2),A2*mid*X1(:,3)];
lengthy2=sqrt(y2(:,1).^2+y2(:,2).^2+y2(:,3).^2);

%% PLOTTING ORIGINAL NODES AND EDGES
figure

% PLOT OUTER AND INNER ORIGINAL NODES
scatter3(D1(:,1),D1(:,2),D1(:,3),'b','filled')
hold on
scatter3(D2(:,1),D2(:,2),D2(:,3),'r','filled')
% PLOT EDGES OF ORIGINAL OUTER NODES
for i=1:12;
    plot3([D1(sum([1:6].*(A1(i,:)==-1)),1);D1(sum([1:6].*(A1(i,:)==1)),1)],[D1(sum([1:6].*(A1(i,:)==-1)),2);D1(sum([1:6].*(A1(i,:)==1)),2)],[D1(sum([1:6].*(A1(i,:)==-1)),3);D1(sum([1:6].*(A1(i,:)==1)),3)],'b','linewidth',2);
end
for i=1:length(A2);
    plot3([D2(sum([1:12].*(A2(i,:)==-1)),1);D2(sum([1:12].*(A2(i,:)==1)),1)],[D2(sum([1:12].*(A2(i,:)==-1)),2);D2(sum([1:12].*(A2(i,:)==1)),2)],[D2(sum([1:12].*(A2(i,:)==-1)),3);D2(sum([1:12].*(A2(i,:)==1)),3)],'r','linewidth',2);
end
axis square
view(20,20)
title('Initial nodes and edges','Interpreter','latex', 'FontSize',16)
grid on
%% PLOTTING DISPLACED NODES AND EDGES
figure
hold on
% PLOT OUTER AND INNER DISPLACED NODES
scatter3(X1(:,1),X1(:,2),X1(:,3),'b','filled')
hold on
scatter3(X2(:,1),X2(:,2),X2(:,3),'r','filled')
% PLOT EDGES OF DISPLACED OUTER NODES
for i=1:12;
    plot3([X1(sum([1:6].*(A1(i,:)==-1)),1);X1(sum([1:6].*(A1(i,:)==1)),1)],[X1(sum([1:6].*(A1(i,:)==-1)),2);X1(sum([1:6].*(A1(i,:)==1)),2)],[X1(sum([1:6].*(A1(i,:)==-1)),3);X1(sum([1:6].*(A1(i,:)==1)),3)],'b','linewidth',2);
end
% PLOT EDGES OF DISPLACED INNERER NODES
for i=1:length(A2);
    plot3([X2(sum([1:12].*(A2(i,:)==-1)),1);X2(sum([1:12].*(A2(i,:)==1)),1)],[X2(sum([1:12].*(A2(i,:)==-1)),2);X2(sum([1:12].*(A2(i,:)==1)),2)],[X2(sum([1:12].*(A2(i,:)==-1)),3);X2(sum([1:12].*(A2(i,:)==1)),3)],'r','linewidth',2);
end
axis square
view(20,20)
title('Final nodes and edges','Interpreter','latex', 'FontSize',16)


grid on
%% MAKING SURE 'a_n,i' PARAMETER VALUES ARE CORRECTLY DEFINED
a0FailMessage='Warning: a_oi>|y_i|-|b_i|';
a1FailMessage='Warning: a_1i<|y_i|-|b_i|';
a2FailMessage='Warning: a_2i>|y_i|-|b_i|';
a3FailMessage='Warning: a_3i<|y_i|-|b_i|';
if min(lengthy1-lengthb1-a0)<0;
    error(a0FailMessage)
end
if max(lengthy1-lengthb1-a1)>0;
    error(a1FailMessage)
end
if min(lengthy2-lengthb2-a2)<0;
    error(a2FailMessage)
end
if max(lengthy2-lengthb2-a3)>0;
    error(a3FailMessage)
end