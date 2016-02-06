
%This code is to find the left singular vectors from the perturbation of
%column augmented strain matrix based on the Stewart iterative perturbation
%scheme.

%Here we can mention the singular vector to be perturbed in this code so
%that it will swap accordingly and decompose the matrix and perturb
%accordingly.

%E1Matrix,E2Matrix,E6Matrix ----Strain matrices
%A---Row augmented strain matrix
%m,n---Dimension of A
%U---Matrix whose columns are left singular vectors of A
%S---Diagonal Matrix whose elements are singular values of A. 
%V---Matrix whose columns are right singular vectors of A (Since we take
%transpose this will be our left singular vectors)
clc;
clear all;
%% Obtaining the number of realisations
prompt1='Enter the number of realisations for which you want to run this code : ';
numCopies = input(prompt1);

%%Loading the data and performing SVD on row-augmented matrix
set(0,'DefaultFigureWindowStyle','docked');
load('Eps1.mat');
load('Eps2.mat');
load('Eps6.mat');
A = [E1Matrix, E2Matrix, E6Matrix];
[m,n]=size(A);
iOrder=1;
if m<n
    A=A';
end
[m,n]=size(A);
[U,S,V]=svd(A);

%%Swapping the rows and columns of singular values,singular vectors based on the eigen function input
prompt2='Enter the left eigen function to be used : ';
k=input(prompt2);
swap=S(1,1);
S(1,1)=S(k,k);
S(k,k)=swap;

sw=U(:,1);
U(:,1)=U(:,k);
U(:,k)=sw;

sw=V(:,1);
V(:,1)=V(:,k);
V(:,k)=sw;

%%Stewart perturbation scheme implementation
m_1 = 1;  % number of SVs of interest
m_2 = n - m_1;
% U
U1 = U(:, 1:m_1);
U2 = U(:, m_1+1:n);
U3 = U(:, n+1:m);
% S
S1 = S(1:m_1, 1:m_1);
S2 = S(m_1+1:n , m_1+1:n);
%
V1 = V(:,1:m_1);
V2 = V(:,m_1+1:n);
%
u1 = U(:,1);
v1 = V(:,1);

sigma1 = S1(1:1);
sigNoiseRatio = [25 20 15 10 5  0 -5 -10 -15 -20 -25];
Doc_leftEV = zeros(size(V,2),size(sigNoiseRatio,2),numCopies);
Doc_leftEV_true = zeros(size(V,2),size(sigNoiseRatio,2),numCopies);

for counter_copies = 1 : 1 :numCopies
	
	counter=1;
	for counter=1:size(sigNoiseRatio,2)
		[counter_copies,counter]
		 
        
        [Atilde,E]=add_wgn(A,sigNoiseRatio(counter)); 
		[Utrue,Strue,Vtrue] = svd(Atilde);
        
		% initialize iterates
		x0 = v1;
		gamma0 = sqrt(norm(sigma1 + u1'*E*v1)^2 + norm(U3'*E*v1)^2);
		% gamma1 = gamma0;
		y0 = u1*gamma0;

		%% get higher order perturbation estimates 
		u1 = U(:,1);
		v1 = V(:,1);
		xvec = zeros(n,iOrder+1);
		yvec = zeros(m,iOrder+1);
		uvec = zeros(m,iOrder+1);
		gammavec = zeros(iOrder+1,1);
	   
		%
		xvec(:,1) = x0;
		yvec(:,1) = y0;
		gammavec(1) = gamma0;
		uvec(:,1)=y0/gamma0;
		%
		% figure; 
		%
		for i = 1:iOrder
			%
			j = i + 1;
			x = xvec(:,i);% right
			y = yvec(:,i);% left
			gamma = gammavec(i);
			%
			temp = gamma^2*eye(size(S2)) - S2*S2;
			temp1 = temp\(V2'*E'*y + S2*U2'*E*x);
			VHx = [sqrt(1 - norm(V2'*x)^2); temp1];
			xvec(:,j) = V*VHx;
			%
			temp1 = sqrt(1 - norm(V2'*x)^2)*sigma1 + u1'*E*x;
			temp2 = temp\(S2*V2'*E'*y + gamma^2*U2'*E*x);
			temp3 = U3'*E*x;
			UHy = [temp1; temp2; temp3];
			yvec(:,j) = U*UHy;
			%
			gammavec(j) = norm(Atilde*xvec(:,j));   
			uvec(:,j) = yvec(:,j)/gammavec(j);

		end
		
		
		

		
		Doc_leftEV(:,counter,counter_copies) = xvec(:,2); %Since A=A' we are taking xvec as v.
        Doc_leftEV_true(:,counter,counter_copies) = Vtrue(:,k); %First column of U here is the first right eigen vector of the untransposed row-augmented matrix
		
	end
end
save('userip_copies_11NLevels_1.mat','V','sigNoiseRatio','Doc_leftEV','Doc_leftEV_true');

