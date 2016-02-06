%userip_copies_11NLevels_2.mat ----'U','sigNoiseRatio','Doc_rightEV','Doc_rightEV_true'
%Doc_right 760*11*30
%A_r terms 1,2,6,rhs_RSV1 are computed here


clear all;
set(0,'DefaultFigureWindowStyle','docked')
load('userip_copies_11NLevels_2.mat');
load('Eps1');
load('Eps2');
load('Eps6');

[r_eps,c_eps] = size(E1Matrix);

term2_actual_6 = 0;
term2_actual_1 = 0;
term2_actual_2 = 0;

for k=1:r_eps
    term2_actual_6 = term2_actual_6 + E6Matrix(k,:)*U(:,1);    
    term2_actual_1 = term2_actual_1 + E1Matrix(k,:)*U(:,1);    
    term2_actual_2 = term2_actual_2 + E2Matrix(k,:)*U(:,1);

end

U_mean = zeros(760,11);
U_mean_true=zeros(760,11);
for j=1:11
    rhs = 0;
    rhs_true = 0;
    rhs_abs=0;
    rhs_abs_true =0;
    summ = zeros(760,1);
    sum_true = zeros(760,1);
    sum_diff_norm = 0;
    sum_diff_norm_true = 0;
    
    % term - 2    
    A_first_colum_sum_6 = 0;    
    A_first_colum_sum_true_6 = 0;
    
    A_first_colum_sum_1 = 0;    
    A_first_colum_sum_true_1 = 0;
    
    A_first_colum_sum_2 = 0;    
    A_first_colum_sum_true_2 = 0;
    % term - 2
    std_A_first_column_sum_6 =0;
    std_A_first_column_sum_6_true =0;
        std_A_first_column_sum_1 =0;
    std_A_first_column_sum_1_true =0;
        std_A_first_column_sum_2 =0;
    std_A_first_column_sum_2_true =0;
    
    for i=1:30
        rhs = rhs + sum(Doc_rightEV(:,j,i));
        rhs_true = rhs_true + sum(Doc_rightEV_true(:,j,i));
        
        rhs_abs = rhs_abs + power(   abs (sum(Doc_rightEV(:,j,i)) - sum(U(:,1))) , 2);        
        rhs_abs_true = rhs_abs + power(    abs (sum(Doc_rightEV_true(:,j,i)) - sum(U(:,1))) , 2);
        
        summ = summ + Doc_rightEV(:,j,i);
        sum_true = sum_true + Doc_rightEV_true(:,j,i);
        sum_diff_norm = sum_diff_norm + norm(Doc_rightEV(:,j,i)-U(:,1),2);        
        sum_diff_norm_true = sum_diff_norm_true + norm(Doc_rightEV_true(:,j,i)-U(:,1),2);
        
        sum_across_rows_6 = 0;
        sum_across_rows_6_true =0;
        
        sum_across_rows_1 = 0;
        sum_across_rows_1_true =0;
        
        sum_across_rows_2 = 0;
        sum_across_rows_2_true =0;
        for k=1:r_eps
            
            sum_across_rows_6 = sum_across_rows_6 + E6Matrix(k,:) * Doc_rightEV(:,j,i);            
            sum_across_rows_6_true = sum_across_rows_6_true + E6Matrix(k,:) * Doc_rightEV_true(:,j,i);
            
            sum_across_rows_1 = sum_across_rows_1 + E1Matrix(k,:) * Doc_rightEV(:,j,i);            
            sum_across_rows_1_true = sum_across_rows_1_true + E1Matrix(k,:) * Doc_rightEV_true(:,j,i);
            
            sum_across_rows_2 = sum_across_rows_2 + E2Matrix(k,:) * Doc_rightEV(:,j,i);            
            sum_across_rows_2_true = sum_across_rows_2_true + E2Matrix(k,:) * Doc_rightEV_true(:,j,i);
        
        end
        A_first_colum_sum_6 = A_first_colum_sum_6 + sum_across_rows_6;        
        A_first_colum_sum_true_6 = A_first_colum_sum_true_6 + sum_across_rows_6_true;    
        
        A_first_colum_sum_1 = A_first_colum_sum_1 + sum_across_rows_1;        
        A_first_colum_sum_true_1 = A_first_colum_sum_true_1 + sum_across_rows_1_true;    
        
        A_first_colum_sum_2 = A_first_colum_sum_2 + sum_across_rows_2;        
        A_first_colum_sum_true_2 = A_first_colum_sum_true_2 + sum_across_rows_2_true; 
        
        std_A_first_column_sum_6 = std_A_first_column_sum_6 + (sum_across_rows_6  - term2_actual_6)^2;  
        std_A_first_column_sum_6_true = std_A_first_column_sum_6_true + (sum_across_rows_6_true  - term2_actual_6)^2;  
        
        std_A_first_column_sum_1 = std_A_first_column_sum_1 + (sum_across_rows_1  - term2_actual_1)^2;  
        std_A_first_column_sum_1_true = std_A_first_column_sum_1_true + (sum_across_rows_1_true  - term2_actual_1)^2;  

        std_A_first_column_sum_2 = std_A_first_column_sum_2 + (sum_across_rows_2  - term2_actual_2)^2;  
        std_A_first_column_sum_2_true = std_A_first_column_sum_2_true + (sum_across_rows_2_true  - term2_actual_2)^2;   
        
    end
    bias_rhs(j) = abs(rhs/30 - sum(U(:,1))) ;
    bias_rhs_true(j) = abs(rhs_true/30 - sum(U(:,1)));
    
    std_rhs(j) = sqrt(rhs_abs/30);
    std_rhs_true(j) = sqrt(rhs_abs_true/30);
    
    bias_term2_6(j) = abs(A_first_colum_sum_6/30  - term2_actual_6);
    bias_term2_6_true(j) = abs(A_first_colum_sum_true_6/30 - term2_actual_6);
    bias_term2_1(j) = abs(A_first_colum_sum_1/30  - term2_actual_1);
    bias_term2_1_true(j) = abs(A_first_colum_sum_true_1/30 - term2_actual_1);
    bias_term2_2(j) = abs(A_first_colum_sum_2/30  - term2_actual_2);
    bias_term2_2_true(j) = abs(A_first_colum_sum_true_2/30 - term2_actual_2);
    
    mean1(j)=rhs/30; %RSV
    mean2(j)=A_first_colum_sum_1/30; %Ar1
    mean3(j)=A_first_colum_sum_2/30; %Ar2
    mean4(j)=A_first_colum_sum_6/30; %Ar6
    
    std_term2_6(j) = sqrt(std_A_first_column_sum_6/30);    
    std_term2_6_true(j) = sqrt(std_A_first_column_sum_6_true/30);
    
    std_term2_1(j) = sqrt(std_A_first_column_sum_1/30);    
    std_term2_1_true(j) = sqrt(std_A_first_column_sum_1_true/30);
    
    std_term2_2(j) = sqrt(std_A_first_column_sum_2/30);    
    std_term2_2_true(j) = sqrt(std_A_first_column_sum_2_true/30);
    
    
    std(j) = sum_diff_norm/30;
    std_true(j) = sum_diff_norm_true / 30;
    U_mean(:,j)=summ/30;
    U_mean_true(:,j)=sum_true/30;
end
bias = zeros(1,11);
bias_true=zeros(1,11);
for i=1:11
    bias(i)=norm(U(:,1)-U_mean(:,i),2);
    bias_true(i) = norm(U(:,1)-U_mean_true(:,i),2);
end
%% Plot bias of the REV1

figure; semilogy(sigNoiseRatio,bias,'b*-');hold on; semilogy(sigNoiseRatio,bias_true,'r*-');
% daspect([1 1 1]);
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Bias (based on the norm of the difference-vector)', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained bias with the true bias', 'interpreter', 'latex', 'FontSize', 22);
legend({'Bias (analytical)','Bias (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');

%% Plot std. dev. of REV1
figure; semilogy(sigNoiseRatio,std,'b*-');hold on; semilogy(sigNoiseRatio,std_true,'r*-');
% daspect([1 1 1]);
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Standard deviation (based on the norm of the difference-vector)', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained Std. Dev. with the true Std. Dev.', 'interpreter', 'latex', 'FontSize', 22);
legend({'Std. Dev. (analytical)','Std. Dev. (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');

%% Plot std. dev. and bias of REV1
figure; semilogy(sigNoiseRatio,bias,'b*-');hold on; semilogy(sigNoiseRatio,bias_true,'r*-'); semilogy(sigNoiseRatio,std,'ks-');hold on; semilogy(sigNoiseRatio,std_true,'gs-');
% daspect([1 1 1]);
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Bias/Standard Deviation', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained bias and std. dev with the true ones', 'interpreter', 'latex', 'FontSize', 22);
legend({'Bias (analytical)','Bias (true)','Std. Dev. (analytical)','Std. Dev. (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');

%% Plot the RHS term (sum of the REV1) - bias

figure; semilogy(sigNoiseRatio(2:11),bias_rhs(2:11),'b*-');hold on; semilogy(sigNoiseRatio(2:11),bias_rhs_true(2:11),'r*-');
% daspect([1 1 1]);
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Bias of summation REV1', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained bias for the RHS term with the true one', 'interpreter', 'latex', 'FontSize', 22);
legend({'Bias - RHS term (analytical)','Bias - RHS (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');

%% Plot the RHS term (sum of the REV1) - std. dev.

% figure; semilogy(sigNoiseRatio(2:11),std_rhs(2:11),'b*-');hold on; semilogy(sigNoiseRatio(2:11),std_rhs_true(2:11),'r*-');
figure; semilogy(sigNoiseRatio,std_rhs,'b*-');hold on; semilogy(sigNoiseRatio,std_rhs_true,'r*-');

% daspect([1 1 1]);
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Standard Deviation of summation REV1', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained std. dev. for the RHS term with the true one', 'interpreter', 'latex', 'FontSize', 22);
legend({'Std. Dev. - RHS term (analytical)','Std. Dev. - RHS (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');


%% Plot the RHS term - std. dev. and bias

figure; semilogy(sigNoiseRatio(2:11),std_rhs(2:11),'b*-');hold on; semilogy(sigNoiseRatio(2:11),std_rhs_true(2:11),'rd-');  
semilogy(sigNoiseRatio(2:11),bias_rhs(2:11),'ms-');hold on; semilogy(sigNoiseRatio(2:11),bias_rhs_true(2:11),'bs-');

% figure; semilogy(sigNoiseRatio,std_rhs,'b*-');hold on; semilogy(sigNoiseRatio,std_rhs_true,'r*-');

% daspect([1 1 1]);
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Standard Deviation/Bias of the RHS term', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained std. dev. and bias for the RHS term with the true valeus', 'interpreter', 'latex', 'FontSize', 22);
legend({'Std. Dev. - RHS term (analytical)','Std. Dev. - RHS (true)','Bias - RHS term (analytical)','Bias - RHS (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');








%% Plot the bias in term2 - 6
figure; semilogy(sigNoiseRatio,bias_term2_6,'b*-');hold on; semilogy(sigNoiseRatio,bias_term2_6_true,'r*-');
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Bias of the term-2', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained bias for the term-2 with the true value', 'interpreter', 'latex', 'FontSize', 22);
legend({'Bias - term-2 (analytical)','Bias - term-2 (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');



%% Plot the bias in term2 - 1
figure; semilogy(sigNoiseRatio,bias_term2_1,'b*-');hold on; semilogy(sigNoiseRatio,bias_term2_1_true,'r*-');
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Bias of term-2', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained bias for the term-2 with the true value', 'interpreter', 'latex', 'FontSize', 22);
legend({'Bias - term-2 (analytical)','Bias - term-2 (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');

%% Plot the bias in term2 - 2
figure; semilogy(sigNoiseRatio,bias_term2_2,'b*-');hold on; semilogy(sigNoiseRatio,bias_term2_2_true,'r*-');
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Bias of term-2', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained bias for the term-2 with the true value', 'interpreter', 'latex', 'FontSize', 22);
legend({'Bias - term-2 (analytical)','Bias - term-2 (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');







%% Plot the std in term2 - 6
figure; semilogy(sigNoiseRatio,std_term2_6,'b*-');hold on; semilogy(sigNoiseRatio,std_term2_6_true,'r*-');
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Std. Dev. of the term-2', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained std. dev. for the term-2 with the true value', 'interpreter', 'latex', 'FontSize', 22);
legend({'std. dev. - term-2 (analytical)','Std. dev. - term-2 (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');



%% Plot the std in term2 - 1
figure; semilogy(sigNoiseRatio,std_term2_1,'b*-');hold on; semilogy(sigNoiseRatio,std_term2_1_true,'r*-');
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Std. Dev. of term-2', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained std. dev. for the term-2 with the true value', 'interpreter', 'latex', 'FontSize', 22);
legend({'std. dev. - term-2 (analytical)','Std. dev. - term-2 (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');


%% Plot the std in term2 - 2
figure; semilogy(sigNoiseRatio,std_term2_2,'b*-');hold on; semilogy(sigNoiseRatio,std_term2_2_true,'r*-');
set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Std. Dev. of term-2', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained std. dev. for the term-2 with the true value', 'interpreter', 'latex', 'FontSize', 22);
legend({'std. dev. - term-2 (analytical)','Std. dev. - term-2 (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');





%% Plot the std in term2 - 1,2 and 6
figure;semilogy(sigNoiseRatio,std_term2_1,'b*-');hold on; semilogy(sigNoiseRatio,std_term2_1_true,'r*-');
semilogy(sigNoiseRatio,std_term2_2,'ks-');hold on; semilogy(sigNoiseRatio,std_term2_2_true,'ms-');
semilogy(sigNoiseRatio,std_term2_6,'gd-');hold on; semilogy(sigNoiseRatio,std_term2_6_true,'y*-');

set(gca,'FontSize',20);
xlabel('SNR (Signal-to-noise ratio)', 'interpreter', 'latex', 'FontSize', 22);
ylabel('Std. Dev. of the $\sum\mathbf{A^r}$ term ', 'interpreter', 'latex', 'FontSize', 22);
% zlabel('zlabel', 'interpreter', 'latex', 'FontSize', 22);
title('Comparing the analytically obtained std. dev. for the summation-$\mathbf{A^r}$ term with the true value', 'interpreter', 'latex', 'FontSize', 22);
legend({'$\sum\mathbf{A^r}$ term from $\varepsilon_1$ (analytical)','$\sum\mathbf{A^r}$ term from $\varepsilon_1$ (true)','$\sum\mathbf{A^r}$ term from $\varepsilon_2$ (analytical)','$\sum\mathbf{A^r}$ term from $\varepsilon_2$ (true)','$\sum\mathbf{A^r}$ term from $\varepsilon_6$ (analytical)','$\sum\mathbf{A^r}$ term from $\varepsilon_6$ (true)'},'FontSize', 20);
h = legend;
set(h, 'interpreter', 'latex');
