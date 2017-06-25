%OH_eigenvalues(B,E,beta)

figure
range = 0:500:500000;

Aseq = zeros(8,8,length(range));
for i = 1:length(range)
    Aseq(:,:,i) = OH_Ham_Lab_Fixed(0,0,0,range(i),0,0)/(6.626e-28);
end
[Vseq,Dseq] = eigenshuffle(Aseq);

plot(range/100,Dseq'-733.7)

%plot(range,cell2mat(arrayfun(@(B) OH_eigenvalues_HF(B,1000,pi/10),range,'UniformOutput',0))')
% %%
% figure
% all = cell2mat(arrayfun(@(B) OH_eigenvalues_HF(B,350,pi/2),range,'UniformOutput',0))';
% plot(range,all(:,1)-all(:,3),range,all(:,3)-all(:,5),range,all(:,3)-all(:,6))
% 
% % figure
% % subplot(3,1,1)
% % range = 1400:0.1:1600;
% % plot(range,cell2mat(arrayfun(@(B) OH_eigenvalues(B,350,pi),range,'UniformOutput',0))')
% % ylim([-0.1,0.1])
% % subplot(3,1,2)
% % plot(range,cell2mat(arrayfun(@(B) OH_eigenvalues(B,350,pi/10),range,'UniformOutput',0))')
% % ylim([-0.1,0.1])
% % subplot(3,1,3)
% % plot(range,cell2mat(arrayfun(@(B) OH_eigenvalues(B,350,pi/2),range,'UniformOutput',0))')
% % ylim([-0.1,0.1])
% 
% %okay first we have to find the location of the min gap
% gap = @(B,beta) min(diff(OH_eigenvalues(B,100,beta)));
% mingap = @(beta) min(cell2mat(arrayfun(gap,1489:.01:1493,repmat([beta],1,401),'UniformOutput',0)));
% 
% figure
% subplot(1,2,1)
% range = 0:pi/100:pi/2;
% plot(180*range/pi,1000*cell2mat(arrayfun(@(beta) mingap(beta),range,'UniformOutput',0))')
% title('Gap at X_{1/2} Crossing','FontSize',14)
% xlabel('EB angle (deg)','FontSize',12)
% ylabel('Energy (MHz)','FontSize',12)
% 
% 
% subplot(1,2,2)
% range = 0:pi/10000:pi/400;
% plot(180*range,cell2mat(arrayfun(@(beta) 1-exp(-2*pi*mingap(beta)^2/(8.91*10^-7)),range,'UniformOutput',0))')
% title('Adiabatic Probability v Angle','FontSize',14)
% xlabel('EB angle (deg)','FontSize',12)
% ylabel('Hopping Probability','FontSize',12)
