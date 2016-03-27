clc;clear;
load('seismo.dat');
plot(seismo(:,1),seismo(:,2),'k','LineWidth',2.0)
xlabel('t/s','FontSize',16,'Color','r')
ylabel('v/m/s','FontSize',16,'Color','r')
legend('fd1d_stg4th')