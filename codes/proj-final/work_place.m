clear
clc
close all


%whole_rectangle(1,2,1); %you can use i to determine the mesh sparsebess, the number of node on one direction is 20*i

%calculate the error in whole plate
% error_data = zeros(4,3);
% for i = 1:4
%     error_data(i, :) = whole_rectangle(i, 1, 1); 
% end
% 
% figure;
% plot(error_data(:,3),error_data(:,1)); %the H1 error
% hold on
% plot(error_data(:,3),error_data(:,2));
% legend('e in H1','e in L2');
% title('error convergence rate')


with_hole(1,4,2); 

%calculate the error in plate with hole
% error_data = zeros(4,3);
% for i = 1:4
%     error_data(i, :) = with_hole(i, 1, 1); % the max value of i is 4
% end
% 
% figure;
% plot(error_data(:,3),error_data(:,1)); %the H1 error
% legend('e in H1');
% title('error convergence rate')


%with_hole_last(3, 1); %you can only input 2, 3, 4 at the first input place



