% Activity for the Data Analysis Course
% Week 1: Introduction to Matlab

% Student name:
stname = 'Bob Marley';
% Student codigo:
codigo = '0000000';
disp(['This is the work of ' stname ' with codigo ' codigo])
disp('')
disp('Activity 1: Introduction to Matlab')

%% Commands
disp('I - Commands')
% 1 - Check the help of a function of your choice using the command 'help'

% 2 - Create a variable 'x' with a range from 0 to 10 with increments of
% 0.01

% 3 - Check the class of variable 'x' is and its size with the command
% 'whos'

disp('The variable x has a size of XX and its class is XX')

% 4 - What is the class of the variable 'stname'
disp('The class of the variable stname is XX')

% 5 - Can I combine variables of different classes usually?
disp('Can I combine variables of difference classes: Yes/No')

%% Matrix operations
disp('II - Matrix operations')
% 1 - Create a ROW vector 'r' ranging from 0 to 5 with increments of 1

% 2 - Create a COLUMN vector 'c' ranging from 5 to 0 with increments of 1

% 3 - Try to do c*c. Why this doesn't work?
disp('This doesn t work because: ...')

% 4 - How can you solve this problem if you want to multiply the 2 vectors
% value by value?

% 5 - How can you solve this problem if you want to multiply the 2 vectors
% and sum the values at the same time with only one math operation?

% 6 - Create a ROW vector 'r1' manually with these values: 4 10 56

% 7 - Create a COLUMN vector 'c1' manually with these values: 4 10 56

% 8 - In Matlab you can't multiply arrays without checking their
% dimensions, but multiplying with a scalar is simpler. Multiply variable
% 'r1' by 10.

% 9 - What is the function to know the dimensions of an array?
disp('The function to know array dimensions is: ...')

%% Simple plotting
disp('III - Simple plotting')
clear % Clear everything
% 1 - Create a sine signal using the 'x' variable you created earlier and
% the function 'sin'. Call this variable 'sigs'.

% 2 - Plot the variable 'sig' as a function of 'x'. Don't forget to put
% labels on each axis and a title to the figure.
figure;

% 3 - Create a sine signal using the 'x' variable you created earlier and
% the function 'cos'. Call this variable 'sigc'.

% 4 - Create a exponential signal using the 'x' variable you created earlier and
% the function 'exp' (exp(-x)). Call this variable 'sige'.

% 5 - Do a 'subplot' following these indications:
% + The first subplot spans two columns and only one row and plot the
% following signal in there: sig4 = sige.*sigs;
% Then use the 'hold on' command to plot the exponential signal 'sige'
% + The 2nd subplot is on the 1st column and 2nd row and plot the sigc in
% there
% + The 3rd subplot is on the 2nd column and 2nd row and plot the sigs in
% there
figure;

% If you have some time: change the color and the linewidth of the last
% plot you did using the following properties: 'LineWidth', 'k'

% 6 - Click on the plot containing both variables 'sig4' and 'sige'. Using
% the command 'legend', put a legend to this plot.

% 7 - Plotting and handling a 2D plot in Matlab:
% + Load the clown image (it creates the 'X' variable): load clown
% + Plot this image with the imagesc command: imagesc(1:320,1:200,X)
% + Change the colormap to another one (jet is default): colormap()
% + Show the colorbar of the plot with the command: colorbar
% + Play with the caxis command to change the color bounds
figure; 

%% Operations with loops and logical operators
disp('IV - Operations, loops and logical operators'); clear
% 1 - Do a simple if elseif else end sequence
% + Generate a random integer with the 'randi' function between 1 and 10
% + Make if elseif else end sequence so that
%   - It says 'My random number is lower than 4' if the random number is
%   lower than 4 (strictly lower)
%   - It says 'My random number is between 4 and 8' if the random number is
%   between 4 and 8 (including 4 and 8)
%   - It says 'My random number is greater than 8' if the random number is
%   greater than 8 (strictly greater)


% 2 - Do a 'for' loop following these instructions
% + Do a first 'for' loop with indices from 1 to 10 with increments of 1
% + Do a nested 'for' loop with indices from 1 to 15 with increments of 1
% + Within the nested loop, have a variable 'lop1' calculating the
% multiplication between both indices and assigning this value to the 'lop1'
% variable as 'lop1(ind1,ind2)'.


% 3 - Do another 'for' loop but this time the increment isn't one as
% + Do a 'for' loop with indices from 10 to 100 with increments of 3
% + Calculate a variable 'lop2' calculating the following multiplication
% 'lop2(ind1) = indices^2'. Use another variable 'k' as indices 'ind1' 
% since you can't use the indices of the 'for' loop this time.


%% Working with different classes of variables
disp('V - Working with different classes of variables'); clear
% 1 - Make a structure array 'student' with the following fields:
% + Add a field called 'name' with your name in it (first and last names)
% + Add a field called 'age' with your age in it
% + Add a field called 'streets' with the 5 street numbers


student
% 2 - Call the third value of your 'streets' field in your 'student'
% structure array

% 3 - Call your first name from the 'name' field in your 'student'
% structure array


% 4 - Make a cell array called 'stucell' with the same fields as before 
% with each field in a different cell


stucell

