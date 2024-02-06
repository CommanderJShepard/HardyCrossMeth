clc
clear

u = input('What is the kinematic visocisity of the fluid?: ');
p = input ('What is the density of the fluid? (in slugs/ft^3): ');
E = input('What is the surface roughness of the pipes (Epsilon; in feet)?: ');
loops = input( 'How many loops are present in the system?: ');
Q = zeros(1,loops); %Build # of rows for flow
L = zeros(1,loops);%Build # of rows for length
d = zeros(1,loops); %Build # of rows for diameter


for i = 1:loops
    promptString = sprintf('How many lines are in loop %d ?: ', i);
    input2 = input(promptString);
    Q=zeros(i,input2);
    L=zeros(i,input2);
    d=zeros(i,input2);
    Re=zeros(i, input2);
    Ve=zeros(i, input2);
    f=zeros(i, input2);
    clindicator=zeros(i, input2);
    press=zeros(i,input2);
end

for j = 1:height(Q)
    for k = 1:width(Q)
        promptString = sprintf('What is the guessed flowrate of loop %d, line %d? (In CFS. Enter 0 if line does not exist): ', j, k);
        input2 = input(promptString);
        Q(j,k) = input2;
    end
end

for j = 1:height(Q)
    for k = 1:width(Q)
        promptString = sprintf('What is the length of loop %d, line %d? (in feet): ', j, k);
        input2 = input(promptString);
        L(j,k) = input2;
    end
end

for j = 1:height(Q)
    for k = 1:width(Q)
        promptString = sprintf('What is the diameter of loop %d, line %d? (in FEET): ', j, k);
        input2 = input(promptString);
        d(j,k) = input2;
    end
end

clnumber = input('How many common lines are there in the system?');
a = 1;
for j = 1:clnumber
        
       promptString = sprintf('What is the loop number where the common line occurs?: ');
       input1 = input(promptString);
       promptString2 = sprintf('What is the line number of that loop?: ');
       input2 = input(promptString2);
       promptString3 = sprintf('What is the loop number that shares a line with loop %d?: ', input1);
       input3 = input(promptString3);
       promptString4 = sprintf('What is the line number of that loop?: ');
       input4 = input(promptString4);
       clindicator(input1,input2) = a;
       clindicator(input3,input4) = a;
       a = a+1;
end


g=32.2;
%%%All Data has been collected at this point. Begin iteration process
%% 
e = 5;
iteration=0;
while e>0.0005
    
Ve=Vel(Q,d);
Re=abs(Rey(p, Ve, d, u));
f = FF(E, d, abs(Re));
Al = alpha(L, d);
Head = Hd(Al, f, Q);
%%Convert Head to proper negative values rep direction of flow
[m n]=size(Head);
for i=1:m
    for j=1:n
        if Q(i,j)<0 && Head(i,j)>0
            Head(i,j) = -1*Head(i,j);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Head sum of each loop%
Headsum = sum(Head,2)

%%%%%Find df/dq
dfdq = DFDQ(E, d, Re, Q);
%%%%%%Find dh/dq using df/dq
dhdq = abs(DHDQ(Al, Q, f, dfdq));

dhdqsum = sum(dhdq,2);
deltaqi=-(Headsum)./dhdqsum
    
[m n] = size(Q);
b=1
%Deltaqi for uncommon lines
for i = 1:m
    for j=1:n
        if clindicator (i,j) == 0
            Q(i,j) = Q(i,j) + deltaqi(i,1);
        end
    end
end
%deltaqi for common lines
for b = 1:clnumber
    for i = 1:m
        for j=1:n
            if clindicator (i,j) == b
                for k = 1:m
                    for l = 1:n
                        if clindicator (k,l) == b && clindicator (k,l)~clindicator(i,j);
                            Q(i,j) = Q(i,j) + (deltaqi(i,1)-deltaqi(k,1));
                        end
                    end
                end
            end
        end
    end
end

iteration = iteration + 1;
e = max(abs(deltaqi));
end
fprintf('Convergence factor: %d \n', e)
fprintf('It took %d iterations to fall within the convergence factor. \n', iteration)
[m n]=size(Q);
for i=1:m
    for j=1:n 
fprintf('The flow rate for loop %d, line %d is: %.3f ft^3/s \n', i,j,Q(i,j))
    end
end
for i=1:m
    for j=1:n 
fprintf('The head loss for loop %d, line %d is: %.3f ft \n', i,j,abs(Head(i,j)))
    end
end
%%Specific for Test Case 2 Pressure Calculations
% %%1000kPa = 20885.4 lbf/ft2
% press(1,1) = 20885.4 - 1.94*32.2*abs(Head(1,1))
% press(1,2) = 20885.4 - 1.94*32.2*abs(Head(1,2))
% press(1,3) = press(1,2) - 1.94*32.2*abs(Head(1,3))
% press(1,4) = press(1,1) - 1.94*32.2*abs(Head(1,4))
% press(2,1) = press(1,1) - 1.94*32.2*abs(Head(2,1))
% press(2,2) = press(1,1) - 1.94*32.2*abs(Head(2,2))
% press(2,3) = press(2,2) - 1.94*32.2*abs(Head(2,3))
% press(2,4) = press(2,1) - 1.94*32.2*abs(Head(2,4))
% press(3,1) = press(1,2) - 1.94*32.2*abs(Head(3,1))
% press(3,2) = press(1,2) - 1.94*32.2*abs(Head(3,2))
% press(3,3) = press(3,2) - 1.94*32.2*abs(Head(3,3))
% press(3,4) = press(3,1) - 1.94*32.2*abs(Head(3,4))
% press(4,1) = press(1,3) - 1.94*32.2*abs(Head(4,1))
% press(4,2) = press(1,4) - 1.94*32.2*abs(Head(4,2))
% press(4,3) = press(3,3) - 1.94*32.2*abs(Head(4,3))
% press(4,4) = press(2,4) - 1.94*32.2*abs(Head(4,4))

%Specific For Test Case 1. Head Loss Dist.
% NodeA = 120 - abs(Head(1,4))
% NodeH = NodeA - abs(Head(1,3))
% NodeE = 120 - abs(Head(1,1))
% NodeC = 120 - abs(Head(2,1))
% NodeD = NodeE - abs(Head(2,3))
% NodeG = NodeH - abs(Head(3,3))
% NodeF = NodeE - abs(Head(3,1))


function [ALPHA] = alpha(L, d)
ALPHA = (8.*L)./(pi^2*32.2.*d.^5);
end

function [HEAD] = Hd(Al, f, Q)
HEAD = Al.*f.*(Q.^2);
end

function [dfdq] = DFDQ(E, d, Re, Q)
dfdq = (13.69*(E./(3.7.*d)+(5.74./(Re.^0.9))).^-1)./((Re.^0.9).*Q.*(log(E./(3.7.*d)+5.74./(Re.^0.9))).^3);
end

function [dhdq] = DHDQ(Al, Q, f, dfdq)
dhdq = 2*Al.*Q.*f + Al.*(Q.^2).*dfdq;
end

function velocity = Vel(Q,d)
velocity = Q./((d./2).^2*pi);
end

function Reynolds = Rey(p, Ve, d, u) 
Reynolds = (p.*Ve.*d)./u;
end

function friction = FF(E, d, Re)
friction = (1.325)./((log(((E./d)/3.7)+(5.74./(Re.^0.9)))).^2);
end






