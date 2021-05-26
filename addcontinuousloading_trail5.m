% General Code For all 1D problems
% Ayush Kumar Dixit, FEM , MIN 324
NoElements=input('No of elements:')
x=sym('x');
L=input('Enter the length:');
A=input('Enter the area :');
E=input('Enter the youngs modulus :');
element_type = input ('Enter 1 for linear element and 2 for quadratic element: ');
switch element_type
case 1
  % Calculation of Number of Nodes
 nodeCoordinates=linspace(0,L,NoElements+1);
 xx=nodeCoordinates;
ii=1:NoElements;
elementNodes=zeros(NoElements,2);
elementNodes(:,1)=ii;
elementNodes(:,2)=ii+1;
NoNodes=size(nodeCoordinates,2);
% Define displacements,force and stiffness matrices
displacements=zeros(NoNodes,1);
force=zeros(NoNodes,1);
stiffness=zeros(NoNodes);
force=input('Enter the force applied in vector form  at each node in [] :')
l=L/NoElements;
B=[-1/l  1/l];
  b=input('Enter the value of body force per unit length: ');
  % Body force matrix calculation
BodyForceGlobal=zeros(NoNodes,1);
for z=1:NoElements ;
    N=[z-(x/l) (x/l)+1-z];
    BodyForce=zeros(2,1);
    M1=b*N';
fun = @(x)M1;
BodyForce(:,1)=int(fun,x,(z-1)*l,z*l);
disp(BodyForce);
    BodyForceGlobal(z:z+1,1)=BodyForceGlobal(z:z+1,1)+ BodyForce(1:2,1);
end
% Stiffness matrix calculations
for e=1:NoElements;
KLocal=zeros(2);
KL=B'*A*E*B;
fun = @(x)KL;
KLocal=int(fun,x,(e-1)*l,e*l);
    % elementDof: element degrees of freedom (Dof)
stiffness(e:e+1,e:e+1)=stiffness(e:e+1,e:e+1)+KLocal;
end
% define activedof since we can not directly put f=kx
NoFixedNodes=input('Enter 1 for fixed-free and 2 for fixed-fixed: ')
switch NoFixedNodes
    case 1 
fixedDof=[1];
    case 2
        fixedDof=[1;NoNodes];
end
Netforce=zeros(NoNodes,1);
Netforce(:,1)=force(:,1)+BodyForceGlobal(:,1);
activeDof=setdiff([1:NoNodes]',fixedDof);
U=stiffness(activeDof,activeDof)\Netforce(activeDof,1);
displacements(activeDof)=U;
fprintf('Displacements %6.12f\n',double(displacements));
forcewithreaction=zeros(NoNodes,1);
forcewithreaction=stiffness*displacements;
fprintf('Net Force Matrix %6.12f\n',double(forcewithreaction));
% Now Strain and stress
% e=input('Insert the element for which you want the strain and stress')
strain=zeros(NoElements,1);
stress=zeros(NoElements,1);
for e=1:NoElements;  
    strain(e,1)= strain(e,1)+ B*displacements(e:e+1,1);
    stress= A*E*strain;
end
fprintf('Strain %6.12f\n',double(strain));
% Add interpolation
xcoordinate=input('Enter the value of length where you want to get the interpolated displacement:');
ElementNoL=ceil(xcoordinate/l);
interpolated_displacement=(1-(xcoordinate+l-ElementNoL*l)/l)*displacements(ElementNoL,1)+((xcoordinate+l-ElementNoL*l)/l)*displacements(ElementNoL+1,1);
fprintf('Value of Interpolated Displacements %6.12f\n',double(interpolated_displacement));;

% Add quadratic element
    case 2
x=sym('x');
l=L/NoElements;
% calculations of Number of Nodes
nodeCoordinates=linspace(0,L,NoElements*2+1);
 xx=nodeCoordinates;
ii=1:NoElements;
elementNodes=zeros(NoElements,3);
elementNodes(:,1)=ii*2-1;
elementNodes(:,2)=ii*2;
elementNodes(:,3)=ii*2+1;
NoNodes=size(nodeCoordinates,2);
% Body Force Calculations
 b=input('Enter the value of body force per unit length: ');
BodyForceGlobal=zeros(NoNodes,1);
for z=1:NoElements ;
    N=[1-3*((x-(z-1)*l)/l)+2*(x^2 +(z-1)^2*l*l -2*x*(z-1)*l)*(1/l^2) 4*((x-(z-1)*l)/l)-4*(x-(z-1)*l)^2*(1/l^2) (-((x-(z-1)*l)/l)+2*(x-(z-1)*l)^2*(1/l^2))];
    BodyForce=zeros(3,1);
    M1=b*N';
fun = @(x)M1;
BodyForce(:,1)=int(fun,x,(z-1)*l,z*l);
fprintf('Body Force %6.12f\n',double(BodyForce));
    BodyForceGlobal((2*z)-1:(2*z)+1,1)=BodyForceGlobal((2*z)-1:(2*z)+1,1)+ BodyForce(1:3,1);
end
disp(BodyForceGlobal);
% Defining displacements, forces and stiffness
displacements=zeros(NoNodes,1);
force=zeros(NoNodes,1);
stiffness=zeros(NoNodes);
force=input('Enter the force applied in vector form  at each node in [] :')
NoFixedNodes=input('Enter the number of nodes that is fixed: ')
switch NoFixedNodes
    case 1 
fixedDof=[1];
    case 2
fixedDof=[1;NoNodes];
end
%Stiffness Matrix Calculations
for z=1:NoElements;
    N=[1-3*((x-(z-1)*l)/l)+2*(x^2 +(z-1)^2*l*l -2*x*(z-1)*l)*(1/l^2) 4*((x-(z-1)*l)/l)-4*(x-(z-1)*l)^2*(1/l^2) (-((x-(z-1)*l)/l)+2*(x-(z-1)*l)^2*(1/l^2))];
B=diff(N);
    KL=B'*A*E*B;
    fun = @(x)KL;
KLocal=int(fun,x,(z-1)*l,z*l);
stiffness(2*z-1:2*z+1,2*z-1:2*z+1)=stiffness(2*z-1:2*z+1,2*z-1:2*z+1)+KLocal;
end
% define activedof since we can not directly put f=kx
activeDof=setdiff([1:NoNodes]',fixedDof);
Netforce=zeros(NoNodes,1);
Netforce(:,1)= force(:,1) + BodyForceGlobal(:,1);
U=stiffness(activeDof,activeDof)\Netforce(activeDof,1);
displacements(activeDof)=U;
forcewithreaction=zeros(NoNodes,1);
forcewithreaction=stiffness*displacements;
fprintf('displacements %6.12f\n'),double(displacements)
% Now Strain and stress
strain=zeros(NoElements,1);
stress=zeros(NoElements,1);
i=input('The number of elements where you want the stress:');
for e=1:i;
 z=input('Enter the element on which you want the strain and stress:');  
    N=[1-3*((x-(z-1)*l)/l)+2*(x^2 +(z-1)^2*l*l -2*x*(z-1)*l)*(1/l^2) 4*((x-(z-1)*l)/l)-4*(x-(z-1)*l)^2*(1/l^2) (-((x-(z-1)*l)/l)+2*(x-(z-1)*l)^2*(1/l^2))];
   B=diff(N);
    strain=B*displacements(2*z-1:2*z+1,1);
    stress= A*E*strain;
    disp(strain);
    disp(stress);
end
% Add interpolation
xc=input('Enter the value of length where you want to get the interpolated displacement:');
ElementNo=ceil(xc/l);
% modified x xm
xm=xc-(ElementNo-1)*l;
interpolated_displacement=(1-3*(xm/l)+2*(xm/l)^2)*(displacements(2*ElementNo-1)) +(4*(xm/l)-4*(xm/l)^2)*((displacements(2*ElementNo)))+ ((-(xm/l)+2*(xm/l)^2))*(displacements(2*ElementNo+1));
disp(interpolated_displacement);
end
