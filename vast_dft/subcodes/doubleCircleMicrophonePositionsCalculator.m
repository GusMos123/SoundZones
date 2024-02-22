function mpos = doubleCircleMicrophonePositionsCalculator(num_microphones,inner_circle_radius,outer_circle_radius,center_point,rotation_angle)
% function mpos = dcirc(M,radi,rado,cent,rot)
% Double circle microphone positions
% M        : the number of microphone positions, default is 6
% radi     : the radius of the inner circle, default is 0.2
% rado     : the radius of the outer circle, default is 0.5
% cent     : the center of two circles, concentric, default is [0 0]
% rot      : rotation angle in radian, default is 0
if nargin == 0
    num_microphones = 6;
    inner_circle_radius = 0.2;
    outer_circle_radius = 0.5;
    center_point = [0 0];
    rotation_angle = 0;
elseif nargin == 1
    inner_circle_radius = 0.2;
    outer_circle_radius = 0.5;
    center_point = [0 0];
    rotation_angle = 0;
elseif nargin == 2
    outer_circle_radius = inner_circle_radius+0.2;
    center_point = [0 0];
    rotation_angle = 0;
elseif nargin == 3
    center_point = [0 0];
    rotation_angle = 0;
elseif nargin == 4
    rotation_angle = 0;
elseif nargin == 5
else
    error('See help')
end

% radi = 0.2;
mindx = 0:num_microphones/2-1;
rtheta = inner_circle_radius.*exp(1j*mindx'*2*pi/(num_microphones/2)+1j*rotation_angle);

% rado = 0.5;
rtheta2 = outer_circle_radius.*exp(1j*(mindx*2+1)'*pi/(num_microphones/2)+1j*rotation_angle);

mpos = zeros(num_microphones,2);
mpos(:,1) = center_point(1)+[real(rtheta);real(rtheta2)];
mpos(:,2) = center_point(2)+[imag(rtheta);imag(rtheta2)];
end

%{
circ1 = radius.*exp(1j*(0:5:360)'*2*pi/360);
circ2 = radius2.*exp(1j*(0:5:360)'*2*pi/360);


figure
plot(real(circ1),imag(circ1))
hold on
plot(real(circ2),imag(circ2))
scatter(mpos1x,mpos1y)
scatter(mpos2x,mpos2y);axis equal
xlim([-1 1]);ylim([-1 1])
grid minor
%}
