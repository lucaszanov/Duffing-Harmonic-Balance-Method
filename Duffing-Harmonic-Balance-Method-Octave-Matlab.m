% Response computing via Harmonic Balance Method

% Frequency array
wi = 0 % Initial frequency [rad/s]
wf = 10 % Final frequency [rad/s]
dw = 0.01 % Frequency increment [rad/s]
w_hb = wi:dw:wf % Frequency array [rad/s]
w_hz_hb = w_hb/(2*pi) % Frequency array [rad/s]

% Mass
m = 1; % [kg]

% Linear stiffness
k = 25; % [N/m]

% Nonlinear stiffness
knl = 25; % [N/m^3]

% Damping
c = 0.2; % [N.s/m]

% Excitation amplitude
Y_low = 0.1; % [N]
Y_high = sqrt(2); % [N]

Y = Y_high;

a = 1;
b = 1;
XX = zeros(6,length(w_hb));
X = zeros(6,length(w_hb));

for j = 1:length(w_hb)
    % Polynomial solution (roots computing)
    coefs = [(9/16)*knl.^2,0,(3*k*knl-3*knl*m*w_hb(j).^2)/2,0,(m*m*(w_hb(j).^4)-2*k*m*(w_hb(j).^2)+c*c*(w_hb(j).^2)+k*k),0,-Y*Y];
    X(:,j) = roots(coefs);
    for i = 1:size(X,1)
        if isreal(X(i,j)) == 1
            XX(a,b) = X(i,j);
            a = a+1;
        endif
    endfor   
    if a>1
        b = b+1;
    endif
    a = 1;
endfor

XX = abs(XX);

clear X
X = zeros(6,length(w_hb));

for z =1:length(w_hb)
    aux = unique(XX(:,z));
    X(1:length(aux),z) = aux;
endfor

H = abs(X)/Y;

% Visualization

figure()
for i =1:6
  plot(transpose(w_hb)./(2*pi),(abs(H(i,:)))./max(abs(H(2,:))),'ob');
  hold on;
endfor  
ylim([0.01 1.2])
xlim([0.4 1.4])
grid on;
hold off;
xlabel('Frequency [Hz]')
ylabel('Amplitude [Normalized]')
%xlim([0.4 1.4])