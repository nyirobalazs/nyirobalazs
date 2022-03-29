clear
%% Read csv
filename = 'budapest_connectome_3.0_209_0_median_coun_and_length.csv';
budapest_dataset = readtable(filename);
N = 83;
WEIGHT_MATRIX = zeros([N,N]);

% Create weight matrix
table_size = size(budapest_dataset); 
rows = table_size(1); 

for row = 1:rows
    WEIGHT_MATRIX(budapest_dataset{row,5},budapest_dataset{row,6}) = budapest_dataset{row,10};
    WEIGHT_MATRIX(budapest_dataset{row,6},budapest_dataset{row,5}) = budapest_dataset{row,10};
end


%% Calculate reatve change
%define variables
c_array         = zeros([1,N]);
c_array(26)     = 0.025;
c_array(68)     = 0.025;
q_array         = zeros([1,N]);
tspan           = 1:30;
beta            = [0.25, 0, 4];

% ODE function
input           = [c_array; q_array;WEIGHT_MATRIX];
RESULTS         = {};

for BETA = beta
    [t, out]        = ode45(@(t,C0) yearly_chang_ode(t,C0,BETA),tspan,input);
    RESULTS{end+1}  = out;
end


%% Create plots

% FIGURE 2
ind = 1;
linestyles = {'-'; ':'; '--'};

for OUT = RESULTS
    %Separate the outputs
    C_OUT                = OUT{1,1}(:,1:N);
    Q_OUT                = OUT{1,1}(:,N+1:2*N);
    W_OUT                = OUT{1,1}(:,(2*N)+1:2*N+N^2);
    
    %Calculate and plot average concentration
    C_T                  = (1/N)*sum(normalize(C_OUT,'range'),2);
    fig                  = plot(t,C_T)
    fig.LineStyle        = linestyles{ind}
    fig.Color            = 'b'
    fig.LineWidth        = 1
    hold on
    
    %Calculate and plot average damage
    Q_T                  = (1/N)*sum(normalize(Q_OUT,'range'),2);
    fig                  = plot(t,Q_T)
    fig.LineStyle        = linestyles{ind}
    fig.Color            = 'r'
    fig.LineWidth        = 1
    hold on

    %Calculate and plot scaled average connection weight
    w_t                  = (1/(N^2)).*sum(normalize(W_OUT,'range'),2);
    w_0                  = (1/(N^2)).*sum(sum(normalize(WEIGHT_MATRIX,'range'),2),1);
    W_T                  = normalize(w_t/w_0,'range');
    fig                  = plot(t,W_T)
    fig.LineStyle        = linestyles{ind}
    fig.Color            = 'g'
    fig.LineWidth        = 1
    hold on
    
    ind = ind + 1;
end

% few extra setting for fig. 2
title('Evolution of averaged toxic concentration and demage')
xlabel('Time T (yr)')
xlim([1 30])
ylim([0.005 1.05])
legend(append('C(\beta=',num2str(beta(1)),',\gamma=',num2str(beta(1)/2),')'),...
       append('Q(\beta=',num2str(beta(1)),',\gamma=',num2str(beta(1)/2),')'),...
       append('C(\beta=',num2str(beta(1)),',\gamma=',num2str(beta(1)/2),')'),...
       append('C(\beta=',num2str(beta(2)),',\gamma=',num2str(beta(2)/2),')'),...
       append('Q(\beta=',num2str(beta(2)),',\gamma=',num2str(beta(2)/2),')'),...
       append('C(\beta=',num2str(beta(2)),',\gamma=',num2str(beta(2)/2),')'),...
       append('C(\beta=',num2str(beta(3)),',\gamma=',num2str(beta(3)/2),')'),...
       append('Q(\beta=',num2str(beta(3)),',\gamma=',num2str(beta(3)/2),')'),...
       append('C(\beta=',num2str(beta(3)),',\gamma=',num2str(beta(3)/2),')'));

% FIGURE 3

%Defining brain areas
frontal         = [1, 4:11, 42, 45:54];
parietal        = [16:20, 22, 55:61];
occipital       = [23, 63:65];
temporal        = [25:27, 29:32, 34:35,40, 66:76, 81];
limbic          = [12:15, 82];
basa_ganglia    = [36:39, 77:80];
brain_stem      = 83;

% Calculate Q in every brain regions
Q_FRONT         = normalize(sum(RESULTS{1,1}(:,N + frontal),2),'range');
Q_PARI          = normalize(sum(RESULTS{1,1}(:,N + parietal),2),'range');
Q_OCC           = normalize(sum(RESULTS{1,1}(:,N + occipital),2),'range');
Q_TEMP          = normalize(sum(RESULTS{1,1}(:,N + temporal),2),'range');
Q_LIMB          = normalize(sum(RESULTS{1,1}(:,N + limbic),2),'range');
Q_BASG          = normalize(sum(RESULTS{1,1}(:,N + basa_ganglia),2),'range');
Q_BS            = normalize(sum(RESULTS{1,1}(:,N + brain_stem),2),'range');

% Make the plot
figure(2)
plot(tspan,Q_FRONT,tspan,Q_PARI,tspan,Q_OCC,tspan,Q_TEMP,tspan,Q_LIMB,tspan,Q_BASG,tspan,Q_BS)
title('Evolution of damage in different brain regions')
xlabel('Time T (yr)')
ylabel('Q_{s}')
xlim([1 30])
ylim([0.005 1.05])
legend('Frontal','Parietal','Occipital','Temporal','Limbic','Basal ganglia','Brain stem')


function output = yearly_chang_ode(t,input,beta)

    % Constants
    alpha = 0.75;
    gamma = beta/2;
    p = 0.01;

    N = 83;
    C = input(1:N);
    Q = input(N+1:2*N);
    W = (2*N+1:2*N+N^2);
    Q_PLUS = Q + Q';
    W=reshape(W,N,N);
    
    %Calculate Laplacian matrix
    diag_matrix = diag(sum(W,2));
    L = p*(diag_matrix-W);

    %Calculations  

      CDOT     = -L*C+alpha*C.*(1-C); 
      QDOT     = beta*C.*(1-Q);
      WDOT     = -gamma*W.*Q_PLUS;

    output     = [CDOT;QDOT;reshape(WDOT,N^2,1)];
    
end