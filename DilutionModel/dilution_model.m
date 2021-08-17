N_T = 1e8;
N = N_T*[ 0.3, 0.7 ]; % pop size {cells}
s = [ 0.72 0.7 ]; % fitness (absolute) {per gen}
nu = [ 0.1 0.7 ]; % prob filtration

QbyV = (N*s')/(N*nu'); % filtration amount {cell frac per gen} 

Ngen = 5;

QbyV0 = QbyV;
Nss = 1000;
dt = 1/Nss; % time in {gen}

N0 = N;

K = s*[ 1; -1 ] - QbyV*(nu*[ 1; -1 ]);



Nframes = 60;
tot_it = Nss*Ngen;
it_per_frame = floor(tot_it/Nframes);

Nt = zeros(2, Nframes);
QbyVt = zeros(1, Nframes);
gent = zeros(1, Nframes);

for k = 1:Ngen
    
    
    for u = 1:Nss % defined as derivatives so care must be taken to integrate small steps
        
        it = (k - 1)*Nss + u;
        
        dN1 = (N(1)*(1 - N(1)/N_T)*K)*dt;
        
        if max(N(2) - dN1, 0) == 0
            dN1 = max(N(2), 0);
        elseif max(N(1) + dN1, 0) == 0
            dN1 = -max(N(1), 0);
        end
        
        N(2) = N(2) - dN1;         
        N(1) = N(1) + dN1;        
        
        QbyV = (N*s')/(N*nu');
        K = s*[ 1; -1 ] - QbyV*(nu*[ 1; -1 ]);
        
        if mod(it,it_per_frame) == 0
            it_id = it/it_per_frame;
            Nt(:, it_id) = N';
            QbyVt(it_id) = QbyV;
            gent(it_id) = (k - 1) + u*dt;
        end
    end
%     
%     QbyVt(k) = QbyV;    
%     Nt(:,k) = N';
end

figure(1)
plot([0, gent], [ N0', Nt ])
xlabel('Time (gens)')
ylabel('Population count')
legend('Mutant pop', 'Last fixed pop')
figure(2)
plot([0, gent], [ QbyV0, QbyVt ])
xlabel('Time (gens)','interpreter', 'latex', 'fontsize', 14)
ylabel('$\frac{Q}{V}$, filtration (cell frac per gen)', 'interpreter', 'latex', 'fontsize', 14)

% test analytical solution
% h = Nt(1,:)./Nt(2,:);
% plot(1:Ngen, log(h))
% hold on
% plot(1:Ngen, log(N0(1)/N0(2)*exp(QbyVt.*(1:Ngen)))) % actually need to
% integrate
% hold off
