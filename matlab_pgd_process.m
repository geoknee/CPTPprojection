for d=2:5   
    close all; clear A; clear n;
    % num_qubits = 1;
    % initialise ground truth
    choi_ground     = rand(d*d,d*d)-rand(d*d,d*d)+1.0j*rand(d*d,d*d)-1.0j*rand(d*d,d*d);
    % choi_ground     = 0.5*[1 0 0 1 ;0 1 1 0; 0 1 1 0; 1 0 0 1];
    % choi_ground = kron(U,U);
    % choi_ground= diag([0.4,0.4,0.2,1]); % this is CP but not TP (i.e. it is T I !!)
    % choi_ground= diag([0.4,0.4,0.4,0.4]); % this is CP but not TP (i.e. it is TD !!)
    choi_ground_vec = reshape(choi_ground,[],1);
    choi_ground_vec = CPTP_project(choi_ground_vec);
    choi_ground     = reshape(choi_ground_vec,[],d*d)

    addpath('./QETLAB-0.9')
    addpath('./QETLAB-0.9/helpers')
    % initialise preparations

    % preparations     = {};
    % num_preparations = 10*d^2;
    % for i =1:num_preparations
    %     ket       = rand(d,1)+1.0j*rand(d,1);
    %     ket       = ket/norm(ket);
    %     projector = ket*ket';
    %     preparations{i}=projector;
    % end

    % initialise measurements
    % measurements      = {};
    % num_measurements = 10*d^2;
    % i = 1;
    % for i=1:num_measurements
    %     ket       = rand(d,1)+1.0j*rand(d,1);
    %     ket       = ket/norm(ket);
    %     projector = ket*ket';
    %     measurements{i}=projector;
    % end

    % activate pauli measurements
    paulis{1} = [1 0;0 0];
    paulis{2} = [0 0;0 1];
    paulis{3} = 0.5*[1 1;1 1];
    paulis{4} = 0.5*[1 -1;-1 1];
    paulis{5} = [0.5 -0.5j ;0.5j 0.5];
    paulis{6} = [0.5 0.5j ;-0.5j 0.5];

    preparations = {};
    i = 1;
    for j=0:d-1
        for k=0:d-1
            if j==0 && k==0
                continue % don't need identity operator?? maybe do for trace decreasing
            end
            G = GenGellMann(j,k,d);
            [V,D] = eig(G);
            for l=1:d
                preparations{i} = V(:,l)*V(:,l)';
                i = i + 1;
            end
        end
    end
    % preparations = {1};
    % for i=1:num_qubits
    %     preparations = additional_paulis(preparations);
    % end


    measurements = preparations;
    % 
    % 
    % if num_qubits ==2
    %     l = 1;
    %     for i =1:6
    %         for j =1:6
    %             measurements{l} = kron(paulis{i},paulis{j});
    %             l = l +1;
    %         end
    %     end
    %     preparations = measurements;
    % end

    % build prepare and measure matrix
    % M = num_measurements * num_preparations;
    %A = zeros([M,d*d*d*d]);
    i = 1;
    num_measurements = length(measurements);
    num_preparations = length(preparations);
    for e=1:num_measurements
        for r=1:num_preparations
            E   = measurements{e};
            rho = preparations{r};
            row = reshape(kron(E',conj(rho)),[],1)'; % this was wrong before!
            A(i,:) = row;
            i = i+1;
        end
    end

    p = A*choi_ground_vec;
    p = p/sum(p);
    % n = mnrnd(1e4,p)';
    % n = n/sum(n);
    n = p;
    % [choi_ML, outside_solution, inside_solution, outside_costs, inside_costs] = gdap(A,n);
    tic;
    % [choi_ML, ~, solution, ~, costs] = gdap(A,n);
    [choi_ML, solution, costs] = gdapB(A,n);
    time_taken = toc;
    sprintf('condition number is')
    cond(A)
    choi_ML = reshape(choi_ML,[],d*d)
    error = trace_dist(choi_ML/trace(choi_ML),choi_ground/trace(choi_ground))
    scrsz = get(groot,'ScreenSize');
    figure('Position',[1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2])
    subplot(3,2,1);
    b=bar3(real(choi_ground));
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
    title('ground real')
    subplot(3,2,2);
    b=bar3(imag(choi_ground));
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
    title('ground imag')
    subplot(3,2,[3,4]);
    plot(costs);
    xlabel('iteration')
    ylabel('cost')
    subplot(3,2,5);
    b=bar3(real(choi_ML));
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
    title('solution real')
    subplot(3,2,6);
    b=bar3(imag(choi_ML));
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
    title('solution imag')
    %z = [inside_costs(:),outside_costs(:)].'; z=z(:); plot(linspace(1,length(inside_solution),length(z)),z)
    %%

    M = zeros([d*d,d*d*d*d]);
    for i=1:d
        e = zeros(1,d);
        e(i)  = 1;
        B = kron(eye(d),e);
        M = M + kron(B,B);
    end
    b = reshape(eye(d),[],1);

    tic;     
    % cvx_solver sdpt3
    cvx_solver mosek
    cvx_begin
    variable cvx_choi(d*d,d*d) hermitian semidefinite
    variable P(d,d) hermitian semidefinite
    % todo TP constraint
    choi_vec_cvx = reshape(cvx_choi,d*d*d*d,[]);
    n = real(n);
    p_cvx        = real(A*choi_vec_cvx);
    maximize(n'*log(p_cvx))
    subject to
    % eye(d)-reshape(M*choi_vec_cvx,d,d) ==  P; % this is for TNI
    M*choi_vec_cvx == b; % this is for TP
    cvx_end
    cvx_time = toc
    gdap_to_ground = trace_dist(choi_ML/trace(choi_ML),choi_ground/trace(choi_ground))
    gdap_to_cvx    = trace_dist(choi_ML/trace(choi_ML),cvx_choi/trace(cvx_choi))
    cvx_to_ground  = trace_dist(cvx_choi/trace(cvx_choi),choi_ground/trace(choi_ground))

    % tic;
    % [DIA_choi, ~, DIAcosts] = DIA(A,n);
    % dia_time = toc
    % DIA_choi = reshape(DIA_choi,[],d*d);
    % DIA_to_ground  = trace_dist(DIA_choi/trace(DIA_choi),choi_ground/trace(choi_ground))



    subplot(3,2,[3,4]);
    text(length(costs)*0.1,(max(costs)+min(costs))/2,[sprintf('d = %d process',d),char(10),sprintf('GDAP returned in %0.3f s',time_taken),char(10), sprintf('distance to ground %0.5f',error),char(10),sprintf('cvx returned in %0.5f s',cvx_time),char(10),sprintf('cvx to ground %0.5f',cvx_to_ground),char(10),sprintf('cvx to gdap %0.5f',gdap_to_cvx)])
    saveas(gcf, sprintf('d_%drandom_case.png',d));
end

