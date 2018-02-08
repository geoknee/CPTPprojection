function [ choi_ml_vec, outside_solution, inside_solution, outside_costs, inside_costs ] = gdap( A,n )
%gdap projected gradient descent algorithm
%   Detailed explanation goes here
    d = sqrt(sqrt(size(A)));
    d = d(2);
    choi_init = eye(d*d)/d;
    choi_init = reshape(choi_init,[],1);
    inside_solution  = {choi_init};
    outside_solution = {choi_init};
    stepsize      = 1.0/(1e1*d);
    for i=1:3e3
%         cost(A,n,inside_solution{i});
        inside_costs(i)     = cost(A,n,inside_solution{i});
        outside_costs(i)     = cost(A,n,outside_solution{i});
        outside_solution{i+1}  = outside_solution{i}-stepsize*gradient(A,n,outside_solution{i});
        inside_solution{i+1}   = CPTP_project(outside_solution{i+1});
%         solution{i+1}  = solution{i}-stepsize*gradient(A,n,solution{i});
        if i>20
            if var(inside_costs(i-20:i)) < 1e-5
                break
            end
    end
%     
%     end = time.time()
    choi_ml_vec = inside_solution{end};
end

