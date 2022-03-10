function [tau,i] = custom_solver(equation, var, i_initial, step, max_time)
    t = var;
    tau = [];
    i = 0;
    while isempty(tau)
        i = i + 1;
        if step*i > max_time
            tau = inf;
            return
        end
        tau = vpasolve(equation, t, [i_initial*step, i_initial*step + step*i]);
    end
    tau = double(tau);
end

