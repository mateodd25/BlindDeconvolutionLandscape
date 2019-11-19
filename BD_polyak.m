function [ w, x, func_vals ] = BD_polyak( L, R, y, w0, x0 )
%BD_polyak Runs Polyak's subgradient method for the blind deconvolution
% problem initialized at (w0, x0)

max_iter = 100000;
w = w0; x = x0;
func_vals = zeros(max_iter,1);
m = size(L,1);
for kk = 1:max_iter 
    if kk == 125
        fprintf('');
    end
    diff = (L*w).*(R*x) - y;
    signs = sign(diff);
    val_fun = norm(diff,1)/m;
    if val_fun < 1e-10
        break;
    end
    if val_fun > 1e30 
        break
    end
    gw = L'*((R*x).*signs);
    gx = R'*((L*w).*signs);
    
    w = w - val_fun*gw/norm(gw)^2;
    x = x - val_fun*gx/norm(gx)^2;
%     if mod(kk,1) == 0
%         fprintf('Iter %d objective val_fun %7.2e \n', kk, val_fun);
%     end
end
func_vals = func_vals(1:kk);

end

