function x = newton_raphson (fun, jacobian, x, eps)
% Remember to use x0 in the function call!!!!

fun_value = fun(x);
fun_norm = norm(fun_value);

while fun_norm > eps
    delta = jacobian(x)\(-fun_value);
    x = x + delta;
    fun_value = fun(x);
    fun_norm = norm(fun_value);
end

end