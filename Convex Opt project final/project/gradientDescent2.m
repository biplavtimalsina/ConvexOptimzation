function [ret_x,ret_itr,ret_diff_f]=gradientDescent2(param_Q,param_q)

x_act=  -param_Q\param_q;
f_act=(1/2) *x_act' *param_Q *x_act + param_q' *x_act;

inta=0.000001
epsilon=0.000001
size_Q=size(param_Q);
size_q=size(param_q);
eigen_Q=eig(param_Q);

if ~(all(eigen_Q)>0)
    disp('not a positive definite matrix')
end

eigen_max=max(eigen_Q);
eigen_min=min(eigen_Q);
cond_num=eigen_max/eigen_min;

if(size_Q(1)~=size_Q(2) || size_Q(1) ~=size_q(1))
    disp('Parameter Q and q have incorrect dimension')
end
x=rand(size_q(1),1);
grad_f=param_Q*x+param_q;
itr =0;

while (norm(grad_f))>inta
    grad_f=param_Q*x + param_q;
    search_dir=-grad_f;
    
    t_low=0;
    t_up=1;
    t=0;
    
    track_f=param_Q*(x+t_up .* search_dir) + param_q;
    
    while((search_dir' * track_f)<0)
        t_low=t_up;
        t_up=t_up*2;
    end
    
    while((t_up-t_up)>epsilon)
        t_mid=(1/2)*(t_low+t_up)
        track_f=param_Q * (x+t_up.*search_dir)+param_q;
        if((search_dir' *track_f)<0)
            t_low=t_mid;
        else
            t_up=t_mid;
        end
        t=(1/2)*(t_low+t_up);
    end
    
    x=x+t*search_dir;
    itr=itr+1
    
    f_itr=(1/2) *x' * param_Q *x+param_q'*x
    diff_f(itr)=f_itr -f_act;
end

ret_x=x;
ret_itr=itr;
ret_diff_f=diff_f;
