function sp_message = sp_coding(message,delta)
%% this function is used to compute a sparse coding for 2:4 spreading rate
    [K,L] = size(message);
    if(delta == 0.25)
        if(mod(L,2))
            error("the length of message is wrong!");
        end
        block = L/2;
        sp_message = zeros(K,2*L);
        for i = 1:K
            for j =1:block
                switch 2*message(i,2*j-1)+message(i,2*j)
                    case 0
                        sp_message(i,4*j-3:4*j) = [0,0,0,1];
                    case 1
                        sp_message(i,4*j-3:4*j) = [0,0,1,0];
                    case 2
                        sp_message(i,4*j-3:4*j) = [0,1,0,0];
                    case 3
                        sp_message(i,4*j-3:4*j) = [1,0,0,0];
                end
            end
        end
    elseif(delta == 0.125)
        if(mod(L,3))
            error("the length of message is wrong!");
        end
        block = L/3;
        sp_message = zeros(K,8*L/3);
        for i = 1:K
            for j =1:block
                switch 4*message(i,3*j-2)+2*message(i,3*j-1)+message(i,3*j)
                    case 0
                        sp_message(i,8*j-7:8*j) = [0,0,0,0,0,0,0,1];
                    case 1
                        sp_message(i,8*j-7:8*j) = [0,0,0,0,0,0,1,0];
                    case 2
                        sp_message(i,8*j-7:8*j) = [0,0,0,0,0,1,0,0];
                    case 3
                        sp_message(i,8*j-7:8*j) = [0,0,0,0,1,0,0,0];
                    case 4
                        sp_message(i,8*j-7:8*j) = [0,0,0,1,0,0,0,0];
                    case 5
                        sp_message(i,8*j-7:8*j) = [0,0,1,0,0,0,0,0];
                    case 6
                        sp_message(i,8*j-7:8*j) = [0,1,0,0,0,0,0,0];
                    case 7
                        sp_message(i,8*j-7:8*j) = [1,0,0,0,0,0,0,0];
                end
            end
        end
    else
        error("the rates of code is wrong!");
    end
end