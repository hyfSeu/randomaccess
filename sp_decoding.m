function [orignal_code,hamming_distance,corres_idx] = sp_decoding(Hhat,orignal_message,delta)
    [K,L] = size(Hhat);
    temp_code = zeros(K,L);
    hamming_distance = [];
    corres_idx = [];
    if delta == 0.25
        message = zeros(K,L/2);
        P = L/4;
        for i = 1:K
            for j = 1:P
                biggest_idx = find(abs(Hhat(i,(4*j-3):4*j)) == max(abs(Hhat(i,(4*j-3):4*j))));
                if(length(biggest_idx)>1)
                    biggest_idx = biggest_idx(1);
                end
                temp_code(i,4*(j-1)+biggest_idx) = 1;
                switch 2^(4-biggest_idx)
                    case 8
                        message(i,2*j-1:2*j) = [1,1];
                    case 4
                        message(i,2*j-1:2*j) = [1,0];
                    case 2
                        message(i,2*j-1:2*j) = [0,1];
                    case 1
                        message(i,2*j-1:2*j) = [0,0];
                end
            end
            current_distance = pdist2(message(i,:),unique(orignal_message,'rows'),'hamming');
            temp_idx = find(current_distance==min(current_distance));
            corres_idx = [corres_idx;temp_idx(1)];
            hamming_distance = [hamming_distance;min(current_distance)];
        end
    elseif delta == 0.125
        message = zeros(K,3*L/8);
        P = L/8;
        for i = 1:K
            for j = 1:P
                biggest_idx = find(abs(Hhat(i,(8*j-7):8*j)) == max(abs(Hhat(i,(8*j-7):8*j))));
                if(length(biggest_idx)>1)
                    biggest_idx = biggest_idx(1);
                end
                temp_code(i,8*(j-1)+biggest_idx) = 1;
                switch 2^(8-biggest_idx)
                    case 128
                        message(i,3*j-2:3*j) = [1,1,1];
                    case 64
                        message(i,3*j-2:3*j) = [1,1,0];
                    case 32
                        message(i,3*j-2:3*j) = [1,0,1];
                    case 16
                        message(i,3*j-2:3*j) = [1,0,0];
                    case 8
                        message(i,3*j-2:3*j) = [0,1,1];
                    case 4
                        message(i,3*j-2:3*j) = [0,1,0];
                    case 2
                        message(i,3*j-2:3*j) = [0,0,1];
                    case 1
                        message(i,3*j-2:3*j) = [0,0,0];
                end
            end
            current_distance = pdist2(message(i,:),unique(orignal_message,'rows'),'hamming');
            temp_idx = find(current_distance==min(current_distance));
            corres_idx = [corres_idx;temp_idx(1)];
            hamming_distance = [hamming_distance;min(current_distance)];
        end
    end
    correct_idx = hamming_distance<=0.25;
    message = message(correct_idx,:);
    hamming_distance = hamming_distance(correct_idx,:);
    orignal_code = message;
end