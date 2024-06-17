%2a
impulse_response = [1,0,0;1,0,1;1,1,1];
binary_data = [1,0,1,1,0];
encoded_data = conv_enc(binary_data, impulse_response);
disp('2a,encoded data');
disp(encoded_data);

%2b
d=[0,1,0,0,0,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,0];
decoded_data = conv_dec(d, impulse_response);
disp('2b,decoded data');
disp(decoded_data);

%2c
ll = 100000;
simbitstr = zeros(1,ll);
%generate input bitstring
for i =1:ll
    seed = rand;
    if seed > 0.5
        simbitstr(1,i) = 1;
    else
        simbitstr(1,i) = 0;
    end
end
%disp(simbitstr);

%encoding
simenc = conv_enc(simbitstr, impulse_response);


p = [];
p(1) = 0;
for i = 1:10
    p(1+i) = 0.1*i;
end

%simulate memoryless binary symmetric channel
BER_list = zeros(1,11);
noised_o = zeros(1,ll*size(impulse_response,1));
for i = 1:11
    for j = 1:ll*size(impulse_response,1)
        seed = rand;
        if seed < p(i)
            if simenc(j) == 0
                noised_o(j) = 1;
            else
                noised_o(j) = 0;
            end
        else
            noised_o(j) = simenc(j);
        end
    end
    %decoding
    simdec = conv_dec(noised_o,impulse_response);
    %calculate BER
    tempBER = 0;
    for k = 1:ll
        if simdec(k)~=simbitstr(k)
            tempBER = tempBER + 1;
        end
    end
    tempBER = tempBER / ll;
    BER_list(i) = tempBER;
end
disp('2c, BER list');
disp(BER_list);
figure;
plot(BER_list);

%3
impulse_response_c = [1,1,0;1,0,1];
%ll = 10;
simbitstr = zeros(1,ll);
for i =1:ll
    seed = rand;
    if seed > 0.5
        simbitstr(1,i) = 1;
    else
        simbitstr(1,i) = 0;
    end
end
%disp('simbitstr');
%disp(simbitstr);

simenc = conv_enc(simbitstr, impulse_response_c);

%disp('simenc');
%disp(simenc);

BER_list_c = zeros(1,11);
noised_o = zeros(1,ll*size(impulse_response_c,1));
for i = 1:11
    for j = 1:ll*size(impulse_response_c,1)
        seed = rand;
        if seed < p(i)
            if simenc(j) == 0
                noised_o(j) = 1;
            else
                noised_o(j) = 0;
            end
        else
            noised_o(j) = simenc(j);
        end
    end
    %disp('noised');
    %disp(noised_o);
    simdec = conv_dec(noised_o,impulse_response_c);
    tempBER = 0;
    for k = 1:ll
        if simdec(k)~=simbitstr(k)
            tempBER = tempBER + 1;
        end
    end
    tempBER = tempBER / ll;
    BER_list_c(i) = tempBER;
end
disp('3, BER list');
disp(BER_list_c);
figure;
plot(BER_list_c);

function encoded_data = conv_enc(binary_data, impulse_response)
    %k = 1;
    n = size(impulse_response,1);
    buflen = size(impulse_response,2);
    poss = power(2,buflen);
    %generate FSM_table
    for i = 1:poss
        bi = dec2bin((i-1),buflen);
        tem_buf = bi;

        FSM_table{i,1} = bi(1);
        FSM_table{i,2} = bi(2:buflen);

        for j = 1:n
            g = impulse_response(j,:);
            FSM_table{i,2+j} = 0;
            for k = 1:size(g,2)
               if(g(k) == 1)
                   FSM_table{i,2+j} = xor(FSM_table{i,2+j},str2double(tem_buf(k)));
               else
                   continue
               end
            end
        end
        FSM_table{i,6} = bi(1:2);
    end
    %disp('FSM');
    %disp(FSM_table);

    %traversing and encoding
    buf = zeros(1,buflen);
    encoded_data = [];
    for i = 1:size(binary_data,2)
        buf(1) = binary_data(i);
        for j = 1:poss
            if(strcmp(int2str(buf(1)),FSM_table{j,1}) && strcmp(strcat(int2str(buf(2)),int2str(buf(3))),FSM_table{j,2}))
                encoded_data = [encoded_data , FSM_table{j,3} , FSM_table{j,4} ,FSM_table{j,5}];
                buf(3) = buf(2);
                buf(2) = buf(1);
                break
            end
        end   
    end
end


function decoded_data = conv_dec(binary_data,impulse_response)
    buflen = size(impulse_response,2);
    
    %generate FSM_table
    n = size(impulse_response,1);
    poss = power(2,buflen);
    for i = 1:poss
        bi = dec2bin((i-1),buflen);
        tem_buf = bi;
        %disp(bi(1));
        FSM_table{i,1} = bi(1);
        FSM_table{i,2} = bi(2:buflen);
        for j = 1:n
            g = impulse_response(j,:);
            FSM_table{i,2+j} = 0;
            for k = 1:size(g,2)
               if(g(k) == 1)
                   FSM_table{i,2+j} = xor(FSM_table{i,2+j},str2double(tem_buf(k)));
               else
                   continue
               end
            end
        end
        FSM_table{i,6} = bi(1:2);
    end

    %generate state transition diagram
    for i = 1:4
        state_transition{i,1} = dec2bin((i-1),buflen-1);
        old = dec2bin((i-1),buflen-1);
        current = dec2bin((i-1),buflen-1);
        for j = 1:(buflen-2)
            current(buflen-j) = current(buflen-j-1);
        end
        current(1) = '0';
        state_transition{i,2} = bin2dec(current)+1;
        %state_transition{i,3} = bin2dec(current);
        for j = 1:poss
            if(strcmp(old,FSM_table{j,2}) && strcmp(current,FSM_table{j,6}))
                state_transition{i,4} = [FSM_table{j,3} , FSM_table{j,4} ,FSM_table{j,5}];
                break
            end
        end   
        current(1) = '1';
        state_transition{i,3} = bin2dec(current)+1;
        for j = 1:poss
            if(strcmp(old,FSM_table{j,2}) && strcmp(current,FSM_table{j,6}))
                state_transition{i,5} = [FSM_table{j,3} , FSM_table{j,4} ,FSM_table{j,5}];
                break
            end
        end   
    end

    %generate trellis_path table 'cost,from'
    for i  = 1:poss/2
        %trellis_table{i,1} = dec2bin((i-1),buflen);
        for j = 1:size(binary_data,2)/n
            trellis_table{i,j+1} = '';
        end
    end
    trellis_table{1,1} = [0,-1];
    

    for j = 1:size(binary_data,2)/n 
        for i = 1:poss/2
            if(strcmp(int2str(trellis_table{i,j}),''))
                continue
            else
                curr_y = binary_data(1,n*(j-1)+1:n*j);
                cost = trellis_table{i,j}(1);
                next_0 = state_transition{i,2};
                next_1 = state_transition{i,3}(1);
                code_0 = state_transition{i,4};
                code_1 = state_transition{i,5};
                error_0 = 0;
                error_1 = 0;
                for k = 1:size(curr_y,2)
                    if curr_y(k) ~= code_0(k)
                        error_0 = error_0 +1;
                    end
                    if curr_y(k) ~= code_1(k)
                        error_1 = error_1 +1;
                    end
                end
                % select 0 path
                cost_0 = cost + error_0;
              
               
                if(strcmp(int2str(trellis_table{next_0,j+1}),''))
                    trellis_table{next_0,j+1} = [cost_0,i];
                else
                    if(cost_0 < trellis_table{next_0,j+1}(1))
                        trellis_table{next_0,j+1} = [cost_0,i];
                    end
                end
                %select 1 path
                if j > size(binary_data,2)/n-2
                    continue
                else
                    cost_1 = cost + error_1;
                    if(strcmp(int2str(trellis_table{next_1,j+1}),''))
                        trellis_table{next_1,j+1} = [cost_1,i];
                    else
                        if(cost_1 < trellis_table{next_1,j+1}(1))
                            trellis_table{next_1,j+1} = [cost_1,i];
                        end
                    end
                end
            end
        end
    end
    %disp("trellis_table");
    %disp(trellis_table);
    
    %generate decoded data
    path = [1];
    curri = 1;
    pathlen = size(trellis_table,2);
    for i = 1:pathlen-1
        path(1+i) = trellis_table{curri,pathlen-i+1}(2);
        curri = trellis_table{curri,pathlen-i+1}(2);
    end
    path = flip(path);
    decoded_co = [];
    decoded_data = [];
    for i = 1:pathlen-1
        start = path(i);
        desti = path(i+1);
        if state_transition{start,2} == desti
            decoded_co = cat(2,decoded_co,state_transition{start,4});
            decoded_data(i) = 0;
        else
            decoded_co = cat(2,decoded_co,state_transition{start,5});
            decoded_data(i) = 1;
        end
    end
    %fin_cost = trellis_table{1,pathlen}(1);
end