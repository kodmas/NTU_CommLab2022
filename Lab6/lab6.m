%% 1 symbol mapping
sym_seq = symbol_mapper('0',2,1,'PSK');
disp(sym_seq);

%% 2a
ll = 10000;
simbitstr = '';
for i =1:ll
    seed = rand;
    if seed > 0.5
        simbitstr = strcat(simbitstr,'1');
    else
        simbitstr = strcat(simbitstr,'0');
    end
end

sim_sym_seq = symbol_mapper(simbitstr,4,1,'PSK');
a2_x = zeros(1,5000);
a2_y = zeros(1,5000);
for i = 1:5000
    a2_x(i) = sim_sym_seq{i}(1);
    a2_y(i) = sim_sym_seq{i}(2);
end
figure;
h = histogram2(a2_x,a2_y);

%%receiver side
%Eb/N0 = 0dB
for i = 1:5000
    nx = normrnd(0,1/4);
    ny = normrnd(0,1/4);

    a2_x(i) = sim_sym_seq{i}(1) + nx;
    a2_y(i) = sim_sym_seq{i}(2) + ny;
end

figure;
h1 = histogram2(a2_x,a2_y);

%Eb/N0 = 10dB
for i = 1:5000
    nx = normrnd(0,1/40);
    ny = normrnd(0,1/40);

    a2_x(i) = sim_sym_seq{i}(1) + nx;
    a2_y(i) = sim_sym_seq{i}(2) + ny;
end

figure;
h2 = histogram2(a2_x,a2_y);
for i = 1:12
    morebins(h2);
end



%Eb/N0 = 20dB
for i = 1:5000
    nx = normrnd(0,1/400);
    ny = normrnd(0,1/400);
  
    a2_x(i) = sim_sym_seq{i}(1) + nx;
    a2_y(i) = sim_sym_seq{i}(2) + ny;
end

figure;
h3 = histogram2(a2_x,a2_y);
for i = 1:25
    morebins(h3);
end

%% 2b
try_l = 10000;
try_bin = generate_bin(try_l);
try_sym = symbol_mapper(try_bin,4,1,'PSK');
try_o = receiver(try_sym,1/4,4);
try_dec = MD_symbol_demapper(try_o,4,1,'PSK');
disp(bin_seq(1:10));
SER = calculate_SER(try_bin,try_dec,4);
disp(SER);

%% 3
Eb_N0 = [1:10];
in_dB = [power(10,Eb_N0/10)];
%disp(in_dB);
test_l = 10000;
test_bin = generate_bin(test_l);

% PAM
M_list_1 = [2,4,8,16];
M_list_2 = [4,16,64];
ii1 = 1;
for temp_M = M_list_1
    theo_PAM = 2*qfunc(6*log2(temp_M)/(power(temp_M,2)-1)*Eb_N0);
    temp_Eb = (power(temp_M,2)-1)/(12*log2(temp_M));
    test_sym = symbol_mapper(test_bin,temp_M,1,'PAM');
    for temp_dB = in_dB
        temp_N0 = temp_Eb/temp_dB;
        temp_o  = receiver(test_sym,temp_N0,temp_M);
        temp_dec = MD_symbol_demapper(temp_o,temp_M,1,'PAM');
        temp_SER = calculate_SER(test_bin,temp_dec,temp_M);
        SER_list_PAM(ii1) = temp_SER;
        ii1 = ii1 + 1;
    end
end
PAM_2 = SER_list_PAM(1:10);
PAM_4 = SER_list_PAM(11:20);
PAM_8 = SER_list_PAM(21:30);
PAM_16 = SER_list_PAM(31:40);

figure();
theor_PAM_2 = 2*qfunc(sqrt(6*log2(2)/(power(2,2)-1)*in_dB));
semilogy(Eb_N0,theor_PAM_2,'ro-.');
hold on
semilogy(Eb_N0,PAM_2,'-o','MarkerFaceColor','red');
hold on 
theor_PAM_4 = 2*qfunc(sqrt(6*log2(4)/(power(4,2)-1)*in_dB));
semilogy(Eb_N0,theor_PAM_4,'bo-.');
hold on
semilogy(Eb_N0,PAM_4,'-o','MarkerFaceColor','blue');
hold on 
theor_PAM_8 = 2*qfunc(sqrt(6*log2(8)/(power(8,2)-1)*in_dB));
semilogy(Eb_N0,theor_PAM_8,'go-.');
hold on
semilogy(Eb_N0,PAM_8,'-o','MarkerFaceColor','green');
hold on 
theor_PAM_16 = 2*qfunc(sqrt(6*log2(16)/(power(16,2)-1)*in_dB));
semilogy(Eb_N0,theor_PAM_16,'yo-.');
hold on
semilogy(Eb_N0,PAM_16,'-o','MarkerFaceColor','yellow');
hold on 
xlabel("Eb/N0 : SNR value")
ylabel("SER(N0): the symbol error rate")
legend('theoretical curve for M = 2','simulation result for M = 2','theoretical curve for M = 4', ...
    'simulation result for M = 4','theoretical curve for M = 8','simulation result for M = 8' ...
    ,'theoretical curve for M = 16','simulation result for M = 16','Location','best');
% PSK
ii2 = 1;
for temp_M = M_list_1
    %theo_PAM = 2*qfunc(6*log2(temp_M)/(power(temp_M,2)-1)*Eb_N0);
    temp_Eb = 1/(power(sin(pi/temp_M),2)*log2(temp_M)*4);
    test_sym = symbol_mapper(test_bin,temp_M,1,'PSK');
    for temp_dB = in_dB
        temp_N0 = temp_Eb/temp_dB;
        temp_o  = receiver(test_sym,temp_N0,temp_M);
        temp_dec = MD_symbol_demapper(temp_o,temp_M,1,'PSK');
        temp_SER = calculate_SER(test_bin,temp_dec,temp_M);
        SER_list_PSK(ii2) = temp_SER;
        ii2 = ii2 + 1;
    end
end
PSK_2 = SER_list_PSK(1:10);
PSK_4 = SER_list_PSK(11:20);
PSK_8 = SER_list_PSK(21:30);
PSK_16 = SER_list_PSK(31:40);

figure();
theor_PSK_2 = 2*qfunc(sqrt(2*log2(2)*(power(sin(pi/2),2))*in_dB));
semilogy(Eb_N0,theor_PSK_2,'ro-.');
hold on
semilogy(Eb_N0,PSK_2,'-o','MarkerFaceColor','red');
hold on 
theor_PSK_4 = 2*qfunc(sqrt(2*log2(4)*(power(sin(pi/4),2))*in_dB));
semilogy(Eb_N0,theor_PSK_4,'bo-.');
hold on
semilogy(Eb_N0,PSK_4,'-o','MarkerFaceColor','blue');
hold on 
theor_PSK_8 = 2*qfunc(sqrt(2*log2(8)*(power(sin(pi/8),2))*in_dB));
semilogy(Eb_N0,theor_PSK_8,'go-.');
hold on
semilogy(Eb_N0,PSK_8,'-o','MarkerFaceColor','green');
hold on 
theor_PSK_16 = 2*qfunc(sqrt(2*log2(16)*(power(sin(pi/16),2))*in_dB));
semilogy(Eb_N0,theor_PSK_16,'yo-.');
hold on
semilogy(Eb_N0,PSK_16,'-o','MarkerFaceColor','yellow');
hold on 
xlabel("Eb/N0 : SNR value")
ylabel("SER(N0): the symbol error rate")
legend('theoretical curve for M = 2','simulation result for M = 2','theoretical curve for M = 4', ...
    'simulation result for M = 4','theoretical curve for M = 8','simulation result for M = 8' ...
    ,'theoretical curve for M = 16','simulation result for M = 16','Location','best');

% QAM
ii3 = 1;
for temp_M = M_list_2
    %theo_PAM = 2*qfunc(6*log2(temp_M)/(power(temp_M,2)-1)*Eb_N0);
    temp_Eb = (temp_M-1)/(6*log2(temp_M));
    test_sym = symbol_mapper(test_bin,temp_M,1,'QAM');
    for temp_dB = in_dB
        temp_N0 = temp_Eb/temp_dB;
        temp_o  = receiver(test_sym,temp_N0,temp_M);
        temp_dec = MD_symbol_demapper(temp_o,temp_M,1,'QAM');
        temp_SER = calculate_SER(test_bin,temp_dec,temp_M);
        SER_list_QAM(ii3) = temp_SER;
        ii3 = ii3 + 1;
    end
end
QAM_4 = SER_list_QAM(1:10);
QAM_16 = SER_list_QAM(11:20);
QAM_64 = SER_list_QAM(21:30);

figure();
theor_QAM_4 = 4*qfunc(sqrt(3*log2(4)/(4-1))*in_dB);
semilogy(Eb_N0,theor_QAM_4,'ro-.');
hold on
semilogy(Eb_N0,QAM_4,'-o','MarkerFaceColor','red');
hold on 
theor_QAM_16 = 4*qfunc(sqrt(3*log2(16)/(16-1))*in_dB);
semilogy(Eb_N0,theor_QAM_16,'go-.');
hold on
semilogy(Eb_N0,QAM_16,'-o','MarkerFaceColor','green');
hold on 
theor_QAM_64 = 4*qfunc(sqrt(3*log2(64)/(64-1))*in_dB);
semilogy(Eb_N0,theor_QAM_64,'yo-.');
hold on
semilogy(Eb_N0,QAM_64,'-o','MarkerFaceColor','yellow');
hold on 
xlabel("Eb/N0 : SNR value")
ylabel("SER(N0): the symbol error rate")
legend('theoretical curve for M = 4', 'simulation result for M = 4' ...
    ,'theoretical curve for M = 16','simulation result for M = 16', ...
    'theoretical curve for M = 64','simulation result for M = 64','Location','best');

%% 4
%a
impulse_response_c = [1,0,1;1,1,1];
BPSK_l = 10000;
BPSK_bin= zeros(1,BPSK_l);
%generate input bitstring for 4a
for i =1:BPSK_l
    seed = rand;
    if seed > 0.5
        BPSK_bin(1,i) = 1;
    else
        BPSK_bin(1,i) = 0;
    end
end

encoded_data = conv_enc(BPSK_bin, impulse_response_c);
disp(encoded_data(1:10));
for i = 1:size(encoded_data,2)
    if encoded_data(i) == 0
        BPSK_soft{i} = [0,1/2]; 
    else
        BPSK_soft{i} = [0,-1/2]; 
    end
end

new_Eb_N0 = linspace(1,20,15);

in_dB = [power(10,new_Eb_N0/10)];
ii4 = 1;
for d = in_dB
    mm_N0 =0.25/d;

    BPSK_o = receiver(BPSK_soft,mm_N0,2);
    
    BPSK_temp_dec = MD_symbol_demapper(BPSK_o,2,1,'PSK');
    for i = 1:2*BPSK_l
        BPSK_dec(1,i) = str2double(BPSK_temp_dec(1,i));
    end
    
    BPSK_fin = conv_dec(BPSK_dec,impulse_response_c);
    
    tempBER = 0;
    for k = 1:BPSK_l
        if BPSK_bin(k)~=BPSK_fin(k)
            tempBER = tempBER + 1;
        end
    end
    tempBER = tempBER / BPSK_l;
    temp_BER_list(ii4) = tempBER;
    ii4 = ii4+1;
end
semilogy(new_Eb_N0,temp_BER_list,'go-.');
xlabel("Eb/N0 : SNR value")
ylabel("BER(N0): the bit error rate")

%b
%BPSK_o = receiver(BPSK_soft,1/400,2);
%disp('-----');
%disp(BPSK_o(1:10));

%% -------------------------------------------- %%

%% generate bin
function test_bin = generate_bin(ll)
    %ll = 10000;
    test_bin = '';
    for i =1:ll
        seed = rand;
        if seed > 0.5
            test_bin = strcat(test_bin,'1');
        else
            test_bin = strcat(test_bin,'0');
        end
    end
end

%% receive noise
function output = receiver(sim_sym_seq,var,M)
    le = size(sim_sym_seq,2)/log2(M);
    for i = 1:le
        nx = normrnd(0,var);
        ny = normrnd(0,var);
      
        noised_x(i) = sim_sym_seq{i}(1) + nx;
        noised_y(i) = sim_sym_seq{i}(2) + ny;
    end
    for i = 1:size(noised_x,2)
        output{1,i} = [noised_x(i),noised_y(i)];    
    end
end

%% calculate SER
function SER = calculate_SER(simbitstr,bin_seq,M)
    l = log2(M);
    SER = 0;
    for i = 1:size(bin_seq,2)
        temp_str = simbitstr((i-1)*l+1:i*l);
        if(~strcmp(temp_str,bin_seq(i)))
            SER = SER + 1;
        end
    end
    SER = SER/size(bin_seq,2);
end

%% gray code function
function g = graycoding(b)
    g(1) = b(1);
    for i = 2 : length(b)
        x = xor(str2num(b(i-1)), str2num(b(i)));
        g(i) = num2str(x);
    end
end

%% symbol_mapper
function sym_seq = symbol_mapper(bin_seq,M,d,name)
    if(strcmp(name,'PAM'))
        sqrt_Ep = d/2;
        l = log2(M);
        for i = 1:strlength(bin_seq)/l
            temp = bin_seq((i-1)*l+1:i*l);
            %disp(temp)
            gray_encoded_m = bin2dec(graycoding(temp));
            Am = 2*(gray_encoded_m+1)-1-M;
            sym_seq{i} = [Am*sqrt_Ep,0];
        end
    end
    if(strcmp(name,'PSK'))
        sqrt_Ep = d/(2*sin(pi/M));
        l = log2(M);
        for i = 1:strlength(bin_seq)/l
            temp = bin_seq((i-1)*l+1:i*l);
            %disp(temp)
            gray_encoded_m = bin2dec(graycoding(temp));
            %disp(gray_encoded_m);
            Amx = cos((2*pi*(gray_encoded_m)+pi)/M);
            Amy = sin((2*pi*(gray_encoded_m)+pi)/M);
            
            sym_seq{i} = [Amx*sqrt_Ep,Amy*sqrt_Ep];
        end
    end
    if(strcmp(name,'QAM'))
        sqrt_Ep = d/2;
        l = log2(M);
         for i = 1:strlength(bin_seq)/l
            temp = bin_seq((i-1)*l+1:i*l);
            %disp(temp)
            gray_encoded_m = bin2dec(graycoding(temp));
            %disp(gray_encoded_m);
            %[r,q] = quorem(sym(gray_encoded_m),sym(sqrt(M))); 
            r = floor(gray_encoded_m/sqrt(M));
            q = mod(gray_encoded_m,sqrt(M));
            Ami = 2*r-1-M;
            Amq = 2*q-1-M;
            
            sym_seq{i} = [Ami*sqrt_Ep,Amq*sqrt_Ep];
        end
    end
end

function bin_seq = MD_symbol_demapper(sym_seq, M,d,name)
    l = log2(M);
    bin_seq = "";
    %generate gray_code list / constellation base points
    for i = 1:M
        tmp_str = num2str(graycoding(dec2bin(i-1)));
        while(strlength(tmp_str)<l)
            tmp_str = strcat('0',tmp_str);
        end
        %disp(tmp_str);
        gray_list{i} = tmp_str;
        const_base{i} = symbol_mapper(tmp_str,M,d,name);
    end
    %%calculate minimum distance
    for i = 1:size(sym_seq,2)
        for j = 1:size(const_base,2)
            dist(j) = power((sym_seq{i}(1) - const_base{j}{1}(1)),2) + power((sym_seq{i}(2) - const_base{j}{1}(2)),2);
        end
        mini = min(dist);
        ind = find(dist == mini);
        bin_seq(i) = gray_list{1,ind};
        %bin_seq = strcat(bin_seq,gray_list{1,ind});
    end
end

%% ------------function from lab5--------------- %%
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