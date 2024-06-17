%2a
symbols = {   'a',   'b',   'c',   'd',   'e' ,'f','g','h'  };
prob    = [   0.2,   0.05,   0.005,   0.2,   0.3 ,0.05,0.045,0.15  ];
dict = huffman_dict( symbols, prob );
disp(dict);

%2b
sym_seq = 'gacab';
bin_seq = huffman_enc(sym_seq, dict);
disp(bin_seq);

%2c
sym_seq = huffman_dec(bin_seq, dict);
disp(sym_seq);

%3a
disp('Q3');
alphaseq = '';
for i=1:10
    seed = rand;
    x = sum(seed >= cumsum([0, prob]));
    alphaseq = strcat(alphaseq,symbols{x});
end
disp(alphaseq);
new_bin = huffman_enc(alphaseq,dict);
disp(new_bin);
disp(strlength(new_bin));

%3b
L = zeros(1,200);
avg = 0;
for j=1:200    
    alphaseq = '';
    for i=1:10
        seed = rand;
        x = sum(seed >= cumsum([0, prob]));
        alphaseq = strcat(alphaseq,symbols{x});
    end
    new_bin = huffman_enc(alphaseq,dict);
    L(j) = strlength(new_bin);
    avg = avg + L(j);
end
avg = avg/200;
%disp(avg);
%disp(L);
figure();
histogram(L);
title("average codeword length = "+ num2str(avg));

%3c
entropy = 2.53;
cal_avg = 2.6;
R = [10, 20, 50, 100, 200, 500, 1000];
n = [10,50,100];

lr_result = zeros(1,21);

for m = 1:size(n,2)
    curr_n = n(m);
    for k = 1:size(R,2)
        curr_R = R(k);
        L = zeros(1,curr_R);
        avg = 0;
        for j=1:curr_R    
            alphaseq = '';
            for i=1:curr_n
                seed = rand;
                x = sum(seed >= cumsum([0, prob]));
                alphaseq = strcat(alphaseq,symbols{x});
            end
            new_bin = huffman_enc(alphaseq,dict);
            L(j) = strlength(new_bin);
            avg = avg + L(j);
        end
        %disp(L);
        avg = avg/curr_R;
        avg = avg/curr_n;
        lr_result(7*(m-1)+k) = avg;
    end
end
%disp(lr_result);
lr_10n = lr_result(1:7);
lr_50n = lr_result(8:14);
lr_100n = lr_result(15:21);
%disp(lr_10n);
%disp(lr_50n);
%disp(lr_100n);
figure();
semilogx(R,lr_10n,'-o','MarkerFaceColor','red');
hold on 
semilogx(R,lr_50n,'-o','MarkerFaceColor','blue');
hold on 
semilogx(R,lr_100n,'-o','MarkerFaceColor','green');
hold on 
yline(entropy);
hold on 
yline(cal_avg);
xlabel("R : number of the Monte–Carlo runs")
ylabel("avgLn(R) : experimental average codeword length")
legend('n = 10','n = 50','n = 100','entropy','calculated-average codeword length');
function dict = huffman_dict(symbols,prob)
    for i = 1:size(symbols,2)
        table{i,1} = symbols{i};
        table{i,2} = prob(i);
    end
    %table = { 'a', 0.2; 'b',  0.05; 'c', 0.005;'d',0.2 ;'e', 0.3;'f' ,0.05;'g',0.045;'h',0.15  };
    nos  = size(table,1);
    weight_list = zeros(nos,1);
    for i = 1:1:nos
        weight_list(i) = table{i,2};
    end
    
    dict = cell(2*nos-1,5);
    %Huffman_tree{2*nos-1,6} = 0;
    
    for i = 1:nos
        dict{i,1} = table{i,1};
        dict{i,2} = table{i,2};
    end
    
    buffer2 = weight_list;
    
    for i = 1:2*nos-1
        parentlist{i} = '';
    end
    %disp(parentlist);
    for i = 1:nos-1
        buffer1 = buffer2(:,1);
        [probabilities,index] = sort(buffer1,'descend');
        summ  = probabilities(nos-i+1) + probabilities(nos-i);
        buffer2(nos+i,1) = summ;
        buffer2(index(nos-i+1),1) = 0;
        buffer2(index(nos-i),1) = 0;
    
        dict{nos+i,1} = strcat(dict{index(nos-i+1),1},dict{index(nos-i),1});
        dict{nos+i,2} = summ;
        parentlist{index(nos-i+1)} = nos+i;
        parentlist{index(nos-i)} = nos+i;
        %Huffman_tree{index(nos-i+1),6} = nos+i;
        %Huffman_tree{index(nos-i),6} = nos+i;
        dict{nos+i,3} = index(nos-i+1);
        dict{nos+i,4} = index(nos-i);
    end
    %disp(parentlist);
    %disp(Huffman_tree);
    for i = 1:2*nos-1
        child = i;
        parent = parentlist{i};
        while(parent~=0)
            if dict{parent,3} == child
                dict{i,5} = strcat('1',dict{i,5});
            else
                dict{i,5} = strcat('0',dict{i,5});
            end
    
            child = parent;
            parent = parentlist{parent};
        end
    end
end

function bin_seq = huffman_enc(sym_seq, dict)
    bin_seq = '';
    for i=1:strlength(sym_seq)
        curr = sym_seq(i);
        ind = find([dict{1:(size(dict,1)+1)/2,1:1}] == curr);
        %disp(ind);
        bin_seq = strcat(bin_seq,dict{ind,5});
    end
end

function sym_seq = huffman_dec(bin_seq, dict)
    sym_seq = '';
    currstart = 1;
    currend = 1;
    root = size(dict,1);
    currnode = root;
    while(currstart < strlength(bin_seq) && currend <= strlength(bin_seq))
        currstr = bin_seq(currstart:currend);
        if(currnode <= (size(dict,1)+1)/2)
            currstart = currend;
            %currend = currend;
            sym_seq = strcat(sym_seq,dict{currnode,1});
            currnode = root;
        else
            a = dict{dict{currnode,3},5};
            if(strcmp(a,currstr))
                if(currend == strlength(bin_seq))
                   b = dict{dict{currnode,3},1};
                   sym_seq = strcat(sym_seq,b);
                   break;
                end
                currnode = dict{currnode,3};
                currend  = currend + 1;
            else
                if(currend == strlength(bin_seq))
                   b = dict{dict{currnode,4},1};
                   sym_seq = strcat(sym_seq,b);
                   break;
                end
                currnode = dict{currnode,4};
                currend  = currend + 1;
            end
        end
    end
end


