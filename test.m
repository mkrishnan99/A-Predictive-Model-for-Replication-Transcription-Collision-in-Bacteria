tdfread('Bacillus_subtilis_rep1.txt','\t')
L0 = input('Length of replichore');
v_rep = input('Speed of replisome');
v_trans= input('Speed of transcription');
t_resol = input('Enter time taken to resolve collision');

[row,colm]=size(Class);
head='"Lagging"';
a0=0;
b0=0;
t0=zeros(row,1);
n_col0=zeros(row,1);
pos_0=zeros(row,1);
pos1_0=zeros(row,1);
J = zeros(row,1);
for i =1:1:row
    flag=strcmp(head,Status(i,:));
    L=End(i)-a0;
    t1=b0+((Start(i)-a0)/v_rep);
    Lt=End(i)-Start(i);
    r_fire=Br_fire(i);
    if flag==1
        tr_start= End(i)-a0;
       [t0(i,1),tot_time0(i,1),n_col0(i,1),pos_0(i,1), J(i,1),pos1_0(i,1)]=headon(L,v_rep,v_trans,tr_start,Lt,t_resol,r_fire,t1);
       
    else
        tr_start=Start(i)-a0;
       [t0(i,1),tot_time0(i,1),n_col0(i,1)]= codirect(L,v_rep,v_trans,tr_start,Lt,t_resol,r_fire,t1);
    end
    a0=End(i);
    b0=t0(i,1);
 
end
t = sum(t0)+((L0-End(row))/v_rep)
tot_time  = sum(tot_time0)+((L0-End(row))/v_rep)
n_col=sum(n_col0)
