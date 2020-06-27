%Codirectional Collision Model
%In this model we have a simple linear dna strand, and we make several
%assumptions such as that once a collision occurs, the accessory items fall
%off and the replisome takes a certain amount of time to reinitiate. The
%model takes into consideration frequency of firing of the gene and its
%length other than the two obvious factors of the speeds of transcripion
%and replication. Since the transcription speed is lower than that of replication
%once a collision occurs, there cannot be a second one. It gives us 
%at the end the number of collisions that
%occurs and the total time taken for the replisome to go the full length of
%the replichore.
function [t,tot_time,n_col]=codirect(L,v_rep,v_trans,tr_start,Lt,t_resol,r_fire,t1)
tr_end= tr_start+Lt;
t0= tr_start/v_rep;

last_shot= floor(t1/r_fire);
if last_shot == 0
    pos=tr_start+v_trans*t1;
else
    pos = tr_start+ v_trans*(rem(t1,r_fire*last_shot));
end
if pos<=tr_end
    for time = 0:0.0001:((tr_end-pos)/v_trans)
        if (tr_start+v_rep*time)>= (pos+v_trans*time)
            t2= time+t_resol;
            time;
            n_col=1;
            t3 = (L-(tr_start+v_rep*time))/v_rep;
            break
        else
            t2= Lt/v_rep;
            n_col = 0;
            t3= (L-tr_end)/v_rep;
        end
    end
else
    t2= Lt/v_rep;
    n_col = 0;
    t3= (L-tr_end)/v_rep;
end
    
t = t0+t2+t3;
n_col;
tot_time= n_col*t_resol+(L/v_rep);
end
        