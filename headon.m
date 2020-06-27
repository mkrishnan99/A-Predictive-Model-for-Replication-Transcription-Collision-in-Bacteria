%Head-on Collision Model
%In this model we have a simple linear dna strand, and we make several
%assumptions such as that once a collision occurs, the accessory items fall
%off and the replisome takes a certain amount of time to reinitiate. The
%model takes into consideration frequency of firing of the gene and its
%length other than the two obvious factors of the speeds of transcripion
%and replication. Unlike the codirectional case, here there is a
%possibility of multiple collisions even after the accessory items fall off
%after the first collision so we need to take that into account as well.
%The model will return total number of collisions and time taken by the
%replisome to cover the length of the replichore.
function [t,tot_time,n_col,pos_rep,j,pos_rep1]=headon(L,v_rep,v_trans,tr_start,Lt,t_resol,r_fire,t1)
tr_end= tr_start-Lt;
%Time taken to reach the end of the transcription unit from origin
t0= tr_end/v_rep;
%The replisome will collide with the first transcription unit which hasn't
%reached the end of the transciption unit. 
time_end=Lt/v_rep;
last_shot= floor((t1+time_end)/r_fire)+1; 
n_col = 0;
pos_rep=0;
j=0;
pos_rep1=0;
if last_shot==1
    n_col = 0;
else
    for i=1:1:last_shot
        if t1>=r_fire*(i-1) 
           if  v_trans*(t1-r_fire*(i-1))<=Lt %if first colliding shot is in gene when replisome reaches end
               shot =i;%Colliding shot
               pos = tr_start-v_trans*(t1-r_fire*(i-1));%position of colliding RNAP
               pos_r = tr_end;%pos of replisome
               n_col=1;
               break
            else
            n_col = 0;%If there are no transcription complexes present while the replisome passes there will be no collision
            end
        else
            pos=tr_start;
            pos_r=tr_end+v_rep*(r_fire*(i-1)-t1);
            t0=t0+r_fire*(i-1)-t1;
            t1=r_fire*(i-1);
            n_col=1;
            break
        end
    end
end
%If there is a collision then we need to consider the possibility of
%multiple collisions.
if n_col>0
    for time = 0:0.0005:((pos-pos_r)/v_rep)
        if (pos_r+v_rep*time)>= (pos-v_trans*time)
            t2= time+t_resol;
            time;
            pos_rep1 = v_rep*(t0+time);
            break
        else
            
        end
       
    end
    t2=time+t_resol;
    pos_rep1 = v_rep*(t0+time);
    T = t1+time;
    T1=t1+t2;
    next_shot= ceil(T/r_fire);%the next shot to be fired after everything is knocked off
    while tr_start-v_trans*(T1-r_fire*next_shot)<=pos_rep1%if its gone past the replisome
        next_shot = next_shot+1;
    end
    if T1<=(next_shot*r_fire)%if no RNAP has started after dwell time
        pos_rep= pos_rep1+v_rep*(next_shot*r_fire-t1-t2);
        pos_trans=tr_start;
        T1= next_shot*r_fire;
    else
        pos_rep=pos_rep1;
        pos_trans=tr_start-v_trans*(t1+t2-next_shot*r_fire);
    end
    
    while pos_rep< tr_start-0.5
          j=j+1;
        
         for time_next= 0:0.0001:((tr_start-pos_rep)/v_rep)
            if pos_rep + v_rep*time_next>= pos_trans-v_trans*time_next
                T = T+t_resol+time_next;
                t2 = t2+time_next+t_resol;
                next_shot= ceil(T/r_fire);
                n_col = n_col +1;
                pos_rep1= pos_rep + v_rep*time_next;
                T1=t1+t2;
                while tr_start-v_trans*(T1-r_fire*next_shot)<=pos_rep1
                      next_shot = next_shot+1;
                end
                    if T1<=(next_shot*r_fire)
                       pos_rep= pos_rep1+v_rep*(next_shot*r_fire-t1-t2);
                       pos_trans=tr_start;
                       t0=t0+next_shot*r_fire-t1-t2;
                       T1= next_shot*r_fire;
                       
                    else
                       pos_rep=pos_rep1;
                       pos_trans=tr_start-v_trans*(t1+t2-next_shot*r_fire);
                    end
                break
            else
                
            end
        end
    end
    if pos_rep>=tr_start-0.5
       t3=(tr_start-pos_rep1)/v_rep;
    end
    t_end = (L-tr_start)/v_rep;
    t = t0+t2+t3+t_end;
else
    t=(L)/v_rep;
end

tot_time= n_col*t_resol+(L/v_rep);

end