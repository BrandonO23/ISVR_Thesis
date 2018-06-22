function sphplot(Cs,dpos,bpos)
scatter3(dpos(:,1),dpos(:,2),dpos(:,3),'k','fill'),hold on
scatter3(bpos(:,1),bpos(:,2),bpos(:,3),'k')
scatter3(Cs(:,1),Cs(:,2),Cs(:,3),150,'*','r')
end

