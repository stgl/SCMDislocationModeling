function centroids=PatchCentroid(c,v)

centroids=zeros(size(v));

for i=1:length(v)
    
    i1=v(i,1);
    i2=v(i,2);
    i3=v(i,3);
    
    centroids(i,1)=(c(i1,1)+c(i2,1)+c(i3,1))/3;
    centroids(i,2)=(c(i1,2)+c(i2,2)+c(i3,2))/3;
    centroids(i,3)=(c(i1,3)+c(i2,3)+c(i3,3))/3;
        
end

end





