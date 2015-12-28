function  [Zsamp] = unif5(zlow,zhigh,k)


n=5;

z = zeros(n,k);
for i = 1:n
    for j = 1:k
        z(i,j) = zlow(i) + ((j-1)/(k-1))*(zhigh(i)-zlow(i));
    end
end

Zsamp = zeros(n,k^n);

cnt = 0;
for a = 1:k
    for b = 1:k
        for c = 1:k
            for d = 1:k
                for e = 1:k
                    cnt = cnt+1;
                    Zsamp(1,cnt) = z(i,a);
                    Zsamp(2,cnt) = z(i,b);
                    Zsamp(3,cnt) = z(i,c);
                    Zsamp(4,cnt) = z(i,d);
                    Zsamp(5,cnt) = z(i,e);
                end
            end
        end
    end
end

end