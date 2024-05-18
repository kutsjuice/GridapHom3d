N = 10000000
k = zeros(Threads.nthreads())
using Polyester
norm(a, b) = sqrt(a^2 + b^2)
@batch for i in 1:N
    x, y, z = rand(3)
    if norm(x, y) < 1 && norm(x, z) < 1 && norm(z, y) < 1
        k[Threads.threadid()]+=1
    end
end
print(sum(k)/N)
2 - sqrt(2)