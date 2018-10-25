using PyPlot

const L=64; const J=1.0; const h=0.0; const N=L*L
const MCSPS=80000;
const M=Ref{Float64}(0.0)
const E=Ref{Float64}(0.0)

function initialize()
    s=[Int8(rand(0:1)==0?-1:1) for i in 1:L^2]
    M[]=sum(s)/N
    s=reshape(s,L,L)
    for i in 1:L
        for j in 1:L
            sp=s[i,j]
            E[]-=(J*sp*(s[(i+L-2)%L+1,j]+s[(i+L)%L+1,j]+s[i,(j+L-2)%L+1]+s[i,(j+L)%L+1]))/2.0+h*sp
        end
    end
    E[]/=N
    return s
end       

function update(temp::Float64,s::Array{Int8})
    # Metropolis
    # update for a monte carlo step per sites (a sweep)
    for step in 1:N
        i, j=rand(1:L),rand(1:L)
        sp=s[i,j]
        delE=2*J*sp*(s[(i+L-2)%L+1,j]+s[(i+L)%L+1,j]+s[i,(j+L-2)%L+1]+s[i,(j+L)%L+1])+2*h*sp
        if (rand()<exp(-delE/temp))
            s[i,j]=-sp
            E[]+=delE/N
            M[]-=2sp/N
        end
    end
end

function main()
    #temps=linspace(3, 0.1, 30)
    temps=linspace(2.3, 2.2, 21)
    values=Array{Float64}(MCSPS,length(temps)*2)
    s=initialize()

    for t in 1:length(temps)
        temp=temps[t]
        for i in 1:MCSPS
            update(temp,s)
            values[i,2*t-1]=M[]
            values[i,2*t]=E[]
        end
        println(temp)
    end
    writedlm(@sprintf("CMPL%dM%d",L,MCSPS),values)
end
tic()
main()
toc()
