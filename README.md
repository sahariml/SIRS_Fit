###Fit the parameters of discrete SIRS model
##\begin{equation}
\begin{cases}
S_{n+1}=F(S_n,I_n,R_n,):=S_{n}+h(A-\mu S_{n}-\lambda S_{n}I_{n}+\beta R_{n}),\\
I_{n+1}=G(S_n,I_n,R_n,):=I_{n}+h(\lambda S_{n}I_{n}-(\mu+r)I_{n}),\\
R_{n+1}=H(S_n,I_n,R_n,):=R_{n}+h(rI_{n}-(\mu+\beta)R_{n}),\\
N_{n+1}=(1-h\mu)N_{n}+hA,
\end{cases}
\end{equation}
