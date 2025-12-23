## Task3

### 伯努利分布

简而言之，一个事件只有两种可能，假如成功的概率值为$p$，失败则为则为$q(或是1-p)$，这是显而易见的，其归一性也可由$1-p+p=1$证得。

![](.\media\image34.GIF)

其概率质量函数可表示为：
$$
P(x)=p^x(1-p)^{1-x}=\begin{cases}
p	& \text{if x=1} \\
q	& \text{if x=0}
\end{cases}
$$
由于其是离散分布，其期望可表示为：
$$
E(X)=\sum{xP(x)}=0\times{q}+1\times{p}=p
$$
其方差可表示为：
$$
Var(X)=E((X-E(X))^2)=\sum(X-p)^2P(x)=(1-p)^2p+p^2q=pq^2+p^2q=pq
$$

### 二项分布

重复多次伯努利试验，便成了二项分布，成功和失败概率依旧用$p$和$q$表示，因为重复进行了多次，此时要讨论的问题是，进行了$n$次试验后，成功的次数到底为多少次？

这个可以由小学的组合数的概念得到，$n$次实验中成功了$r$次，那么我们可以把这$r$次任意的穿插到这$n$次中进行组合，因此，我们可以简单的通过组合数的公式，得到二项分布的概率密度为：
$$
P(X=m)={n\choose{m}}p^mq^{n-m}
$$
由二项式定理可证得其归一性：
$$
{n\choose{m}}p^mq^{n-m}=(p+q)^n=1
$$
![](.\media\image37.GIF)

其期望为：
$$
\begin{align}
E(X)&=\sum^n_{k=0}{n\choose{k}}p^kq^{n-k}=\sum^n_{k=0}n{{n-1}\choose{k-1}}p^kq^{n-k} \\
&=np\sum^n_{k=1}{{n-1}\choose{k-1}}p^{k-1}q^{n-k}	\\
&= np
\end{align}
$$
其方差为：
$$
Var(X)=E(X^2)-E^2(X)=n(n-1)p^2 + np- n^2p^2	=np(1-p)
$$
从上图可以看到，当抽样次数逐渐增加时，二项分布貌似逐渐向正态分布靠近。事实上，分布之间都有相互的联系，二项分布是超几何分布(不放回抽样)的极限情况，而泊松分布又是二项分布的极限情况，我们先来看看二项分布的极限情况,记$np=\lambda$：
$$
\begin{align}
{n\choose{k}}p^k(1-p)^{n-k}&=\frac{n(n-1)\dots(n-k+1)}{k!}(\frac{\lambda}{n})^k(1-\frac{\lambda}{n})^{n-k}	\\
&=\frac{\lambda^k}{k!}(1-\frac{1}{n})(1-\frac{2}{n})\cdots(1-\frac{k+1}{n})(1-\frac{\lambda}{n})^{n-k}	\\
\end{align}
$$
当$n\to\infty$时，可得$(1-\frac{1}{n})\cdots{(1-\frac{k+1}{n})}$中的每一项，皆为1，而由重要极限：$(1-\frac{\lambda}{n})^{n-k}=(1-\frac{\lambda}{n})^{\frac{n}{\lambda}(\lambda-\frac{\lambda{k}}n{})}=e^\lambda$

则有
$$
\lim_{n\to\infty}{n\choose{k}}p^k{1-p}^{n-k}=e^\lambda\frac{\lambda^k}{k!}
$$
这正是泊松分布的概率密度。

### Gamma分布和Beta分布

伽马分布和Beta分布与诸多分布都有关联：

称以下函数为伽马函数：
$$
\Gamma(\alpha)=\int^{\infty}_0x^{\alpha-1}e^{-x}dx
$$
这是一个特殊的函数，有：$\Gamma(n)=n!\;,\Gamma(\frac{1}{2})=\sqrt{\pi}$

伽马分布可表示为：
$$
p(x)=\frac{\lambda^n}{\Gamma(\alpha)}x^{\alpha-1}e^{-\lambda{x}}
$$
![](.\media\image74.GIF)

伽马分布表示为系统可以抵挡外来冲击，遇到第$k$次时失败，第$k$次冲击来到的时间$X$服从形状参数为$k$的伽马分布$Ga(k, \lambda)$ 。

特别的，有$Ga(1,\lambda)=Exp(\lambda)$，$Gamma(\frac{n}{2}, \frac{1}{2})=X^2(n)$，如果有了解过参数估计的知识，会了解到泊松分布的期望的共轭先验分布是伽马分布。

![](.\media\image76.png)

称以下函数为Beta函数：
$$
B(a, b)=\frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)} \quad{a>0,b>0}
$$
贝塔分布表示为：
$$
Be(a, b)=\frac{x^{a-1}(1-x)^{b-1}}{B(a,b)}\quad{0<x<1}
$$
![](.\media\image80.GIF)

用于不合格率、机器维修率、市场占有率等各种比率的情况。特别的，有$Be(1, 1)=U(0, 1)$。在后续学习参数估计中，会了解到伯努利试验中成功概率的共轭先验分布是贝塔分布。

![](.\media\image82.png)

### 中心极限定理

那他们和高斯分布又有什么关系呢？

用一个简单的高尔顿钉板感受一下：

![](./media/GB.gif)

实际上，任何分布，只要对它抽样足够多，其均值总会靠近正态分布，其证明一般通过特征函数的方法给出，此处暂不涉及。可由下图直观感受：

![](.\media\image16.GIF)

其公式化表达为：若$\{X_n\}$是独立同分布的随机变量序列，且$E(X_i)=\mu,Var(X_i)=\sigma^2>0$存在，记：
$$
Y^*_n=\frac{\sum^n_{i}X_i-n\mu}{\sigma\sqrt{n}}
$$
则$Y^*\to{N(0, 1)}$，即$\lim_{n\to\infty}P(Y^*_n\le{y})=\Phi(y)$

我们称之为，林德伯格-莱维中心极限定理

![](.\media\image18.GIF)

那么，对于二项分布而言也是一样：
$$
Y^*_n=\frac{S_n-np}{\sqrt{npq}}
$$
我们称之为涅莫夫-拉普拉斯中心极限定理

![](.\media\image21.GIF)