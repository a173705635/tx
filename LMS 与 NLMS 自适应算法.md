## 27.1 本章要解决的问题

上一章中，GSC 辅助支路的最优维纳解为：

$$
\mathbf g_{\mathrm{opt}}
=
\mathbf R_u^{-1}\mathbf p_{ud}
$$

其中：

$$
\mathbf R_u
=
E[\mathbf u[k]\mathbf u^H[k]]
$$

$$
\mathbf p_{ud}
=
E[\mathbf u[k]d^*[k]]
$$

这个闭式解在理论上很清楚，但实际使用时存在几个问题：

1. 真实的统计期望未知；
2. 只能通过有限快拍估计协方差和互相关；
3. 需要求矩阵逆；
4. 当干扰环境变化时，需要重新估计和求解；
5. 数据维数较大时，计算量明显增加。

LMS 和 NLMS 的作用，就是不用显式求逆，而是利用每个新样本逐步更新权值，使其逐渐逼近维纳最优解。

本章重点回答以下问题：

- LMS 是怎样从维纳解和最速下降法推导出来的？
- 为什么复数 LMS 的更新式中出现共轭误差？
- 步长 $\mu$ 的物理意义是什么？
- LMS 为什么会收敛到维纳解附近，而不是严格停在维纳解？
- 输入信号功率变化为什么会影响 LMS 稳定性？
- NLMS 为什么要除以输入向量能量？
- NLMS 的更新公式能否通过优化问题严格推导？
- 在 GSC 中，LMS 和 NLMS 实际在学习什么？
- 阻塞矩阵存在目标泄漏时会发生什么？

---

## 27.2 GSC 中的基本信号定义

GSC 的总权值写成：

$$
\mathbf w[k]
=
\mathbf w_q
-
\mathbf B\mathbf g[k]
$$

其中：

- $\mathbf w_q$：固定权值，负责满足原约束；
- $\mathbf B$：阻塞矩阵，满足 $\mathbf C^H\mathbf B=\mathbf0$；
- $\mathbf g[k]$：辅助支路自适应权值。

固定支路输出为：

$$
d[k]
=
\mathbf w_q^H\mathbf x[k]
$$

阻塞支路输出为：

$$
\mathbf u[k]
=
\mathbf B^H\mathbf x[k]
$$

辅助滤波器输出为：

$$
\widehat{i}[k]
=
\mathbf g^H[k]\mathbf u[k]
$$

GSC 最终输出为：

$$
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
$$

这里使用 $e[k]$ 表示误差信号，同时它也是 GSC 的最终输出。

理想情况下：

$$
\mathbf B^H\mathbf a_0=\mathbf0
$$

所以目标方向信号不进入辅助支路。此时 $\mathbf u[k]$ 主要包含干扰和噪声，辅助滤波器通过 $\mathbf g^H[k]\mathbf u[k]$ 估计固定支路中的干扰成分，再从 $d[k]$ 中减去。

---

## 27.3 维纳滤波问题回顾

定义均方误差代价函数：

$$
J(\mathbf g)
=
E[|e[k]|^2]
$$

其中：

$$
e[k]
=
d[k]
-
\mathbf g^H\mathbf u[k]
$$

省略时间下标，写成：

$$
e=d-\mathbf g^H\mathbf u
$$

则：

$$
|e|^2
=
ee^*
$$

代入：

$$
|e|^2
=
\left(
d-\mathbf g^H\mathbf u
\right)
\left(
d^*-\mathbf u^H\mathbf g
\right)
$$

展开：

$$
|e|^2
=
|d|^2
-
d\mathbf u^H\mathbf g
-
\mathbf g^H\mathbf u d^*
+
\mathbf g^H\mathbf u\mathbf u^H\mathbf g
$$

取期望：

$$
J(\mathbf g)
=
E[|d|^2]
-
E[d\mathbf u^H]\mathbf g
-
\mathbf g^HE[\mathbf u d^*]
+
\mathbf g^HE[\mathbf u\mathbf u^H]\mathbf g
$$

定义：

$$
\sigma_d^2
=
E[|d|^2]
$$

$$
\mathbf R_u
=
E[\mathbf u\mathbf u^H]
$$

$$
\mathbf p_{ud}
=
E[\mathbf u d^*]
$$

并且：

$$
E[d\mathbf u^H]
=
\mathbf p_{ud}^H
$$

所以：

$$
J(\mathbf g)
=
\sigma_d^2
-
\mathbf p_{ud}^H\mathbf g
-
\mathbf g^H\mathbf p_{ud}
+
\mathbf g^H\mathbf R_u\mathbf g
$$

---

## 27.4 对代价函数求梯度

由于 $\mathbf g$ 是复向量，使用 Wirtinger 微分，对 $\mathbf g^*$ 求导。

第一项：

$$
\frac{\partial \sigma_d^2}
{\partial\mathbf g^*}
=
\mathbf0
$$

第二项：

$$
-\mathbf p_{ud}^H\mathbf g
$$

只含 $\mathbf g$，不含 $\mathbf g^*$，所以：

$$
\frac{\partial}
{\partial\mathbf g^*}
\left(
-\mathbf p_{ud}^H\mathbf g
\right)
=
\mathbf0
$$

第三项：

$$
-\mathbf g^H\mathbf p_{ud}
$$

因为：

$$
\mathbf g^H=(\mathbf g^*)^T
$$

所以：

$$
\frac{\partial}
{\partial\mathbf g^*}
\left(
-\mathbf g^H\mathbf p_{ud}
\right)
=
-\mathbf p_{ud}
$$

第四项：

$$
\mathbf g^H\mathbf R_u\mathbf g
$$

求导得到：

$$
\frac{\partial}
{\partial\mathbf g^*}
\left(
\mathbf g^H\mathbf R_u\mathbf g
\right)
=
\mathbf R_u\mathbf g
$$

因此：

$$
\nabla_{\mathbf g^*}J
=
-\mathbf p_{ud}
+
\mathbf R_u\mathbf g
$$

也可以写成：

$$
\nabla_{\mathbf g^*}J
=
\mathbf R_u\mathbf g
-
\mathbf p_{ud}
$$

令梯度为零：

$$
\mathbf R_u\mathbf g_{\mathrm{opt}}
-
\mathbf p_{ud}
=
\mathbf0
$$

于是：

$$
\mathbf R_u\mathbf g_{\mathrm{opt}}
=
\mathbf p_{ud}
$$

得到维纳解：

$$
\mathbf g_{\mathrm{opt}}
=
\mathbf R_u^{-1}\mathbf p_{ud}
$$

---

## 27.5 为什么需要迭代算法

维纳解要求知道：

$$
\mathbf R_u
=
E[\mathbf u\mathbf u^H]
$$

以及：

$$
\mathbf p_{ud}
=
E[\mathbf u d^*]
$$

但实际只能获得有限样本：

$$
\mathbf u[1],\mathbf u[2],\ldots
$$

$$
d[1],d[2],\ldots
$$

而且环境可能随时间变化。

因此一种自然想法是：

> 不直接求 $\mathbf R_u^{-1}\mathbf p_{ud}$，而是沿代价函数下降方向逐步调整 $\mathbf g$。

这就是最速下降法的起点。

---

## 27.6 最速下降法

对于一般实值代价函数 $J(\mathbf g)$，沿负梯度方向更新：

$$
\mathbf g[k+1]
=
\mathbf g[k]
-
\mu
\nabla_{\mathbf g^*}J
\big|_{\mathbf g=\mathbf g[k]}
$$

其中：

$$
\mu>0
$$

称为步长。

代入维纳问题的梯度：

$$
\nabla_{\mathbf g^*}J
=
\mathbf R_u\mathbf g
-
\mathbf p_{ud}
$$

得到：

$$
\mathbf g[k+1]
=
\mathbf g[k]
-
\mu
\left(
\mathbf R_u\mathbf g[k]
-
\mathbf p_{ud}
\right)
$$

即：

$$
\boxed{
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu
\left(
\mathbf p_{ud}
-
\mathbf R_u\mathbf g[k]
\right)
}
$$

这就是复数最速下降算法。

---

## 27.7 最速下降法为什么会走向维纳解

如果迭代已经到达稳定点，则：

$$
\mathbf g[k+1]
=
\mathbf g[k]
=
\mathbf g_\infty
$$

代入更新式：

$$
\mathbf g_\infty
=
\mathbf g_\infty
+
\mu
\left(
\mathbf p_{ud}
-
\mathbf R_u\mathbf g_\infty
\right)
$$

所以：

$$
\mathbf p_{ud}
-
\mathbf R_u\mathbf g_\infty
=
\mathbf0
$$

因此：

$$
\mathbf R_u\mathbf g_\infty
=
\mathbf p_{ud}
$$

得到：

$$
\mathbf g_\infty
=
\mathbf R_u^{-1}\mathbf p_{ud}
$$

也就是：

$$
\mathbf g_\infty
=
\mathbf g_{\mathrm{opt}}
$$

因此，只要步长合适并且统计量固定，最速下降法的稳定点就是维纳解。

---

## 27.8 最速下降法仍然存在的问题

最速下降法虽然不需要直接求逆，但仍然需要精确知道：

$$
\mathbf R_u
$$

和：

$$
\mathbf p_{ud}
$$

实际中这些期望仍然未知。

LMS 的关键变化是：

> 用当前样本构造的瞬时梯度，近似真实的统计平均梯度。

---

## 27.9 瞬时代价函数

真实代价函数为：

$$
J(\mathbf g)
=
E[|e[k]|^2]
$$

LMS 不直接使用期望，而是定义当前时刻的瞬时代价：

$$
\widehat{J}[k]
=
|e[k]|^2
$$

其中：

$$
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
$$

瞬时代价不是长期平均性能，但它可以由当前样本立即计算。

LMS 使用：

$$
\nabla_{\mathbf g^*}\widehat{J}[k]
$$

代替：

$$
\nabla_{\mathbf g^*}J
$$

---

## 27.10 复数 LMS 梯度的完整推导

瞬时误差：

$$
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
$$

其共轭为：

$$
e^*[k]
=
d^*[k]
-
\mathbf u^H[k]\mathbf g[k]
$$

瞬时代价：

$$
\widehat{J}[k]
=
e[k]e^*[k]
$$

对 $\mathbf g^*[k]$ 求导。

由乘积法则：

$$
\frac{\partial \widehat{J}[k]}
{\partial\mathbf g^*[k]}
=
\frac{\partial e[k]}
{\partial\mathbf g^*[k]}
e^*[k]
+
e[k]
\frac{\partial e^*[k]}
{\partial\mathbf g^*[k]}
$$

先看第一项。

因为：

$$
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
$$

且：

$$
\mathbf g^H[k]
=
(\mathbf g^*[k])^T
$$

所以：

$$
\frac{\partial e[k]}
{\partial\mathbf g^*[k]}
=
-\mathbf u[k]
$$

再看第二项。

$$
e^*[k]
=
d^*[k]
-
\mathbf u^H[k]\mathbf g[k]
$$

它只含 $\mathbf g[k]$，不含 $\mathbf g^*[k]$。

在 Wirtinger 微分中，把 $\mathbf g$ 和 $\mathbf g^*$ 形式上视作独立变量，因此：

$$
\frac{\partial e^*[k]}
{\partial\mathbf g^*[k]}
=
\mathbf0
$$

于是：

$$
\frac{\partial \widehat{J}[k]}
{\partial\mathbf g^*[k]}
=
-\mathbf u[k]e^*[k]
$$

即：

$$
\boxed{
\nabla_{\mathbf g^*}\widehat{J}[k]
=
-\mathbf u[k]e^*[k]
}
$$

---

## 27.11 LMS 更新公式

沿瞬时负梯度方向更新：

$$
\mathbf g[k+1]
=
\mathbf g[k]
-
\mu
\nabla_{\mathbf g^*}\widehat{J}[k]
$$

代入：

$$
\nabla_{\mathbf g^*}\widehat{J}[k]
=
-\mathbf u[k]e^*[k]
$$

得到：

$$
\boxed{
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
}
$$

同时：

$$
\boxed{
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
}
$$

这就是复数 LMS 的基本更新式。

---

## 27.12 为什么更新式中是 $e^*[k]$

这是复数 LMS 中最容易混淆的一点。

代价函数是：

$$
|e[k]|^2
=
e[k]e^*[k]
$$

对 $\mathbf g^*$ 求导时：

$$
\frac{\partial e[k]}
{\partial\mathbf g^*}
=
-\mathbf u[k]
$$

而：

$$
\frac{\partial e^*[k]}
{\partial\mathbf g^*}
=
\mathbf0
$$

所以剩余的乘法因子是：

$$
e^*[k]
$$

最终：

$$
\nabla_{\mathbf g^*}|e[k]|^2
=
-\mathbf u[k]e^*[k]
$$

因此更新中自然出现：

$$
\mathbf u[k]e^*[k]
$$

不是人为规定，也不是为了保证维度，而是由复数梯度严格推导得到。

---

## 27.13 LMS 与统计梯度之间的关系

真实梯度：

$$
\nabla_{\mathbf g^*}J
=
\mathbf R_u\mathbf g
-
\mathbf p_{ud}
$$

瞬时负梯度方向为：

$$
\mathbf u[k]e^*[k]
$$

计算它的期望：

$$
E[\mathbf u e^*]
$$

由于：

$$
e^*
=
d^*
-
\mathbf u^H\mathbf g
$$

所以：

$$
E[\mathbf u e^*]
=
E[\mathbf u d^*]
-
E[\mathbf u\mathbf u^H]\mathbf g
$$

即：

$$
E[\mathbf u e^*]
=
\mathbf p_{ud}
-
\mathbf R_u\mathbf g
$$

因此：

$$
E[\mathbf u e^*]
=
-
\nabla_{\mathbf g^*}J
$$

这说明：

> LMS 使用的瞬时方向 $\mathbf u[k]e^*[k]$，在平均意义上正是真实负梯度方向。

所以 LMS 虽然每一步带有随机波动，但总体趋势仍然朝向维纳解。

---

## 27.14 LMS 为什么不会严格停在维纳解

最速下降法使用精确统计梯度：

$$
\mathbf p_{ud}
-
\mathbf R_u\mathbf g[k]
$$

当到达维纳解时，该梯度严格为零。

LMS 使用的是随机瞬时梯度：

$$
\mathbf u[k]e^*[k]
$$

即使：

$$
\mathbf g[k]
=
\mathbf g_{\mathrm{opt}}
$$

当前样本仍然可能满足：

$$
\mathbf u[k]e^*[k]
\neq
\mathbf0
$$

只是它的期望为零：

$$
E[\mathbf u[k]e^*[k]]
=
\mathbf0
$$

因此 LMS 到达最优点附近后，权值仍会在附近随机波动。

这种现象称为：

- 稳态失调；
- steady-state misadjustment；
- 权值抖动。

步长越大，这种波动通常越明显。

---

## 27.15 步长 $\mu$ 的作用

LMS 更新为：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

步长 $\mu$ 决定每次更新幅度。

### 步长较大

优点：

- 初期收敛更快；
- 环境变化后跟踪速度更快。

缺点：

- 稳态波动更大；
- 输出剩余误差更大；
- 可能不稳定甚至发散。

### 步长较小

优点：

- 稳态更加平稳；
- 失调较小；
- 数值安全性更高。

缺点：

- 收敛速度慢；
- 环境变化后跟踪能力弱。

所以 LMS 中始终存在一个基本权衡：

$$
\text{收敛速度}
\quad\leftrightarrow\quad
\text{稳态误差}
$$

---

## 27.16 LMS 的均值收敛分析

定义最优权值：

$$
\mathbf g_{\mathrm{opt}}
=
\mathbf R_u^{-1}\mathbf p_{ud}
$$

定义权值误差向量：

$$
\widetilde{\mathbf g}[k]
=
\mathbf g_{\mathrm{opt}}
-
\mathbf g[k]
$$

LMS 更新：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

两边从 $\mathbf g_{\mathrm{opt}}$ 中相减：

$$
\widetilde{\mathbf g}[k+1]
=
\widetilde{\mathbf g}[k]
-
\mu\mathbf u[k]e^*[k]
$$

为了分析均值，采用自适应滤波中常见的独立性近似，可以得到：

$$
E[\widetilde{\mathbf g}[k+1]]
=
\left(
\mathbf I-
\mu\mathbf R_u
\right)
E[\widetilde{\mathbf g}[k]]
$$

递推后：

$$
E[\widetilde{\mathbf g}[k]]
=
\left(
\mathbf I-
\mu\mathbf R_u
\right)^k
E[\widetilde{\mathbf g}[0]]
$$

要使均值误差趋于零，矩阵：

$$
\mathbf I-
\mu\mathbf R_u
$$

的所有特征值模长必须小于 1。

---

## 27.17 LMS 步长稳定条件的推导

设 $\mathbf R_u$ 的特征值为：

$$
\lambda_1,\lambda_2,\ldots,\lambda_M
$$

其中：

$$
M=N-r
$$

则：

$$
\mathbf I-
\mu\mathbf R_u
$$

的特征值为：

$$
1-
\mu\lambda_i
$$

收敛要求：

$$
|1-
\mu\lambda_i|<1
$$

对于每个 $\lambda_i>0$：

$$
-1
<
1-
\mu\lambda_i
<
1
$$

右侧不等式：

$$
1-
\mu\lambda_i<1
$$

得到：

$$
\mu>0
$$

左侧不等式：

$$
1-
\mu\lambda_i>-1
$$

得到：

$$
\mu\lambda_i<2
$$

即：

$$
\mu<\frac{2}{\lambda_i}
$$

要对所有特征值成立，最严格的是最大特征值：

$$
\boxed{
0<\mu<\frac{2}{\lambda_{\max}(\mathbf R_u)}
}
$$

这是 LMS 均值收敛的经典步长条件。

实际使用中通常还要留出裕量，不会把 $\mu$ 取到理论上限附近。

---

## 27.18 为什么输入功率影响 LMS 稳定性

如果输入整体放大：

$$
\mathbf u'[k]
=
c\mathbf u[k]
$$

那么协方差矩阵变为：

$$
\mathbf R_u'
=
E[\mathbf u'\mathbf u'^H]
$$

所以：

$$
\mathbf R_u'
=
|c|^2\mathbf R_u
$$

最大特征值变为：

$$
\lambda_{\max}(\mathbf R_u')
=
|c|^2\lambda_{\max}(\mathbf R_u)
$$

于是稳定步长上限变为：

$$
\mu
<
\frac{2}
{|c|^2\lambda_{\max}(\mathbf R_u)}
$$

也就是说：

> 输入信号功率越大，LMS 允许的固定步长越小。

这正是 NLMS 出现的直接原因。

---

## 27.19 特征值扩展为什么影响收敛速度

将协方差矩阵特征值按大小排列：

$$
\lambda_{\max}
\ge
\cdots
\ge
\lambda_{\min}
>0
$$

不同特征方向上的均值误差衰减因子为：

$$
1-
\mu\lambda_i
$$

大特征值方向通常变化快，小特征值方向变化慢。

定义特征值扩展：

$$
\kappa
=
\frac{\lambda_{\max}}
{\lambda_{\min}}
$$

当：

$$
\kappa
$$

很大时，不同方向的收敛速度差异明显。

为了保证最大特征值方向稳定，$\mu$ 不能太大；但这样一来，最小特征值方向收敛会非常慢。

因此：

> LMS 的收敛速度不仅取决于输入总功率，还取决于输入协方差矩阵的特征值分布。

高度相关的辅助输入常常会造成较大的特征值扩展。

---

## 27.20 从几何上理解 LMS

代价函数：

$$
J(\mathbf g)
=
\sigma_d^2
-
\mathbf p_{ud}^H\mathbf g
-
\mathbf g^H\mathbf p_{ud}
+
\mathbf g^H\mathbf R_u\mathbf g
$$

对于正定 $\mathbf R_u$，它在权值空间中形成一个碗状二次曲面。

- 若特征值接近，等高线接近圆形；
- 若特征值差异很大，等高线是狭长椭圆。

最速下降和 LMS 都沿负梯度方向前进。

在狭长椭圆中，负梯度方向常常横跨谷底，导致权值在两侧来回摆动，再缓慢向最优点前进。

这就是输入相关性大时 LMS 收敛较慢的几何原因。

---

## 27.21 NLMS 的基本思想

LMS 更新量为：

$$
\Delta\mathbf g[k]
=
\mu\mathbf u[k]e^*[k]
$$

其大小与：

$$
\|\mathbf u[k]\|
$$

直接相关。

当输入瞬时能量很大时，更新量会突然变大，可能造成不稳定。

当输入瞬时能量很小时，更新量又可能过小，导致收敛缓慢。

NLMS 的核心思想是：

> 根据当前输入向量的能量，对更新量进行归一化。

NLMS 更新式为：

$$
\boxed{
\mathbf g[k+1]
=
\mathbf g[k]
+
\frac{\mu}
{\|\mathbf u[k]\|^2+\varepsilon}
\mathbf u[k]e^*[k]
}
$$

其中：

$$
\|\mathbf u[k]\|^2
=
\mathbf u^H[k]\mathbf u[k]
$$

而：

$$
\varepsilon>0
$$

是防止分母过小的正则项。

---

## 27.22 NLMS 不是简单的经验除法

除以：

$$
\|\mathbf u[k]\|^2
$$

并不是随意的经验技巧。

NLMS 可以从一个明确的最小扰动优化问题推导出来：

> 在使当前样本误差被校正的前提下，让本次权值改变量尽可能小。

下面完整推导。

---

## 27.23 NLMS 的最小扰动推导

当前权值为：

$$
\mathbf g[k]
$$

更新后权值为：

$$
\mathbf g[k+1]
$$

希望本次改变量最小：

$$
\min_{\mathbf g[k+1]}
\left\|
\mathbf g[k+1]-\mathbf g[k]
\right\|^2
$$

同时，希望更新后的权值对当前样本满足零误差约束：

$$
d[k]
-
\mathbf g^H[k+1]\mathbf u[k]
=
0
$$

定义：

$$
\Delta\mathbf g[k]
=
\mathbf g[k+1]-\mathbf g[k]
$$

那么：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\Delta\mathbf g[k]
$$

约束变为：

$$
d[k]
-
\left(
\mathbf g[k]+\Delta\mathbf g[k]
\right)^H
\mathbf u[k]
=
0
$$

展开：

$$
d[k]
-
\mathbf g^H[k]\mathbf u[k]
-
\Delta\mathbf g^H[k]\mathbf u[k]
=
0
$$

而当前误差为：

$$
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
$$

所以：

$$
e[k]
-
\Delta\mathbf g^H[k]\mathbf u[k]
=
0
$$

即：

$$
\Delta\mathbf g^H[k]\mathbf u[k]
=
e[k]
$$

取共轭：

$$
\mathbf u^H[k]\Delta\mathbf g[k]
=
e^*[k]
$$

因此优化问题变为：

$$
\min_{\Delta\mathbf g}
\Delta\mathbf g^H\Delta\mathbf g
$$

满足：

$$
\mathbf u^H\Delta\mathbf g
=
e^*
$$

---

## 27.24 使用拉格朗日乘子求最小改变量

构造实值拉格朗日函数：

$$
\mathcal L
=
\Delta\mathbf g^H\Delta\mathbf g
+
\lambda^*
\left(
\mathbf u^H\Delta\mathbf g-e^*
\right)
+
\lambda
\left(
\Delta\mathbf g^H\mathbf u-e
\right)
$$

对：

$$
\Delta\mathbf g^*
$$

求导。

第一项：

$$
\frac{\partial}
{\partial\Delta\mathbf g^*}
\left(
\Delta\mathbf g^H\Delta\mathbf g
\right)
=
\Delta\mathbf g
$$

第二项：

$$
\lambda^*
\left(
\mathbf u^H\Delta\mathbf g-e^*
\right)
$$

只含 $\Delta\mathbf g$，不含 $\Delta\mathbf g^*$，所以导数为零。

第三项：

$$
\lambda
\left(
\Delta\mathbf g^H\mathbf u-e
\right)
$$

对 $\Delta\mathbf g^*$ 求导：

$$
\lambda\mathbf u
$$

因此：

$$
\frac{\partial\mathcal L}
{\partial\Delta\mathbf g^*}
=
\Delta\mathbf g
+
\lambda\mathbf u
$$

令其为零：

$$
\Delta\mathbf g
+
\lambda\mathbf u
=
\mathbf0
$$

所以：

$$
\Delta\mathbf g
=
-\lambda\mathbf u
$$

代入约束：

$$
\mathbf u^H\Delta\mathbf g
=
e^*
$$

得到：

$$
-\lambda\mathbf u^H\mathbf u
=
e^*
$$

因此：

$$
\lambda
=
-\frac{e^*}
{\mathbf u^H\mathbf u}
$$

代回：

$$
\Delta\mathbf g
=
\frac{\mathbf u e^*}
{\mathbf u^H\mathbf u}
$$

于是：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\frac{
\mathbf u[k]e^*[k]
}{
\mathbf u^H[k]\mathbf u[k]
}
$$

这对应一次完全校正当前样本误差的更新。

---

## 27.25 为什么还要引入步长 $\mu$

刚才推导得到的更新相当于：

$$
\mu=1
$$

它试图让更新后的权值严格满足当前样本零误差。

但当前样本中含有噪声，若每次都完全拟合当前样本，会导致较大波动。

因此加入松弛系数：

$$
\mu
$$

得到：

$$
\boxed{
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu
\frac{
\mathbf u[k]e^*[k]
}{
\mathbf u^H[k]\mathbf u[k]
}
}
$$

实际中，为防止：

$$
\mathbf u^H[k]\mathbf u[k]
$$

接近零，再加入：

$$
\varepsilon
$$

得到标准 NLMS：

$$
\boxed{
\mathbf g[k+1]
=
\mathbf g[k]
+
\frac{\mu}
{\mathbf u^H[k]\mathbf u[k]+\varepsilon}
\mathbf u[k]e^*[k]
}
$$

---

## 27.26 NLMS 的另一个直观推导

当前误差：

$$
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
$$

假设更新：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\alpha[k]\mathbf u[k]e^*[k]
$$

其中 $\alpha[k]$ 是当前时刻的有效步长。

更新后的同一样本误差为：

$$
e_{\mathrm{new}}[k]
=
d[k]
-
\mathbf g^H[k+1]\mathbf u[k]
$$

代入更新式：

$$
\mathbf g^H[k+1]
=
\mathbf g^H[k]
+
\alpha[k]e[k]\mathbf u^H[k]
$$

于是：

$$
e_{\mathrm{new}}[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
-
\alpha[k]e[k]\mathbf u^H[k]\mathbf u[k]
$$

前两项就是 $e[k]$：

$$
e_{\mathrm{new}}[k]
=
e[k]
\left(
1-
\alpha[k]\|\mathbf u[k]\|^2
\right)
$$

如果希望一次更新使当前误差归零：

$$
e_{\mathrm{new}}[k]=0
$$

则：

$$
1-
\alpha[k]\|\mathbf u[k]\|^2
=
0
$$

所以：

$$
\alpha[k]
=
\frac1{\|\mathbf u[k]\|^2}
$$

加入松弛步长 $\mu$ 后：

$$
\alpha[k]
=
\frac{\mu}
{\|\mathbf u[k]\|^2}
$$

因此再次得到 NLMS 更新式。

---

## 27.27 NLMS 步长范围

忽略 $\varepsilon$，观察更新后当前样本误差：

$$
e_{\mathrm{new}}[k]
=
e[k](1-
\mu)
$$

因此：

- $\mu=0$：完全不更新；
- $0<\mu<1$：误差同号缩小；
- $\mu=1$：理想条件下当前样本误差被一次消除；
- $1<\mu<2$：出现过校正，误差反号但幅度仍缩小；
- $\mu=2$：误差幅度不变，只反号；
- $\mu>2$：误差幅度增大，趋于不稳定。

因此标准 NLMS 常用理论范围为：

$$
\boxed{
0<\mu<2
}
$$

实际中通常选取比 2 小得多的值，以降低噪声引起的稳态波动。

加入 $\varepsilon$ 后，严格瞬时关系会有所变化，但上述范围仍是常用设计依据。

---

## 27.28 $\varepsilon$ 的作用

NLMS 分母为：

$$
\mathbf u^H[k]\mathbf u[k]+\varepsilon
$$

如果当前输入几乎为零：

$$
\mathbf u^H[k]\mathbf u[k]\approx0
$$

没有 $\varepsilon$ 时，等效步长会变得非常大：

$$
\frac{\mu}
{\mathbf u^H[k]\mathbf u[k]}
$$

可能造成数值不稳定。

因此加入：

$$
\varepsilon>0
$$

其作用包括：

1. 防止除零；
2. 限制低输入能量时的最大更新量；
3. 提供一定正则化；
4. 改善有限精度实现的稳定性。

但 $\varepsilon$ 也不能过大。

若：

$$
\varepsilon
\gg
\mathbf u^H[k]\mathbf u[k]
$$

则更新近似为：

$$
\Delta\mathbf g[k]
\approx
\frac{\mu}{\varepsilon}
\mathbf u[k]e^*[k]
$$

归一化效果减弱，算法更接近一个小步长 LMS。

---

## 27.29 LMS 与 NLMS 的核心区别

LMS：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

NLMS：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\frac{\mu}
{\|\mathbf u[k]\|^2+\varepsilon}
\mathbf u[k]e^*[k]
$$

两者的更新方向相同：

$$
\mathbf u[k]e^*[k]
$$

区别只在于有效步长。

LMS 的有效步长固定：

$$
\mu_{\mathrm{eff}}[k]=\mu
$$

NLMS 的有效步长随输入能量变化：

$$
\mu_{\mathrm{eff}}[k]
=
\frac{\mu}
{\|\mathbf u[k]\|^2+\varepsilon}
$$

输入能量大时，NLMS 自动减小步长；输入能量小时，NLMS 自动增大相对步长。

---

## 27.30 为什么 NLMS 对输入缩放更不敏感

假设输入整体乘以常数 $c$：

$$
\mathbf u'[k]
=
c\mathbf u[k]
$$

NLMS 更新中的分母变为：

$$
\|\mathbf u'[k]\|^2
=
|c|^2\|\mathbf u[k]\|^2
$$

输入幅度放大时，分母按照功率比例放大，从而抑制更新量的增长。

因此相比固定步长 LMS，NLMS 对输入功率变化更稳健。

严格来说，误差信号也会随系统状态变化，所以不能简单认为所有尺度完全抵消，但归一化显著降低了输入幅度对稳定性的影响。

---

## 27.31 LMS 和 NLMS 最终逼近什么

在平稳环境、步长合适、建模条件满足的情况下，LMS 和 NLMS 都试图逼近：

$$
\mathbf g_{\mathrm{opt}}
=
\mathbf R_u^{-1}\mathbf p_{ud}
$$

即 GSC 辅助支路的维纳解。

区别在于：

- 维纳解：直接使用统计量的闭式最优解；
- 最速下降：使用精确统计梯度逐步逼近；
- LMS：使用瞬时随机梯度逐步逼近；
- NLMS：在 LMS 基础上按当前输入能量调整有效步长。

---

## 27.32 GSC 中 LMS 更新到底在学习什么

GSC 固定支路输出：

$$
d[k]
=
\mathbf w_q^H\mathbf x[k]
$$

阻塞支路输出：

$$
\mathbf u[k]
=
\mathbf B^H\mathbf x[k]
$$

LMS 更新：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

只有当：

$$
\mathbf u[k]
$$

与：

$$
e[k]
$$

相关时，平均更新才不为零。

算法不断调整 $\mathbf g[k]$，直到：

$$
E[\mathbf u[k]e^*[k]]
=
\mathbf0
$$

这就是维纳正交原理。

其含义是：

> 最终输出中，所有能够由阻塞支路输入线性预测的成分都已经被辅助滤波器消除。

在理想 GSC 中，这些可预测成分主要是干扰。

---

## 27.33 正交原理与 LMS 稳态

维纳最优条件：

$$
E[\mathbf u e^*]
=
\mathbf0
$$

LMS 的瞬时更新方向：

$$
\mathbf u[k]e^*[k]
$$

因此 LMS 本质上是在不断尝试使其长期平均趋于零。

当平均相关不为零时，权值会继续调整。

当：

$$
E[\mathbf u e^*]
=
\mathbf0
$$

时，平均更新停止，但瞬时更新仍可能非零。

所以稳态下：

$$
E[\Delta\mathbf g[k]]
\approx
\mathbf0
$$

但通常：

$$
\Delta\mathbf g[k]
\neq
\mathbf0
$$

这就是平均稳定和瞬时抖动之间的区别。

---

## 27.34 GSC 中的符号约定

本章使用：

$$
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
$$

所以 LMS 更新是：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

因为总权值定义为：

$$
\mathbf w[k]
=
\mathbf w_q-
\mathbf B\mathbf g[k]
$$

某些教材使用：

$$
e[k]
=
d[k]
+
\mathbf g^H[k]\mathbf u[k]
$$

或者：

$$
\mathbf w[k]
=
\mathbf w_q+
\mathbf B\mathbf g[k]
$$

这时更新公式中的符号可能发生变化。

因此看不同教材时，首先检查：

1. 最终输出是加还是减；
2. 总权值中是 $+\mathbf B\mathbf g$ 还是 $-\mathbf B\mathbf g$；
3. 误差信号如何定义。

不要只看更新公式表面的正负号。

---

## 27.35 目标泄漏问题

理想阻塞条件为：

$$
\mathbf B^H\mathbf a_0
=
\mathbf0
$$

真实目标导向矢量可能是：

$$
\mathbf a_{\mathrm{true}}
$$

如果存在模型误差：

$$
\mathbf B^H\mathbf a_{\mathrm{true}}
\neq
\mathbf0
$$

那么目标信号会进入辅助支路：

$$
\mathbf u[k]
=
\text{目标泄漏}
+
\text{干扰}
+
\text{噪声}
$$

LMS 和 NLMS 不知道哪部分是目标、哪部分是干扰。

它们只会降低：

$$
E[|e[k]|^2]
$$

因此，只要目标泄漏与固定支路目标相关，算法就可能把目标也估计并减去。

这称为：

- 目标自消除；
- signal cancellation；
- signal leakage problem。

---

## 27.36 为什么步长大时目标自消除可能更快

当存在目标泄漏时，辅助支路中出现与主支路目标高度相关的分量。

LMS 更新方向：

$$
\mathbf u[k]e^*[k]
$$

会快速捕获这种相关性。

步长越大：

- 学习目标泄漏的速度越快；
- 目标自消除可能更快发生；
- 稳态权值波动也更明显。

因此，GSC 中步长设计不能只考虑干扰收敛速度，还必须考虑模型误差和目标保护。

---

## 27.37 常见改进方法

针对目标泄漏和稳定性问题，可以采用：

1. 更准确的阵列校准；
2. 鲁棒阻塞矩阵；
3. 多点约束；
4. 导数约束；
5. 对角加载；
6. 权值范数限制；
7. 泄漏 LMS；
8. 可变步长 LMS；
9. 目标活动检测；
10. 对辅助支路进行子空间限制。

---

## 27.38 泄漏 LMS 的基本思想

普通 LMS：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

泄漏 LMS 在每次更新时对旧权值进行轻微衰减：

$$
\mathbf g[k+1]
=
(1-
\mu\gamma)
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

其中：

$$
\gamma>0
$$

为泄漏系数。

它等价于在代价函数中加入权值范数惩罚：

$$
J_{\mathrm{reg}}
=
E[|e[k]|^2]
+
\gamma\|\mathbf g\|^2
$$

这样可以限制辅助权值过度增长，提高病态输入或模型误差下的稳定性。

但它也会引入偏差，使结果不再严格等于原始维纳解。

---

## 27.39 可变步长思想

固定步长必须在收敛速度和稳态误差之间折中。

可变步长算法通常采用：

- 初期误差大时使用较大步长；
- 接近稳态后逐渐减小步长；
- 环境突变时再次增大步长。

其目标是同时获得：

- 较快收敛；
- 较小稳态失调；
- 较好的跟踪能力。

NLMS 已经根据输入能量调整有效步长，但它不一定根据误差状态自动调整基本参数 $\mu$。

---

## 27.40 一个标量例子

假设：

$$
d[k]
=
s[k]
+
2v[k]
$$

阻塞支路得到：

$$
u[k]=v[k]
$$

希望通过：

$$
g^*[k]u[k]
$$

估计固定支路中的：

$$
2v[k]
$$

若所有量暂时取实数，LMS 变为：

$$
g[k+1]
=
g[k]
+
\mu u[k]e[k]
$$

其中：

$$
e[k]
=
d[k]-g[k]u[k]
$$

如果忽略目标并令某时刻：

$$
u[k]=1
$$

$$
d[k]=2
$$

初始：

$$
g[0]=0
$$

则：

$$
e[0]=2
$$

更新：

$$
g[1]
=
0+2\mu
$$

若：

$$
\mu=0.2
$$

则：

$$
g[1]=0.4
$$

下一次若样本仍相同：

$$
e[1]
=
2-0.4
=
1.6
$$

更新：

$$
g[2]
=
0.4+0.2\times1\times1.6
=
0.72
$$

权值会逐渐接近：

$$
g_{\mathrm{opt}}=2
$$

这说明 LMS 通过误差反馈，逐渐学习固定支路中干扰的比例系数。

---

## 27.41 同一标量例子下的 NLMS

NLMS：

$$
g[k+1]
=
g[k]
+
\frac{\mu}
{|u[k]|^2+\varepsilon}
u[k]e^*[k]
$$

忽略 $\varepsilon$，若：

$$
u[k]=1
$$

则 NLMS 与步长为 $\mu$ 的 LMS 相同。

但若输入放大到：

$$
u[k]=10
$$

普通 LMS 的更新量包含：

$$
10e[k]
$$

可能非常大。

NLMS 分母为：

$$
|u[k]|^2=100
$$

所以更新比例变为：

$$
\frac{10e[k]}{100}
=
0.1e[k]
$$

输入幅度增大后，归一化会自动减小有效更新量。

---

## 27.42 LMS、NLMS 与 RLS 的区别

### LMS

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

特点：

- 计算量低；
- 实现简单；
- 对输入功率敏感；
- 特征值扩展大时收敛慢。

### NLMS

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\frac{\mu}
{\|\mathbf u[k]\|^2+\varepsilon}
\mathbf u[k]e^*[k]
$$

特点：

- 对输入功率变化更稳健；
- 参数更容易选择；
- 仍会受输入相关性影响；
- 计算量只比 LMS 略高。

### RLS

RLS 通过递推方式近似求解加权最小二乘问题。

特点：

- 收敛快；
- 对特征值扩展不太敏感；
- 计算量明显更高；
- 对数值稳定性要求更高。

---

## 27.43 常见误区

### 误区一：LMS 每一步都在减小误差

不一定。

LMS 使用随机瞬时梯度，单次更新可能使下一时刻误差变大。

它追求的是平均意义上的下降趋势。

### 误区二：到达维纳解后更新会完全停止

不一定。

维纳解处满足：

$$
E[\mathbf u e^*]=\mathbf0
$$

但瞬时：

$$
\mathbf u[k]e^*[k]
$$

通常不为零，所以权值仍会抖动。

### 误区三：步长越大越好

步长大虽然初期快，但可能造成：

- 震荡；
- 发散；
- 稳态失调增大；
- 目标自消除加快。

### 误区四：NLMS 完全解决相关输入问题

NLMS 主要解决输入能量尺度变化问题。

它不能彻底消除协方差矩阵特征值扩展带来的慢收敛。

### 误区五：阻塞支路中一定没有目标

只有在导向矢量完全准确、阵列完全校准且模型成立时才近似成立。

实际系统必须考虑目标泄漏。

---

## 27.44 公式之间的完整关系

维纳解：

$$
\mathbf g_{\mathrm{opt}}
=
\mathbf R_u^{-1}\mathbf p_{ud}
$$

维纳－霍普方程：

$$
\mathbf R_u\mathbf g_{\mathrm{opt}}
=
\mathbf p_{ud}
$$

真实梯度：

$$
\nabla_{\mathbf g^*}J
=
\mathbf R_u\mathbf g
-
\mathbf p_{ud}
$$

最速下降：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu
\left(
\mathbf p_{ud}
-
\mathbf R_u\mathbf g[k]
\right)
$$

瞬时负梯度：

$$
\mathbf u[k]e^*[k]
$$

LMS：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

NLMS：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\frac{\mu}
{\mathbf u^H[k]\mathbf u[k]+\varepsilon}
\mathbf u[k]e^*[k]
$$

GSC 总权值：

$$
\mathbf w[k]
=
\mathbf w_q
-
\mathbf B\mathbf g[k]
$$

GSC 输出：

$$
e[k]
=
\mathbf w^H[k]\mathbf x[k]
$$

同时也可写成：

$$
e[k]
=
d[k]
-
\mathbf g^H[k]\mathbf u[k]
$$

---

## 27.45 本章核心理解

LMS 不是凭空出现的更新公式，它的完整来源是：

```text
维纳最优解
        ↓
均方误差代价函数
        ↓
对复权值求 Wirtinger 梯度
        ↓
最速下降法
        ↓
用当前样本替代未知统计期望
        ↓
LMS
```

即：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\mu\mathbf u[k]e^*[k]
$$

NLMS 则进一步解决输入能量变化导致的步长敏感问题：

```text
LMS 更新量受输入能量影响
        ↓
根据当前输入能量归一化
        ↓
NLMS
```

即：

$$
\mathbf g[k+1]
=
\mathbf g[k]
+
\frac{\mu}
{\|\mathbf u[k]\|^2+\varepsilon}
\mathbf u[k]e^*[k]
$$

在 GSC 中，LMS 和 NLMS 的本质任务是：

> 在不破坏约束的零空间内，逐步学习一个干扰估计器，使最终输出与阻塞支路输入在平均意义下正交。

理想稳态条件为：

$$
E[\mathbf u[k]e^*[k]]
=
\mathbf0
$$

这表示所有能够由辅助支路线性预测的输出成分，都已经被消除。

但实际系统中必须同时关注：

- 步长稳定性；
- 输入相关性；
- 稳态失调；
- 环境跟踪速度；
- 阻塞矩阵误差；
- 目标泄漏与目标自消除。

这些因素共同决定 LMS 和 NLMS 在真实阵列抗干扰系统中的效果。
