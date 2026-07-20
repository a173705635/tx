# 拉格朗日乘子法：从直观理解到 MVDR、LCMV、GSC

你在波束形成里经常看到这种问题：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

同时要求：

$$
\mathbf w^H\mathbf a=1
$$

或者：

$$
\mathbf C^H\mathbf w=\mathbf f
$$

这类问题不能像普通函数一样，只令梯度等于零，因为权向量还必须满足约束。

拉格朗日乘子法的核心作用是：

> 把“目标函数最小”和“必须满足约束”这两件事，合并到同一个函数中处理。

---

## 1. 什么是优化问题

优化问题一般由两部分组成。

### 1.1 目标函数

也就是希望最小化或最大化的量。

例如：

$$
\min_{\mathbf x} f(\mathbf x)
$$

在波束形成中经常是输出功率：

$$
f(\mathbf w)
=
\mathbf w^H\mathbf R\mathbf w
$$

### 1.2 约束条件

变量不能任意选择，必须满足某些条件。

例如：

$$
g(\mathbf x)=0
$$

在 MVDR 中：

$$
\mathbf w^H\mathbf a-1=0
$$

表示目标方向响应必须保持为 1。

所以 MVDR 是：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf w^H\mathbf a=1
$$

---

## 2. 无约束优化为什么令梯度为零

先考虑没有约束的情况：

$$
\min_{\mathbf x} f(\mathbf x)
$$

在一个光滑函数的最低点，向任何方向稍微移动，都不能继续下降。

因此该点的梯度为零：

$$
\nabla f(\mathbf x)=\mathbf0
$$

一维情况下就是：

$$
\frac{df(x)}{dx}=0
$$

例如：

$$
f(x)=x^2
$$

求导：

$$
f'(x)=2x
$$

令其为零：

$$
2x=0
$$

得到：

$$
x=0
$$

---

## 3. 有约束时，为什么不能直接令梯度为零

考虑：

$$
\min_{x,y}
x^2+y^2
$$

约束为：

$$
x+y=1
$$

如果没有约束，最小值点显然是：

$$
x=0,\qquad y=0
$$

但这个点不满足：

$$
x+y=1
$$

所以不能选。

现在变量只能在线：

$$
x+y=1
$$

上移动。

最优点是这条直线上距离原点最近的点：

$$
x=y=\frac12
$$

在这个点：

$$
\nabla f
=
\begin{bmatrix}
2x\\
2y
\end{bmatrix}
=
\begin{bmatrix}
1\\
1
\end{bmatrix}
$$

梯度并不等于零。

原因是：

> 虽然沿梯度反方向可以继续减小目标函数，但那个方向会离开约束直线，所以不允许走。

---

## 4. 约束曲线的法向量

约束写成：

$$
g(x,y)=x+y-1=0
$$

约束函数的梯度为：

$$
\nabla g
=
\begin{bmatrix}
1\\
1
\end{bmatrix}
$$

梯度：

$$
\nabla g
$$

垂直于约束曲线或约束平面。

因此它是约束面的法向量。

在最优点，目标函数梯度：

$$
\nabla f
$$

也必须垂直于约束面。

否则，如果 $\nabla f$ 还有沿约束面的分量，就可以沿着约束面继续降低目标函数。

所以最优点满足：

$$
\nabla f
=
\lambda\nabla g
$$

或者写成：

$$
\nabla f-\lambda\nabla g=\mathbf0
$$

这里的：

$$
\lambda
$$

就是拉格朗日乘子。

---

## 5. 拉格朗日函数从哪里来

对于问题：

$$
\min_{\mathbf x}f(\mathbf x)
$$

满足：

$$
g(\mathbf x)=0
$$

构造拉格朗日函数：

$$
L(\mathbf x,\lambda)
=
f(\mathbf x)
+
\lambda g(\mathbf x)
$$

也可以写成：

$$
L(\mathbf x,\lambda)
=
f(\mathbf x)
-
\lambda g(\mathbf x)
$$

正负号都可以，只会使最终的拉格朗日乘子符号相反，不会改变最优变量。

对 $\mathbf x$ 求梯度：

$$
\nabla_{\mathbf x}L
=
\nabla f
+
\lambda\nabla g
$$

令其为零：

$$
\nabla f
+
\lambda\nabla g
=
\mathbf0
$$

这表示：

$$
\nabla f
=
-\lambda\nabla g
$$

也就是目标函数梯度与约束面法向量平行。

对 $\lambda$ 求导：

$$
\frac{\partial L}{\partial\lambda}
=
g(\mathbf x)
$$

令其为零：

$$
g(\mathbf x)=0
$$

这又恢复了原来的约束条件。

因此拉格朗日函数同时产生：

1. 最优性条件；
2. 约束条件。

---

## 6. 第一个完整例子

求解：

$$
\min_{x,y}
x^2+y^2
$$

满足：

$$
x+y=1
$$

构造：

$$
L(x,y,\lambda)
=
x^2+y^2
+
\lambda(x+y-1)
$$

分别求偏导。

对 $x$：

$$
\frac{\partial L}{\partial x}
=
2x+\lambda
=
0
$$

对 $y$：

$$
\frac{\partial L}{\partial y}
=
2y+\lambda
=
0
$$

对 $\lambda$：

$$
\frac{\partial L}{\partial\lambda}
=
x+y-1
=
0
$$

前两式相减：

$$
2x-2y=0
$$

所以：

$$
x=y
$$

代入约束：

$$
x+y=1
$$

得到：

$$
x=y=\frac12
$$

---

## 7. 拉格朗日乘子本身是什么意思

从几何上看，$\lambda$ 表示：

> 约束面对目标函数梯度施加了多大的“平衡作用”。

无约束情况下，最优点要求：

$$
\nabla f=\mathbf0
$$

有约束时，目标函数梯度可以不为零，但它必须被约束法向量抵消：

$$
\nabla f+\lambda\nabla g=\mathbf0
$$

可以类比成：

```text
目标函数梯度：
想把解推向更低的位置

约束条件：
不允许解离开约束面

拉格朗日乘子：
表示约束需要施加多大的限制作用
```

在经济学中，拉格朗日乘子还经常解释为“影子价格”：

> 把约束放宽一点，最优目标值会变化多少。

---

## 8. 多个等式约束

如果有多个约束：

$$
g_1(\mathbf x)=0
$$

$$
g_2(\mathbf x)=0
$$

一直到：

$$
g_Q(\mathbf x)=0
$$

分别引入：

$$
\lambda_1,\lambda_2,\ldots,\lambda_Q
$$

拉格朗日函数为：

$$
L
=
f(\mathbf x)
+
\sum_{q=1}^{Q}
\lambda_qg_q(\mathbf x)
$$

矩阵形式可以写成：

$$
L
=
f(\mathbf x)
+
\boldsymbol{\lambda}^T
\mathbf g(\mathbf x)
$$

最优点满足：

$$
\nabla f
+
\sum_{q=1}^{Q}
\lambda_q\nabla g_q
=
\mathbf0
$$

也就是说：

> 目标函数梯度必须位于所有约束法向量张成的空间内。

---

## 9. 二次型优化问题

波束形成中最常见的是二次型：

$$
f(\mathbf x)
=
\mathbf x^T\mathbf Q\mathbf x
$$

假设：

$$
\mathbf Q=\mathbf Q^T
$$

则梯度为：

$$
\nabla_{\mathbf x}
\left(
\mathbf x^T\mathbf Q\mathbf x
\right)
=
2\mathbf Q\mathbf x
$$

如果目标函数写成：

$$
f(\mathbf x)
=
\frac12
\mathbf x^T\mathbf Q\mathbf x
$$

那么：

$$
\nabla f
=
\mathbf Q\mathbf x
$$

所以数学推导中经常在二次型前乘：

$$
\frac12
$$

只是为了消掉求导后出现的 2。

最优解本身不会因此改变。

---

## 10. 一个实数矩阵约束问题

考虑：

$$
\min_{\mathbf x}
\mathbf x^T\mathbf Q\mathbf x
$$

满足：

$$
\mathbf A^T\mathbf x=\mathbf b
$$

构造拉格朗日函数：

$$
L(\mathbf x,\boldsymbol{\lambda})
=
\mathbf x^T\mathbf Q\mathbf x
-
2\boldsymbol{\lambda}^T
\left(
\mathbf A^T\mathbf x-\mathbf b
\right)
$$

这里乘 2 只是为了简化后面的公式。

对 $\mathbf x$ 求梯度：

$$
\nabla_{\mathbf x}L
=
2\mathbf Q\mathbf x
-
2\mathbf A\boldsymbol{\lambda}
$$

令其为零：

$$
\mathbf Q\mathbf x
=
\mathbf A\boldsymbol{\lambda}
$$

如果 $\mathbf Q$ 可逆：

$$
\mathbf x
=
\mathbf Q^{-1}
\mathbf A
\boldsymbol{\lambda}
$$

代入约束：

$$
\mathbf A^T
\mathbf Q^{-1}
\mathbf A
\boldsymbol{\lambda}
=
\mathbf b
$$

所以：

$$
\boldsymbol{\lambda}
=
\left(
\mathbf A^T
\mathbf Q^{-1}
\mathbf A
\right)^{-1}
\mathbf b
$$

最终：

$$
\mathbf x_{\mathrm{opt}}
=
\mathbf Q^{-1}
\mathbf A
\left(
\mathbf A^T
\mathbf Q^{-1}
\mathbf A
\right)^{-1}
\mathbf b
$$

这已经和 LCMV 的形式完全一样。

---

## 11. 为什么这个解看起来像投影

如果：

$$
\mathbf Q=\mathbf I
$$

那么：

$$
\mathbf x_{\mathrm{opt}}
=
\mathbf A
\left(
\mathbf A^T\mathbf A
\right)^{-1}
\mathbf b
$$

这相当于在满足约束的所有向量中，寻找欧氏范数最小的向量。

如果：

$$
\mathbf Q\neq\mathbf I
$$

则目标为：

$$
\mathbf x^T\mathbf Q\mathbf x
$$

这不是普通欧氏距离，而是由 $\mathbf Q$ 定义的加权距离。

所以可以理解为：

> 拉格朗日解是在约束平面上，寻找离原点最近的点，只不过“距离”的定义由 $\mathbf Q$ 决定。

---

## 12. 一个加权二维例子

求：

$$
\min_{x_1,x_2}
x_1^2+4x_2^2
$$

满足：

$$
x_1+x_2=1
$$

这里：

$$
\mathbf Q
=
\begin{bmatrix}
1&0\\
0&4
\end{bmatrix}
$$

构造：

$$
L
=
x_1^2+4x_2^2
-
2\lambda(x_1+x_2-1)
$$

求导：

$$
2x_1-2\lambda=0
$$

所以：

$$
x_1=\lambda
$$

另一式：

$$
8x_2-2\lambda=0
$$

所以：

$$
x_2=\frac{\lambda}{4}
$$

代入约束：

$$
\lambda+\frac{\lambda}{4}=1
$$

得到：

$$
\lambda=\frac45
$$

所以：

$$
x_1=\frac45
$$

$$
x_2=\frac15
$$

为什么不是两个都等于 $1/2$？

因为目标函数认为 $x_2$ 比较“昂贵”：

$$
4x_2^2
$$

所以最优解尽量减少 $x_2$，更多使用 $x_1$。

这和 MVDR 非常相似：

> 协方差矩阵中功率较高的方向代价更大，最优权值会尽量避开这些方向。

---

## 13. 复数变量为什么更麻烦

阵列权值通常为复数：

$$
\mathbf w\in\mathbb C^{N\times1}
$$

目标函数为：

$$
\mathbf w^H\mathbf R\mathbf w
$$

这里：

$$
\mathbf w^H
$$

包含 $\mathbf w$ 的共轭。

因此不能完全按照普通实数微积分处理。

复数优化中常用 Wirtinger 微分思想：

> 暂时把 $\mathbf w$ 和 $\mathbf w^*$ 看成两个独立变量。

对于实值目标函数，最优点可以通过：

$$
\frac{\partial L}
{\partial\mathbf w^*}
=
\mathbf0
$$

求得。

---

## 14. 常用复数矩阵求导公式

如果：

$$
\mathbf R=\mathbf R^H
$$

则：

$$
\frac{\partial}
{\partial\mathbf w^*}
\left(
\mathbf w^H\mathbf R\mathbf w
\right)
=
\mathbf R\mathbf w
$$

另外：

$$
\frac{\partial}
{\partial\mathbf w^*}
\left(
\mathbf w^H\mathbf a
\right)
=
\mathbf a
$$

而：

$$
\frac{\partial}
{\partial\mathbf w^*}
\left(
\mathbf a^H\mathbf w
\right)
=
\mathbf0
$$

因为在 Wirtinger 微分中，$\mathbf w$ 和 $\mathbf w^*$ 被暂时视为独立变量。

---

## 15. MVDR 的拉格朗日推导

MVDR 问题为：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf w^H\mathbf a_0=1
$$

因为这是复数约束，还要同时考虑它的共轭：

$$
\mathbf a_0^H\mathbf w=1
$$

构造实值拉格朗日函数：

$$
L
=
\mathbf w^H\mathbf R\mathbf w
-
\lambda^*
\left(
\mathbf w^H\mathbf a_0-1
\right)
-
\lambda
\left(
\mathbf a_0^H\mathbf w-1
\right)
$$

对 $\mathbf w^*$ 求偏导：

$$
\frac{\partial L}
{\partial\mathbf w^*}
=
\mathbf R\mathbf w
-
\lambda^*\mathbf a_0
$$

令其为零：

$$
\mathbf R\mathbf w
=
\lambda^*\mathbf a_0
$$

所以：

$$
\mathbf w
=
\lambda^*
\mathbf R^{-1}\mathbf a_0
$$

现在使用约束：

$$
\mathbf w^H\mathbf a_0=1
$$

由：

$$
\mathbf w^H
=
\lambda
\mathbf a_0^H\mathbf R^{-1}
$$

得到：

$$
\lambda
\mathbf a_0^H
\mathbf R^{-1}
\mathbf a_0
=
1
$$

因此：

$$
\lambda
=
\frac{
1
}{
\mathbf a_0^H
\mathbf R^{-1}
\mathbf a_0
}
$$

因为分母为实数正值，所以最终：

$$
\mathbf w_{\mathrm{MVDR}}
=
\frac{
\mathbf R^{-1}\mathbf a_0
}{
\mathbf a_0^H
\mathbf R^{-1}
\mathbf a_0
}
$$

---

## 16. MVDR 推导的每一步是什么意思

从：

$$
\mathbf R\mathbf w
=
\lambda^*\mathbf a_0
$$

可以看出：

$$
\mathbf w
\propto
\mathbf R^{-1}\mathbf a_0
$$

也就是先使用：

$$
\mathbf R^{-1}
$$

对目标导向矢量进行加权。

它会压低高功率干扰方向，放大低功率方向。

但：

$$
\mathbf R^{-1}\mathbf a_0
$$

还不一定满足目标响应为 1。

所以再除以：

$$
\mathbf a_0^H
\mathbf R^{-1}
\mathbf a_0
$$

完成归一化。

因此 MVDR 公式可以理解成：

```text
R^{-1}a0：
根据干扰环境调整方向

分母：
重新归一化，使目标方向响应等于 1
```

---

## 17. CBF 也可以由拉格朗日法得到

常规波束形成可以看成：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf w
$$

满足：

$$
\mathbf w^H\mathbf a_0=1
$$

这表示：

> 在满足目标方向无失真的所有权值中，选择权值范数最小的一个。

它相当于 MVDR 中：

$$
\mathbf R=\mathbf I
$$

代入 MVDR 公式：

$$
\mathbf w_{\mathrm{CBF}}
=
\frac{
\mathbf a_0
}{
\mathbf a_0^H\mathbf a_0
}
$$

所以 CBF 和 MVDR 的区别可以理解为：

- CBF 使用普通欧氏范数作为代价；
- MVDR 使用协方差矩阵定义的加权功率作为代价。

---

## 18. LCMV 的拉格朗日推导

LCMV 问题为：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf C^H\mathbf w
=
\mathbf f
$$

其中：

$$
\mathbf C\in\mathbb C^{N\times Q}
$$

$$
\mathbf f\in\mathbb C^{Q\times1}
$$

引入复数拉格朗日乘子：

$$
\boldsymbol{\lambda}
\in\mathbb C^{Q\times1}
$$

构造：

$$
L
=
\mathbf w^H\mathbf R\mathbf w
-
\boldsymbol{\lambda}^H
\left(
\mathbf C^H\mathbf w-\mathbf f
\right)
-
\left(
\mathbf w^H\mathbf C-\mathbf f^H
\right)
\boldsymbol{\lambda}
$$

对 $\mathbf w^*$ 求导：

$$
\frac{\partial L}
{\partial\mathbf w^*}
=
\mathbf R\mathbf w
-
\mathbf C\boldsymbol{\lambda}
$$

令其为零：

$$
\mathbf R\mathbf w
=
\mathbf C\boldsymbol{\lambda}
$$

所以：

$$
\mathbf w
=
\mathbf R^{-1}
\mathbf C
\boldsymbol{\lambda}
$$

代入约束：

$$
\mathbf C^H
\mathbf R^{-1}
\mathbf C
\boldsymbol{\lambda}
=
\mathbf f
$$

得到：

$$
\boldsymbol{\lambda}
=
\left(
\mathbf C^H
\mathbf R^{-1}
\mathbf C
\right)^{-1}
\mathbf f
$$

最终：

$$
\mathbf w_{\mathrm{LCMV}}
=
\mathbf R^{-1}
\mathbf C
\left(
\mathbf C^H
\mathbf R^{-1}
\mathbf C
\right)^{-1}
\mathbf f
$$

---

## 19. MVDR 是 LCMV 的特殊情况

如果只有一个约束：

$$
\mathbf C=\mathbf a_0
$$

以及：

$$
\mathbf f=1
$$

那么 LCMV 公式变成：

$$
\mathbf w
=
\mathbf R^{-1}\mathbf a_0
\left(
\mathbf a_0^H
\mathbf R^{-1}
\mathbf a_0
\right)^{-1}
$$

也就是 MVDR：

$$
\mathbf w_{\mathrm{MVDR}}
=
\frac{
\mathbf R^{-1}\mathbf a_0
}{
\mathbf a_0^H
\mathbf R^{-1}
\mathbf a_0
}
$$

---

## 20. 功率反演也是同一类问题

功率反演可以写成：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf e_1^H\mathbf w=1
$$

其中：

$$
\mathbf e_1
=
\begin{bmatrix}
1\\
0\\
\vdots\\
0
\end{bmatrix}
$$

所以它的解是：

$$
\mathbf w_{\mathrm{PI}}
=
\frac{
\mathbf R^{-1}\mathbf e_1
}{
\mathbf e_1^H
\mathbf R^{-1}
\mathbf e_1
}
$$

因此：

- MVDR 约束目标方向；
- PI 约束参考阵元；
- 两者的数学结构完全相同。

---

## 21. 拉格朗日法和 GSC 的关系

LCMV 直接求解：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf C^H\mathbf w=\mathbf f
$$

GSC 则把所有满足约束的权向量参数化为：

$$
\mathbf w
=
\mathbf w_q-\mathbf B\mathbf g
$$

其中：

$$
\mathbf C^H\mathbf w_q=\mathbf f
$$

并且：

$$
\mathbf C^H\mathbf B=\mathbf0
$$

因此：

$$
\mathbf C^H\mathbf w
=
\mathbf C^H\mathbf w_q
-
\mathbf C^H\mathbf B\mathbf g
$$

得到：

$$
\mathbf C^H\mathbf w
=
\mathbf f
$$

所以无论 $\mathbf g$ 怎么变化，约束始终满足。

---

## 22. GSC 的几何意义

约束：

$$
\mathbf C^H\mathbf w=\mathbf f
$$

定义了一个可行集合。

$\mathbf w_q$ 是可行集合中的一个点。

$\mathbf B$ 的列向量位于：

$$
\operatorname{null}(\mathbf C^H)
$$

所以：

$$
\mathbf B\mathbf g
$$

表示沿着约束平面内部移动。

因此 GSC 可以理解为：

```text
w_q：
先找到一个满足约束的点

B：
找到所有不会破坏约束的移动方向

g：
在这些允许的方向上优化输出功率
```

拉格朗日法是直接求最终最优点。

GSC 是先描述整个可行空间，再在可行空间内部自适应搜索。

二者理论上得到同一个最优解。

---

## 23. 为什么导数等于零只是必要条件

构造拉格朗日函数并令导数为零，只能得到候选点。

还需要判断它是不是最小值。

例如：

$$
f(x)=-x^2
$$

在：

$$
x=0
$$

处也有：

$$
f'(0)=0
$$

但这是最大值，不是最小值。

对于二次型：

$$
\mathbf x^T\mathbf Q\mathbf x
$$

如果：

$$
\mathbf Q\succ0
$$

也就是 $\mathbf Q$ 正定，那么目标函数是严格凸函数。

此时拉格朗日条件得到的可行点就是唯一全局最小值。

在波束形成中，协方差矩阵通常满足：

$$
\mathbf R\succeq0
$$

如果：

$$
\mathbf R\succ0
$$

则 MVDR 或 LCMV 解通常唯一。

---

## 24. 什么是正定

矩阵 $\mathbf R$ 正定是指，对于任意非零向量：

$$
\mathbf z^H\mathbf R\mathbf z>0
$$

协方差矩阵的二次型：

$$
\mathbf w^H\mathbf R\mathbf w
$$

表示输出功率，因此不可能为负。

所以协方差矩阵至少是半正定的：

$$
\mathbf R\succeq0
$$

如果所有方向都有非零功率，则：

$$
\mathbf R\succ0
$$

此时矩阵可逆。

如果快拍数不足或通道线性相关，可能出现：

$$
\mathbf R
$$

秩亏，不可逆。

这时通常采用：

$$
\mathbf R+\delta\mathbf I
$$

进行对角加载。

---

## 25. 不等式约束与 KKT 条件

拉格朗日乘子法还可以处理不等式约束。

一般问题：

$$
\min_{\mathbf x}f(\mathbf x)
$$

满足：

$$
g_i(\mathbf x)=0
$$

以及：

$$
h_j(\mathbf x)\le0
$$

对等式约束使用乘子：

$$
\lambda_i
$$

对不等式约束使用乘子：

$$
\mu_j
$$

拉格朗日函数为：

$$
L
=
f(\mathbf x)
+
\sum_i\lambda_i g_i(\mathbf x)
+
\sum_j\mu_j h_j(\mathbf x)
$$

KKT 条件包括四部分。

### 25.1 驻点条件

$$
\nabla_{\mathbf x}L
=
\mathbf0
$$

### 25.2 原始可行性

$$
g_i(\mathbf x)=0
$$

$$
h_j(\mathbf x)\le0
$$

### 25.3 对偶可行性

$$
\mu_j\ge0
$$

### 25.4 互补松弛

$$
\mu_jh_j(\mathbf x)=0
$$

---

## 26. 互补松弛怎样理解

条件：

$$
\mu_jh_j(\mathbf x)=0
$$

意味着两种可能。

### 约束没有卡住最优解

如果：

$$
h_j(\mathbf x)<0
$$

说明约束还有余量，没有真正起作用。

那么必须：

$$
\mu_j=0
$$

### 约束正好卡住最优解

如果：

$$
h_j(\mathbf x)=0
$$

说明约束处于边界上。

此时：

$$
\mu_j
$$

可以大于零。

所以：

> 只有真正处于激活状态的约束，才会对最优解施加作用。

---

## 27. 一个不等式约束例子

求：

$$
\min_x x^2
$$

满足：

$$
x\ge1
$$

写成标准形式：

$$
1-x\le0
$$

拉格朗日函数：

$$
L(x,\mu)
=
x^2+\mu(1-x)
$$

驻点条件：

$$
2x-\mu=0
$$

原始可行：

$$
x\ge1
$$

对偶可行：

$$
\mu\ge0
$$

互补松弛：

$$
\mu(1-x)=0
$$

最优点显然是：

$$
x=1
$$

所以：

$$
\mu=2
$$

这里约束处于激活状态。

---

## 28. 为什么拉格朗日函数可以把约束合进去

拉格朗日法不是简单地把约束当成“惩罚项”。

普通惩罚法可能写成：

$$
f(\mathbf x)
+
\rho g^2(\mathbf x)
$$

它只是让违反约束变得昂贵，但有限的 $\rho$ 下不一定严格满足约束。

拉格朗日函数则通过：

$$
\lambda g(\mathbf x)
$$

引入一个会自动调整的乘子。

求解时同时要求：

$$
\frac{\partial L}{\partial\lambda}
=
g(\mathbf x)
=
0
$$

因此等式约束会被严格满足。

---

## 29. 拉格朗日对偶的基本思想

对固定的乘子 $\lambda$，先计算：

$$
q(\lambda)
=
\inf_{\mathbf x}
L(\mathbf x,\lambda)
$$

这个函数称为对偶函数。

然后考虑：

$$
\max_{\lambda}q(\lambda)
$$

称为对偶问题。

直观上：

- 原问题直接寻找满足约束的最优变量；
- 对偶问题寻找最合适的“约束价格”；
- 在凸优化和适当条件下，二者最优值相同。

波束形成基础推导通常不需要完整求解对偶问题，但拉格朗日乘子法正是对偶理论的入口。

---

## 30. 常见错误一：漏掉约束条件

只求：

$$
\nabla_{\mathbf x}L=0
$$

是不够的。

还必须同时使用：

$$
g(\mathbf x)=0
$$

否则无法确定拉格朗日乘子，也无法保证结果满足约束。

---

## 31. 常见错误二：纠结拉格朗日函数的正负号

有人写：

$$
L=f+\lambda g
$$

有人写：

$$
L=f-\lambda g
$$

两种都正确。

区别只在于：

$$
\lambda
$$

的符号不同。

最终的最优变量相同。

---

## 32. 常见错误三：复数约束只写一项

对复数约束：

$$
\mathbf w^H\mathbf a=1
$$

为了使拉格朗日函数保持实值，通常同时写约束及其共轭：

$$
-\lambda^*
(\mathbf w^H\mathbf a-1)
$$

以及：

$$
-\lambda
(\mathbf a^H\mathbf w-1)
$$

这样才能规范地对复变量求导。

一些教材会使用更简略的复数拉格朗日写法，最终结论相同，但中间过程可能看起来不同。

---

## 33. 常见错误四：把转置和共轭转置混淆

实数情况使用：

$$
\mathbf x^T\mathbf Q\mathbf x
$$

复数情况使用：

$$
\mathbf w^H\mathbf R\mathbf w
$$

不能写成：

$$
\mathbf w^T\mathbf R\mathbf w
$$

因为功率和复数内积必须使用共轭转置。

---

## 34. 常见错误五：直接计算逆矩阵

公式中经常写：

$$
\mathbf R^{-1}\mathbf a
$$

但程序中不建议显式求逆：

```python
w_temp = np.linalg.inv(R) @ a
```

更建议解线性方程：

```python
w_temp = np.linalg.solve(R, a)
```

因为后者通常：

- 数值更稳定；
- 计算量更小；
- 误差更低。

---

## 35. 常见错误六：矩阵不可逆

LCMV 中需要：

$$
\mathbf R^{-1}
$$

以及：

$$
\left(
\mathbf C^H
\mathbf R^{-1}
\mathbf C
\right)^{-1}
$$

如果：

- 快拍不足；
- 约束重复；
- 约束向量线性相关；
- 协方差矩阵秩亏；

这些矩阵可能不可逆。

处理方法包括：

- 对角加载；
- 删除重复约束；
- 使用伪逆；
- 增加快拍；
- 改善阵列校准。

---

## 36. 为什么 LCMV 中间矩阵是 $Q\times Q$

假设：

$$
\mathbf C\in\mathbb C^{N\times Q}
$$

那么：

$$
\mathbf R^{-1}\mathbf C
\in
\mathbb C^{N\times Q}
$$

所以：

$$
\mathbf C^H
\mathbf R^{-1}
\mathbf C
$$

维度为：

$$
(Q\times N)
(N\times N)
(N\times Q)
=
Q\times Q
$$

它描述不同约束之间在：

$$
\mathbf R^{-1}
$$

度量下的相互关系。

如果约束独立，它通常可逆。

---

## 37. 为什么乘子数量等于约束数量

一个等式约束对应一个独立法向方向。

如果有 $Q$ 个约束，就需要 $Q$ 个乘子：

$$
\boldsymbol{\lambda}
=
\begin{bmatrix}
\lambda_1\\
\lambda_2\\
\vdots\\
\lambda_Q
\end{bmatrix}
$$

每个乘子表示对应约束对最优解施加的作用强度。

这也是为什么 LCMV 中：

$$
\boldsymbol{\lambda}\in\mathbb C^{Q\times1}
$$

---

## 38. 从几何上重新理解 MVDR

MVDR 可行集合是：

$$
\left\{
\mathbf w
\mid
\mathbf w^H\mathbf a_0=1
\right\}
$$

这是权值空间中的一个超平面。

输出功率：

$$
\mathbf w^H\mathbf R\mathbf w
$$

对应一系列椭球等值面。

优化过程就是：

> 不断扩大最小的功率椭球，直到它第一次接触约束超平面。

接触点就是 MVDR 最优权值。

在接触点：

- 椭球法向量为 $\mathbf R\mathbf w$；
- 约束面法向量为 $\mathbf a_0$；
- 两者平行。

所以：

$$
\mathbf R\mathbf w
=
\lambda\mathbf a_0
$$

这正是拉格朗日驻点条件。

---

## 39. 为什么求解结果中一定出现 $\mathbf R^{-1}$

驻点条件为：

$$
\mathbf R\mathbf w
=
\mathbf C\boldsymbol{\lambda}
$$

为了求出 $\mathbf w$，需要消除左侧的 $\mathbf R$：

$$
\mathbf w
=
\mathbf R^{-1}
\mathbf C
\boldsymbol{\lambda}
$$

所以 $\mathbf R^{-1}$ 并不是凭空出现的。

它来自：

> 目标函数二次型的梯度中含有 $\mathbf R\mathbf w$，解驻点方程时自然需要乘 $\mathbf R^{-1}$。

---

## 40. 各算法的统一形式

### CBF

$$
\min_{\mathbf w}
\mathbf w^H\mathbf I\mathbf w
$$

满足：

$$
\mathbf w^H\mathbf a_0=1
$$

### MVDR

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf w^H\mathbf a_0=1
$$

### LCMV

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf C^H\mathbf w=\mathbf f
$$

### 功率反演

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf e_1^H\mathbf w=1
$$

这些算法的区别，主要就是：

1. 目标函数中使用什么代价矩阵；
2. 约束矩阵和期望响应是什么。

---

## 41. 一套固定的拉格朗日推导流程

以后看到约束优化，可以按以下顺序操作。

### 第一步：写清楚目标函数

例如：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

### 第二步：把约束移到等号一边

例如：

$$
\mathbf C^H\mathbf w-\mathbf f=\mathbf0
$$

### 第三步：为每个约束引入乘子

$$
\boldsymbol{\lambda}
$$

### 第四步：构造拉格朗日函数

$$
L
=
\text{目标函数}
-
\text{乘子与约束的内积}
$$

### 第五步：对优化变量求导并令零

$$
\frac{\partial L}
{\partial\mathbf w^*}
=
\mathbf0
$$

### 第六步：用约束求乘子

将 $\mathbf w$ 的表达式代回：

$$
\mathbf C^H\mathbf w=\mathbf f
$$

### 第七步：把乘子代回得到最终解

### 第八步：检查

检查：

- 维度是否正确；
- 约束是否满足；
- 矩阵是否可逆；
- 目标函数是否凸；
- 是否真的是最小值。

---

## 42. 最核心的理解

拉格朗日乘子法并不是一个需要死记的技巧。

它解决的根本问题是：

> 当变量不能自由移动，只能在约束面上移动时，怎样找到目标函数在这个约束面上的最低点。

无约束最低点满足：

$$
\nabla f=\mathbf0
$$

有约束最低点通常不满足梯度为零，而是满足：

$$
\nabla f
+
\sum_i\lambda_i\nabla g_i
=
\mathbf0
$$

也就是：

> 目标函数梯度被约束法向量的线性组合抵消。

在波束形成中：

- 目标函数梯度来自输出功率；
- 约束法向量来自导向矢量或约束矩阵；
- 拉格朗日乘子负责协调二者；
- 最终得到既满足目标保护，又使输出功率最小的权向量。

最值得记住的统一公式是：

$$
\min_{\mathbf w}
\mathbf w^H\mathbf R\mathbf w
$$

满足：

$$
\mathbf C^H\mathbf w=\mathbf f
$$

对应解为：

$$
\mathbf w_{\mathrm{opt}}
=
\mathbf R^{-1}
\mathbf C
\left(
\mathbf C^H
\mathbf R^{-1}
\mathbf C
\right)^{-1}
\mathbf f
$$

这就是 CBF、MVDR、LCMV 和部分功率反演推导背后的共同数学骨架。
