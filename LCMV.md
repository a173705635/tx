# 第十九章 LCMV 线性约束最小方差波束形成

## 19.1 学习目标

本章解决以下问题：

- 为什么 MVDR 只能表达一个约束
- 如何同时保护多个方向
- 如何显式规定零陷方向
- LCMV 的矩阵约束如何构造
- LCMV 权值如何推导
- 多约束为什么会消耗阵列自由度
- GNSS 多卫星为什么难以共用一套小阵列权值

---

## 19.2 从 MVDR 到 LCMV

MVDR 只有一个无失真约束：

$$
\mathbf{w}^{H}\mathbf{a}_0
=
1
$$

如果希望同时满足多个线性约束，可以写成：

$$
\boxed{
\mathbf{C}^{H}\mathbf{w}
=
\mathbf{f}
}
$$

约束矩阵为：

$$
\mathbf{C}
=
\begin{bmatrix}
\mathbf{c}_1
&
\mathbf{c}_2
&
\cdots
&
\mathbf{c}_Q
\end{bmatrix}
$$

约束响应向量为：

$$
\mathbf{f}
=
\begin{bmatrix}
f_1 \\
f_2 \\
\vdots \\
f_Q
\end{bmatrix}
$$

---

## 19.3 LCMV 优化问题

LCMV 的优化问题为：

$$
\boxed{
\min_{\mathbf{w}}
\mathbf{w}^{H}\mathbf{R}_x\mathbf{w}
}
$$

约束为：

$$
\boxed{
\mathbf{C}^{H}\mathbf{w}
=
\mathbf{f}
}
$$

MVDR 是 LCMV 的单约束特例：

$$
\mathbf{C}
=
\mathbf{a}_0
$$

$$
\mathbf{f}
=
1
$$

---

## 19.4 LCMV 权值推导

构造拉格朗日函数：

$$
\mathcal{L}
=
\mathbf{w}^{H}\mathbf{R}_x\mathbf{w}
+
\boldsymbol{\lambda}^{H}
\left(
\mathbf{C}^{H}\mathbf{w}-\mathbf{f}
\right)
+
\left(
\mathbf{w}^{H}\mathbf{C}-\mathbf{f}^{H}
\right)
\boldsymbol{\lambda}
$$

对权向量求导并令其为零：

$$
\mathbf{R}_x\mathbf{w}
+
\mathbf{C}\boldsymbol{\lambda}
=
\mathbf{0}
$$

因此：

$$
\mathbf{w}
=
-\mathbf{R}_x^{-1}
\mathbf{C}\boldsymbol{\lambda}
$$

代入约束：

$$
\mathbf{C}^{H}\mathbf{w}
=
\mathbf{f}
$$

得到：

$$
\boldsymbol{\lambda}
=
-
\left(
\mathbf{C}^{H}\mathbf{R}_x^{-1}\mathbf{C}
\right)^{-1}
\mathbf{f}
$$

最终：

$$
\boxed{
\mathbf{w}_{\mathrm{LCMV}}
=
\mathbf{R}_x^{-1}\mathbf{C}
\left(
\mathbf{C}^{H}\mathbf{R}_x^{-1}\mathbf{C}
\right)^{-1}
\mathbf{f}
}
$$

---

## 19.5 多目标保护约束

若需要同时保护两个方向：

$$
\theta_1
$$

和：

$$
\theta_2
$$

可构造：

$$
\mathbf{C}
=
\begin{bmatrix}
\mathbf{a}(\theta_1)
&
\mathbf{a}(\theta_2)
\end{bmatrix}
$$

$$
\mathbf{f}
=
\begin{bmatrix}
1 \\
1
\end{bmatrix}
$$

则 LCMV 保证：

$$
\mathbf{w}^{H}\mathbf{a}(\theta_1)
=
1
$$

$$
\mathbf{w}^{H}\mathbf{a}(\theta_2)
=
1
$$

同时在剩余自由度中最小化输出功率。

---

## 19.6 显式零陷约束

若希望保护目标方向：

$$
\theta_0
$$

并在已知干扰方向：

$$
\theta_i
$$

形成严格零陷，可构造：

$$
\mathbf{C}
=
\begin{bmatrix}
\mathbf{a}(\theta_0)
&
\mathbf{a}(\theta_i)
\end{bmatrix}
$$

$$
\mathbf{f}
=
\begin{bmatrix}
1 \\
0
\end{bmatrix}
$$

于是：

$$
\mathbf{w}^{H}\mathbf{a}(\theta_0)
=
1
$$

$$
\mathbf{w}^{H}\mathbf{a}(\theta_i)
=
0
$$

这种方式要求干扰方向已知或已估计。

---

## 19.7 主瓣宽化与导数约束

若只约束一个精确方向，目标方向稍有误差就可能出现增益下降。

可以在目标附近增加多个保护点：

$$
\theta_0-\Delta
$$

$$
\theta_0
$$

$$
\theta_0+\Delta
$$

构造：

$$
\mathbf{C}
=
\begin{bmatrix}
\mathbf{a}(\theta_0-\Delta)
&
\mathbf{a}(\theta_0)
&
\mathbf{a}(\theta_0+\Delta)
\end{bmatrix}
$$

并令：

$$
\mathbf{f}
=
\begin{bmatrix}
1 \\
1 \\
1
\end{bmatrix}
$$

也可以使用导向矢量对角度的导数约束：

$$
\mathbf{w}^{H}
\frac{\partial\mathbf{a}(\theta)}
{\partial\theta}
\bigg|_{\theta=\theta_0}
=
0
$$

它表示目标方向附近的一阶响应变化较小，从而提高方向误差容忍度。

---

## 19.8 自由度与约束数量

N 阵元复权向量具有 N 个复数维度。

若有 Q 个线性独立约束，满足约束后的可调整自由度大约为：

$$
\boxed{
N-Q
}
$$

因此：

- 约束越多，可用于自适应抑制干扰的自由度越少
- 约束数量接近阵元数时，权值几乎被完全确定
- 独立约束数量超过阵元数时，通常无法同时满足全部约束

一个目标无失真约束时：

$$
Q=1
$$

剩余自适应自由度约为：

$$
N-1
$$

这也是常说 N 阵元阵列通常最多具有约 N 减 1 个独立零陷自由度的原因之一。

---

## 19.9 GNSS 多卫星场景中的问题

GNSS 接收机同时可见多颗卫星。

每颗卫星方向不同，对应不同导向矢量：

$$
\mathbf{a}(\theta_1),
\mathbf{a}(\theta_2),
\cdots,
\mathbf{a}(\theta_S)
$$

如果使用一套公共阵列权值，并要求所有卫星方向均保持单位增益，则：

$$
\mathbf{C}
=
\begin{bmatrix}
\mathbf{a}(\theta_1)
&
\mathbf{a}(\theta_2)
&
\cdots
&
\mathbf{a}(\theta_S)
\end{bmatrix}
$$

$$
\mathbf{f}
=
\mathbf{1}
$$

当卫星数量大于或接近阵元数时，小阵列无法严格满足全部独立约束，同时也几乎没有剩余自由度抑制干扰。

更实际的方案包括：

- 每颗卫星使用独立权值
- 后相关波束形成
- 只保护少数关键卫星或方向区域
- 使用宽主瓣或不确定集约束
- 在阵列自由度和保护目标之间折中

---

## 19.10 数值实现注意事项

LCMV 公式中包含逆矩阵。

代码中应使用线性方程求解。

先计算：

$$
\mathbf{V}
=
\mathbf{R}^{-1}\mathbf{C}
$$

代码：

```python
V = np.linalg.solve(R_loaded, C)
```

再计算：

$$
\mathbf{G}
=
\mathbf{C}^{H}\mathbf{V}
$$

求解：

$$
\mathbf{G}\mathbf{q}
=
\mathbf{f}
$$

最后：

$$
\mathbf{w}
=
\mathbf{V}\mathbf{q}
$$

代码：

```python
V = np.linalg.solve(R_loaded, C)
q = np.linalg.solve(C.conj().T @ V, f)
w = V @ q
```

还需要检查：

- 约束矩阵是否列满秩
- 约束方向是否过于接近
- 约束数量是否超过阵元能力
- 协方差矩阵条件数是否过大
- 是否需要对角加载

---

## 19.11 LCMV 常见疑惑

### 疑惑一：LCMV 是否必须知道干扰方向

不一定。

若只给出保护约束，LCMV 可以利用协方差自动抑制其他高能量方向。

只有设置显式零陷约束时，才需要输入对应干扰方向。

### 疑惑二：LCMV 与 MVDR 的区别是什么

MVDR 只有一个保护约束。

LCMV 可以设置多个线性约束。

MVDR 是 LCMV 的特例。

### 疑惑三：约束越多是否越安全

不一定。

更多约束可以保护更多方向，但会消耗更多自由度，可能使干扰抑制能力下降，也可能使矩阵求解变得不稳定。

### 疑惑四：四阵元能否保护四颗卫星并同时抑制干扰

若四个保护约束线性独立：

$$
N=4
$$

$$
Q=4
$$

剩余自由度约为零。

权值基本被约束完全确定，几乎没有额外自由度自适应抑制未知干扰。

### 疑惑五：多个保护方向很接近会怎样

对应导向矢量可能高度相关，使约束矩阵条件数变差。

小误差会导致很大的权值变化，需要使用加载、约束合并或宽主瓣设计。

---

# 第二十章 本阶段概念关系总结

## 20.1 从测向到抗干扰的逻辑链

```text
阵列接收数据
    |
    v
样本协方差矩阵
    |
    +-----------------------------+
    |                             |
    v                             v
特征分解与子空间              最小输出功率优化
    |                             |
    v                             v
MUSIC / ESPRIT              MVDR / LCMV
    |                             |
    v                             v
估计来波方向                 形成波束和零陷
```

MUSIC 和 ESPRIT 主要回答：

> 信号从什么方向来。

MVDR 和 LCMV 主要回答：

> 已知需要保护的方向后，怎样组合阵元信号以抑制干扰。

---

## 20.2 MUSIC 与 ESPRIT 的共同基础

二者都依赖：

- 多信号阵列模型
- 样本协方差矩阵
- 特征值分解
- 信号子空间
- 阵列流形结构
- 信号数估计

MUSIC 利用：

$$
\mathbf{E}_n^{H}\mathbf{a}(\theta)
\approx
\mathbf{0}
$$

ESPRIT 利用：

$$
\mathbf{E}_{s2}
\approx
\mathbf{E}_{s1}\mathbf{\Psi}
$$

MUSIC 使用噪声子空间正交性。

ESPRIT 使用信号子空间平移不变性。

---

## 20.3 相干信号位于测向算法与实际传播之间

经典 MUSIC 和 ESPRIT 默认信号协方差矩阵满秩。

多径和转发复制会使信号相干：

$$
\operatorname{rank}
\left(
\mathbf{R}_s
\right)
<M
$$

这会导致信号子空间维数小于物理来波数量。

空间平滑通过重叠子阵恢复有效秩，但需要牺牲孔径和阵列自由度。

---

## 20.4 MVDR 与 LCMV 的共同基础

二者都最小化输出功率：

$$
\mathbf{w}^{H}\mathbf{R}\mathbf{w}
$$

区别在于约束。

MVDR：

$$
\mathbf{w}^{H}\mathbf{a}_0
=
1
$$

LCMV：

$$
\mathbf{C}^{H}\mathbf{w}
=
\mathbf{f}
$$

因此：

$$
\boxed{
\text{MVDR}
\subset
\text{LCMV}
}
$$

---

## 20.5 测向与波束形成的协作关系

典型流程可以是：

1. 使用 MUSIC 或 ESPRIT 估计干扰方向
2. 根据卫星星历、姿态和阵列模型得到卫星导向矢量
3. 使用 MVDR 或 LCMV 设计权值
4. 对多阵元数据进行加权合并
5. 评估输出 SINR、零陷深度和卫星信号损失

但标准 MVDR 不一定需要先显式估计干扰方向。

它可以直接从协方差矩阵学习干扰空间。

MUSIC 或 ESPRIT 更适合需要显示干扰方向、进行态势感知或设置显式零陷约束的场景。

---

## 20.6 本阶段最关键的能力边界

### 阵元数必须大于信号数

MUSIC 需要噪声子空间：

$$
N-M>0
$$

### 相干信号会降低有效信号秩

即使物理上存在多个方向，统计上也可能只有一个独立信号维度。

### 空间平滑需要额外阵元

经典前向空间平滑分辨 M 个完全相干信号时，常见最低条件为：

$$
N\ge 2M
$$

### 多个约束会消耗波束形成自由度

剩余自适应自由度约为：

$$
N-Q
$$

### 导向矢量误差可能导致自消除

工程系统中的阵元位置误差、通道相位误差、互耦和姿态误差都会影响结果。

---

## 20.7 本阶段疑惑点汇总

### 疑惑一：多个信号时，特征向量是否一一对应导向矢量

通常不对应。

正确关系是信号子空间与阵列流形空间相同。

### 疑惑二：MUSIC 为什么要扫描

特征分解给出子空间，不直接给角度。

需要将候选角度映射为导向矢量，再检测其与噪声子空间的正交性。

### 疑惑三：ESPRIT 中的旋转是什么

是两个平移子阵之间的复相位变换。

Psi 不是直接的物理旋转矩阵，但其特征值保留真实空间相位。

### 疑惑四：为什么相干信号导致秩下降

多个信号波形之间存在固定线性关系，无法形成多个独立统计维度。

### 疑惑五：两个相干信号为什么通常至少需要四阵元

空间平滑既需要足够长的子阵保留噪声子空间，也需要足够多的重叠子阵恢复信号秩。

### 疑惑六：MVDR 为什么能自动形成零陷

强干扰在协方差矩阵中贡献较大。

最小输出功率优化会降低权值对这些空间模式的响应。

### 疑惑七：MVDR 为什么可能把目标压掉

约束保护的是假设导向矢量。

若真实导向矢量与假设不一致，目标能量可能被优化器当成需要抑制的输出功率。

### 疑惑八：仿真中为什么还要构造理论目标和干扰协方差

这些矩阵可用于独立计算输出目标功率、干扰功率和 SINR。

它们不一定参与 MVDR 权值设计。

### 疑惑九：LCMV 为什么不能无限增加保护方向

每个独立约束都会占用一个自由度。

约束过多后，阵列没有足够自由度自适应抑制未知干扰。

### 疑惑十：投影公式为什么是加号

方向单位向量为：

$$
\mathbf{u}
=
\begin{bmatrix}
\cos\theta \\
\sin\theta
\end{bmatrix}
$$

二维向量为：

$$
\mathbf{x}
=
\begin{bmatrix}
x_1 \\
x_2
\end{bmatrix}
$$

标量投影为内积：

$$
z
=
\mathbf{u}^{T}\mathbf{x}
$$

所以：

$$
\boxed{
z
=
x_1\cos\theta
+
x_2\sin\theta
}
$$

这与主动旋转向量不是同一个操作。

---

## 20.8 下一阶段衔接

下一阶段适合继续学习：

1. GSC 广义旁瓣相消器
2. LMS 和 NLMS 自适应更新
3. RLS 和递推协方差估计
4. 宽带阵列模型
5. 空时自适应处理
6. GNSS 前相关与后相关抗干扰
7. 阵列标定
8. 稳健波束形成

其中 GSC 是 LCMV 的等价实现结构，可以把抽象约束优化拆分成固定保护支路、阻塞矩阵和自适应相消支路。

---
