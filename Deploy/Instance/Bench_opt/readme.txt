The zip file Bench_opt.zip is containing task graph examples with
known optimal schedule length for given multiprocessor
architecture. The file names are ogra<n>_<r>_<p>.td, where <n>
should be substituted with the number of task in that file, <r>
with edge existence density, while <p> should be replaced with the
number of processors for which the optimal solution is given.


In each file, data are written in following format:

- first row contains the value for the number of tasks n;

- next n rows contain following data: task index i=1,...,n, task
duration t_i, i=1,...,n, number of successors n_{succ} and list of
n_{succ} pairs s_j  c_{ij} representing index of j-th successor
and corresponding communication.

Possible values for <p> are:

- 2 meaning scheduling onto 2-processor system;

- 4 represent 2-dimensional hypercube;

- 6 is used to describe 2X3 mesh of processors;

- 8 represent 3-dimensional hypercube;

- 9 is used to describe 3X3 mesh of processors;

- 12 is used to describe 3X4 mesh of processors;

- 16 represent 4-dimensional hypercube;

Moreover, one can define the completely connected multiprocessor
architecture with p processors, the length of optimal schedule
will not change.

The lengths of optimal schedules are given in the following table.
They depend only on n and not on the values for p and r (density).

---------------------------------------------------------------------
|n       | 50 | 100| 150 | 200 | 250 | 300 | 350 | 400 | 450 |500   |
---------------------------------------------------------------------
|SL_{opt}| 600| 800| 1000| 1200| 1400| 1600| 1800| 2000| 2200| 2400 |
---------------------------------------------------------------------

压缩文件Bench_opt.zip包含了任务图的例子，其中有已知的多处理器最佳时间表长度。
含有已知的最佳计划长度的任务图，用于给定的多处理器
架构的最佳时间表。文件名称为ogra<n>_<r>_<p>.td，其中<n>
应该用该文件中的任务数来代替，<r>
用边缘存在密度代替，而<p>应改为
而<p>应替换为给出最佳解决方案的处理器数量。


在每个文件中，数据都以下列格式写入。

- 第一行包含任务数n的值。

- 接下来的n行包含以下数据：任务索引i=1,...,n, 任务
持续时间t_i, i=1,...,n, 继承者的数量n_{succ}和列表中的
n_{succ}对s_j c_{ij}代表第j个继承者的索引
和相应的通信。

<p>的可能值是。

- 2表示调度到2个处理器系统上。

- 4代表2维超立方体。

- 6用于描述2X3网格的处理器。

- 8代表3维超立方体。

- 9用于描述3X3网格的处理器。

- 12用于描述3X4网格的处理器。

- 16代表4维超立方体。

此外，我们可以定义完全连接的多处理器结构
此外，我们可以定义完全连接的多处理器架构，有p个处理器，最佳时间表的长度
将不会改变。

下表给出了最佳时间表的长度。
它们只取决于n，而不取决于p和r的值（密度）。

---------------------------------------------------------------------
|n | 50 | 100| 150 | 200 | 250 | 300 | 350 | 400 | 450 |500 |
---------------------------------------------------------------------
|SL_{opt}|600|800|1000|1200|1400|1600|1800|2000|2200|2400||
---------------------------------------------------------------------
