### 功能

- 讨论MEKF方法

### 修改记录

- 增加主函数test_mian.m

### Usage Guideline

1. use **genTrig()** or **genRotAxis()** to generate reference motions and angular velocity measurements.
2. use **genMea()** to generate attitude or vector measurements for the reference motion.
3. use the functions in **Filters Without Bias** folder to estimate attitude.
