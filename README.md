cpp代码符合GSL规范：类与结构体基本用std实现，使用封装语句避免if分支嵌套（错误尽早退回）与for判断
apps/下是可执行函数的脚本，cpp默认的main.cpp
中间的链接文件放 build/
执行文件放 (built-in) bin/
核心cpp放入 src/core/
外接函数实现放 src/功能名
src/${project}/下是python的核心脚本，实现某些功能
include/文件下面，${project}名下是核心实现的头文件与注册表头文件
功能名函数的声明头文件在include/功能名/下
大量Shell脚本与python脚本放入tools/，封装了各种功能
测试脚本放入tests/
python脚本的测试推荐使用uv run python，单独的uv run有导入库的错误
高级抽象的交给python ，底层并发实现用cpp



