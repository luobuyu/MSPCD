from fileinput import filename
import os
import platform
import random
import sys

exe_path = "../x64/Debug/"
instance_path = "../Deploy/Instance/Bench/"
solution_path = "../Deploy/Solution/Bench/"
log_path = "../Deploy/Solution/log.csv"
help = """
1. 测试所有样例
./test.py all
2. 程序跑10遍
./test.py 10
3. 测试单个样例
./test.py instance

4. 单个样例测10遍
./test.py instance 10
5. 所有样例每个测10遍
./test.py all 10

"""

# 加载所有的样例
def load_all_instance():
    instances = os.listdir(instance_path)
    instances.sort()
    return instances

# 加载可执行文件
def load_exe():
    for name in os.listdir(exe_path):
        if name[-4:] == '.exe':
            return exe_path + name

# 运行单个样例  exe   timeout   seed      instance_name
def run_single_instance_with_times(file_name, times = 1):
    timeout = 10
    for i in range(times):
        r = random.randint(0, 65536)
        order = load_exe() + ' ' + str(timeout) + ' ' + str(r) + ' ' + file_name[:file_name.rfind('.')] + ' < ' + instance_path + file_name + ' > ' + solution_path + file_name
        print(file_name, "\t第" + str(i + 1) + "次运行")
        os.system(order)

# 运行所有的样例  exe timeout seed  instance_name
def run_all_instance_with_times(times = 1):
    for instance in load_all_instance():
        run_single_instance_with_times(instance, times)


# 只运行可执行程序 exe timeout seed
def run_exe_times(times = 1):
    timeout = 10
    for i in range(times):
        r = random.randint(0, 65536)
        print("第" + str(i + 1) + "次")
        order = load_exe() + ' ' + str(timeout) + ' ' + str(r)
        print(order)
        os.system(order)

def check_platform():
    global exe_path, instance_path, solution_path, log_path
    plat = platform.system()
    if plat == 'Windows':
        exe_path = exe_path.replace('/', '\\')
        instance_path = instance_path.replace('/', '\\')
        solution_path = solution_path.replace('/', '\\')
        log_path = log_path.replace('/', '\\')


if __name__ == '__main__':
    check_platform()
    if len(sys.argv) == 2: # 测试一遍
        if sys.argv[1] == 'all':
            run_all_instance_with_times()       # ./test.py all
        elif sys.argv[1].isdigit():
            run_exe_times(int(sys.argv[1]))     # ./test.py 10
        else:
            run_single_instance_with_times(sys.argv[1])     # ./test.py instance_name
    elif len(sys.argv) == 3:
        if sys.argv[1] == 'all':
            run_all_instance_with_times(int(sys.argv[2]))   # ./test.py all 10
        else:
            run_single_instance_with_times(sys.argv[1], int(sys.argv[2]))   # ./test.py instance_name 10
    else:
        print(help)
