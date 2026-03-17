#!/bin/bash
#SBATCH -J REMD_GBP
#SBATCH -p gpu_4l
#SBATCH -N 1
#SBATCH --gres=gpu:4
#SBATCH -o ys_%j.out
#SBATCH -e ys_%j.err
#SBATCH -A csong_g1
#SBATCH --qos=csongg4c
#SBATCH --exclusive

echo "Running on node: $(hostname)"

# Store the job ID in a variable
JOB_ID=$SLURM_JOB_ID
echo "Job ID is: $JOB_ID"

source /appsnew/source/intel2024.sh
source /appsnew/source/cuda-12.6.2.sh
source /home/csong_pkuhpc/lustre3/ys/SoftWare/miniconda3/bin/activate

export PATH=/appsnew/usr/cuda/cuda-12.6.2/bin:/appsnew/usr/intel/intel2024/vtune/2024.0/bin64:/appsnew/usr/intel/intel2024/vtune/2024.0/bin64:/appsnew/usr/intel/intel2024/mkl/2024.0/bin:/appsnew/usr/intel/intel2024/dpcpp-ct/2024.0/bin:/appsnew/usr/intel/intel2024/dev-utilities/2024.0/bin:/appsnew/usr/intel/intel2024/debugger/2024.0/opt/debugger/bin:/appsnew/usr/intel/intel2024/compiler/2024.0/opt/oclfpga/bin:/appsnew/usr/intel/intel2024/compiler/2024.0/bin:/appsnew/usr/intel/intel2024/advisor/2024.0/bin64:/home/csong_pkuhpc/lustre3/ys/SoftWare/miniconda3/bin:/home/csong_pkuhpc/lustre3/ys/SoftWare/miniconda3/condabin:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/csong_pkuhpc/.dotnet/tools:/usr/lpp/mmfs/bin:/rm/rm_prog/slurm/18.08.7/sbin:/rm/rm_prog/slurm/18.08.7/bin:/data01/oldbin/newbin:/rm/rm_prog/slurm/18.08.7/sbin:/rm/rm_prog/slurm/18.08.7/bin:/data01/oldbin/newbin:/home/csong_pkuhpc/Leon-test/gauopen:/home/csong_pkuhpc/Leon-test/g16/bsd:/home/csong_pkuhpc/Leon-test/g16:/home/csong_pkuhpc/bin

conda activate ctgomartini_REMD_multiGPUs3

python ../build_mpirun_configfile.py "python run_REMD.py"
#mpiexec.hydra -f hostfile -configfile configfile
mv hostfile hostfile1
mv configfile configfile1

python ../build_mpirun_configfile.py "python run_REMD_restart.py"
#mpiexec.hydra -f hostfile -configfile configfile
mv hostfile hostfile2
mv configfile configfile2


#TASK1_CMD="mpiexec.hydra -f hostfile1 -configfile configfile1"
TASK1_CMD="mpiexec.hydra -f hostfile2 -configfile configfile2"
TASK2_CMD="mpiexec.hydra -f hostfile2 -configfile configfile2"

# 任务超时时间（2.5小时）
TIMEOUT=9000

# 标记是否是第一次执行
FIRST_RUN=true

# 循环执行直到任务在5小时内完成
while true; do
    # 根据是否是第一次运行选择命令
    if [ "$FIRST_RUN" = true ]; then
        TASK_CMD=$TASK1_CMD
        FIRST_RUN=false
    else
        TASK_CMD=$TASK2_CMD
    fi

    # 启动任务并记录PID
    $TASK_CMD &
    TASK_PID=$!

    # 启动超时计时器
    sleep $TIMEOUT &
    TIMER_PID=$!

    # 等待任意一个进程结束
    wait -n $TASK_PID $TIMER_PID
    EXIT_STATUS=$?

    # 检查任务是否超时
    if kill -0 $TASK_PID 2>/dev/null; then
        # 任务仍在运行（超时）
        echo "任务超时，正在终止并重启..."
        kill $TASK_PID 2>/dev/null       # 终止任务进程
        wait $TASK_PID 2>/dev/null       # 清理僵尸进程
        kill $TIMER_PID 2>/dev/null      # 终止睡眠进程
    else
        # 任务在超时前完成
        kill $TIMER_PID 2>/dev/null      # 终止睡眠进程
        wait $TIMER_PID 2>/dev/null      # 清理睡眠进程
        exit $EXIT_STATUS                # 退出并返回任务状态码
    fi
done
