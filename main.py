#!/usr/bin/env python3
import sys
import os
import subprocess
import glob
import argparse

def need_recompile(exe_file, source_dir, debug=False):
    """检查是否需要重新编译：比较源文件和可执行文件的修改时间"""
    if not os.path.exists(exe_file):
        return True

    exe_mtime = os.path.getmtime(exe_file)

    # 检查所有 .cpp 和 .h 文件
    for pattern in ["*.cpp", "*.h"]:
        for src_file in glob.glob(os.path.join(source_dir, pattern)):
            if os.path.getmtime(src_file) > exe_mtime:
                return True

    return False

def parse_perf_output(perf_output, exe_name="zbuffer"):
    """解析 perf report 输出，提取热点函数信息，只保留用户代码"""
    lines = perf_output.strip().split('\n')
    results = []

    for line in lines:
        line = line.strip()
        # 跳过注释行和空行
        if not line or line.startswith('#'):
            continue

        parts = line.split()
        if len(parts) >= 5:
            try:
                overhead = parts[0].rstrip('%')
                if overhead.replace('.', '').isdigit():
                    overhead_val = float(overhead)
                    # 函数名通常在最后一个字段
                    func_name = parts[-1]
                    # 共享对象名在第 4 个位置 (index 3)
                    shared_obj = parts[3] if len(parts) > 3 else ""

                    # 只保留用户代码（zbuffer 可执行文件中的函数）
                    # 排除 C 标准库和系统库
                    if exe_name in shared_obj:
                        results.append({
                            'overhead': overhead_val,
                            'func': func_name,
                            'shared_obj': shared_obj
                        })
            except (ValueError, IndexError):
                continue

    return results

def print_profile_results(results, top_n=15):
    """打印格式化的性能分析结果（仅用户代码）"""
    if not results:
        print("\n未检测到用户代码中的热点函数")
        return

    # 计算用户代码总耗时
    total_user_overhead = sum(r['overhead'] for r in results)

    print("\n" + "=" * 70)
    print("           性能分析结果 - 用户代码热点 (Performance Profile)")
    print("=" * 70)
    print(f"{'绝对耗时':<12} {'相对耗时':<12} {'函数名':<40}")
    print("-" * 70)

    for i, r in enumerate(results[:top_n]):
        abs_overhead = f"{r['overhead']:.2f}%"
        # 相对于用户代码的占比
        rel_overhead = f"{r['overhead'] / total_user_overhead * 100:.1f}%" if total_user_overhead > 0 else "N/A"
        # 函数名（demangle 后可能很长）
        func_name = r['func']
        if len(func_name) > 39:
            func_name = func_name[:36] + "..."
        print(f"{abs_overhead:<12} {rel_overhead:<12} {func_name:<40}")

    print("=" * 70)
    print(f"\n用户代码总耗时: {total_user_overhead:.2f}% (排除 C 标准库)")
    print(f"显示函数数: {min(top_n, len(results))} / {len(results)}")

def main():
    parser = argparse.ArgumentParser(description='Z-Buffer 渲染器')
    parser.add_argument('algorithm', help='算法类型: scanline 或 interval')
    parser.add_argument('scene_folder', help='场景文件夹路径')
    parser.add_argument('--profile', default=False, 
                        help='启用性能分析（使用 perf）')
    parser.add_argument('--profile-top', type=int, default=15,
                        help='显示前 N 个热点函数（默认 15）')

    args = parser.parse_args()

    algorithm = args.algorithm
    scene_folder = args.scene_folder
    enable_profile = args.profile
    profile_top = args.profile_top

    if not os.path.isdir(scene_folder):
        print(f"Error: Scene folder '{scene_folder}' does not exist")
        sys.exit(1)

    obj_file = os.path.join(scene_folder, "scene.obj")
    xml_file = os.path.join(scene_folder, "scene.xml")

    if not os.path.exists(obj_file):
        print(f"Error: OBJ file '{obj_file}' does not exist")
        sys.exit(1)

    if not os.path.exists(xml_file):
        print(f"Error: XML file '{xml_file}' does not exist")
        sys.exit(1)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    source_dir = os.path.join(script_dir, algorithm)
    exe_file = os.path.join(source_dir, "zbuffer")

    # 性能分析需要调试符号
    compile_flags = ["-O2", "-std=c++17", "-mbmi2"]
    if enable_profile:
        compile_flags.append("-g")  # 添加调试符号以便 perf 显示函数名

    # 检查是否需要编译
    if need_recompile(exe_file, source_dir):
        cpp_files = glob.glob(os.path.join(source_dir, "*.cpp"))

        print(f"Compiling {algorithm} Z-buffer...")
        if enable_profile:
            print("(包含调试符号用于性能分析)")
        compile_cmd = ["g++"] + compile_flags + ["-I", source_dir, "-o", exe_file] + cpp_files
        result = subprocess.run(compile_cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print("Compilation failed:")
            print(result.stderr)
            sys.exit(1)

        print("Compilation successful!")
    else:
        print("No changes detected, skipping compilation.")

    # Run
    print(f"\nRunning {algorithm} Z-buffer algorithm on '{scene_folder}'...")

    perf_data_file = os.path.join(script_dir, "perf.data")

    if enable_profile:
        # 使用 perf record 运行程序
        print("(性能分析模式: 使用 perf 记录执行信息)")
        run_cmd = ["perf", "record", "-g", "-o", perf_data_file, "--", exe_file, scene_folder]
    else:
        run_cmd = [exe_file, scene_folder]

    process = subprocess.Popen(run_cmd)

    try:
        process.wait()

        if process.returncode != 0:
            print("Execution failed!")
            sys.exit(1)

        print(f"\nOutput saved to: {os.path.join(scene_folder, 'scene.bmp')}")

        # 如果启用了性能分析，生成并显示报告
        if enable_profile and os.path.exists(perf_data_file):
            print("\n正在生成性能分析报告...")

            # 使用 perf report 生成报告
            # --dso 只显示指定可执行文件中的符号
            # --demangle 显示 C++ 函数的可读名称
            report_cmd = [
                "perf", "report",
                "-i", perf_data_file,
                "--stdio",
                "--no-children",
                "-n",
                "--dso", exe_file,  # 只显示 zbuffer 可执行文件中的函数
                "--demangle",       # C++ 函数名反修饰
                "--percent-limit", "0.1"  # 降低阈值以捕获更多用户函数
            ]

            report_result = subprocess.run(report_cmd, capture_output=True, text=True)

            if report_result.returncode == 0:
                results = parse_perf_output(report_result.stdout, "zbuffer")
                if results:
                    print_profile_results(results, profile_top)
                else:
                    # 如果解析失败，直接显示原始输出
                    print("\n性能分析结果:")
                    print(report_result.stdout)
            else:
                print("性能分析报告生成失败:")
                print(report_result.stderr)

            # 清理 perf.data 文件
            try:
                os.remove(perf_data_file)
            except OSError:
                pass

    except KeyboardInterrupt:
        print("\nInterrupted! Killing zbuffer process...")
        process.kill()
        process.wait()
        sys.exit(130)

if __name__ == "__main__":
    main()
