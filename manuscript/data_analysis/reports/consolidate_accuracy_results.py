#!/usr/bin/env python3
"""
合并准确性评估结果脚本
将各个数据集的准确性评估结果合并为统一的汇总文件
"""

import json
import sys
from pathlib import Path

def consolidate_accuracy_results():
    """合并准确性评估结果"""
    
    # 定义路径
    accuracy_dir = Path("manuscript/data_analysis/results/accuracy_assessment")
    output_file = accuracy_dir / "accuracy_assessment_summary.json"
    
    # 查找所有准确性评估文件
    evaluation_files = list(accuracy_dir.glob("accuracy_evaluation_*.json"))
    
    if not evaluation_files:
        print("未找到准确性评估文件")
        return False
    
    # 读取所有评估结果
    dataset_details = []
    
    for eval_file in evaluation_files:
        try:
            with open(eval_file, 'r') as f:
                data = json.load(f)
                dataset_details.append(data)
                print(f"已加载: {eval_file.name}")
        except Exception as e:
            print(f"加载 {eval_file} 时出错: {e}")
            continue
    
    if not dataset_details:
        print("未成功加载任何评估结果")
        return False
    
    # 创建汇总结构
    summary = {
        "metadata": {
            "generation_time": dataset_details[0].get("timestamp", ""),
            "total_datasets": len(dataset_details),
            "evaluation_version": "1.0.0"
        },
        "dataset_details": dataset_details,
        "overall_statistics": {
            "avg_sensitivity": sum(d["gene_metrics"]["sensitivity"] for d in dataset_details) / len(dataset_details),
            "avg_precision": sum(d["gene_metrics"]["precision"] for d in dataset_details) / len(dataset_details),
            "avg_f1_score": sum(d["gene_metrics"]["f1_score"] for d in dataset_details) / len(dataset_details),
            "total_reference_genes": sum(d["gene_metrics"]["num_reference_genes"] for d in dataset_details),
            "total_predicted_genes": sum(d["gene_metrics"]["num_predicted_genes"] for d in dataset_details),
            "total_true_positives": sum(d["gene_metrics"]["true_positives"] for d in dataset_details)
        }
    }
    
    # 保存汇总结果
    try:
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        print(f"汇总结果已保存到: {output_file}")
        return True
    except Exception as e:
        print(f"保存汇总结果时出错: {e}")
        return False

if __name__ == "__main__":
    success = consolidate_accuracy_results()
    sys.exit(0 if success else 1)