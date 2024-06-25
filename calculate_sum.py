import pandas as pd
import sys

def calculate_sum(input_file, output_file):
    # 从文件中读取数据，不使用第一行作为列标题
    data = pd.read_csv(input_file, delimiter='\t', header=None)

    # 获取行数
    n_rows = data.shape[0]
    
    # 检查最后一个单元的行数
    last_group_size = n_rows % 50
    if last_group_size > 0 and n_rows > 50: # 小于50则合并到前一个单元
        # 使最后两个单元合并
        n_rows = n_rows - last_group_size
        data = data.iloc[:n_rows, :]

    # 使用groupby和agg函数计算每个单元的和
    sum_result = data.groupby(data.index // 50).agg('sum')

    # 计算每个单元的行数
    count_result = data.groupby(data.index // 50).agg('count').iloc[:,0]

    # 计算每个单元的平均值
    result = sum_result.div(count_result, axis=0)

    # 如果有不足50行的单元，将其加到最后一个完整单元上
    if last_group_size > 0 and n_rows > 50:
        temp_sum = sum_result.copy()
        temp_sum.loc[temp_sum.shape[0]] = data.iloc[n_rows:, :].sum()
        temp_count = count_result.copy()
        temp_count.loc[temp_count.shape[0]] = last_group_size
        temp_result = temp_sum.div(temp_count, axis=0)
        result.iloc[-1, :] = temp_result.iloc[-2:, :].sum(axis=0)
    elif last_group_size > 0 and n_rows <= 50:
        result = sum_result.div(last_group_size, axis=0)

    # 输出并保存结果
    # print(result)  # if you don't want to see the result on the screen, comment out this line
    result.to_csv(output_file, index=False, sep='\t', header=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python calculate_sum.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    calculate_sum(input_file, output_file)
