import csv
def checkPree(args):
    file_path=args.input_file
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            reader = csv.reader(file)
            rows = list(reader)
            
            # 检查是否为空文件或只有行头
            if len(rows) <= 1:
                return False
            
            # 获取每行的列数
            num_columns = len(rows[0])
            
            # 检查每一行的列数是否一致
            for row in rows:
                if len(row) != num_columns:
                    return False
            return True
    except Exception as e:
        return False

