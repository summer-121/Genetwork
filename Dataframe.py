import pandas as pd

# 원본 튜플 데이터
labels = (1, 0, 1, 0)
Gene1  = (6.6, 3.2, 5.6, 2.8)
Gene2  = (2.2, 3.3, 3.0, 2.9)
Gene3  = (1.0, 4.0, 1.2, 3.8)
Gene4  = (3.3, 3.2, 3.1, 3.2)
edges_tuples = [
    ('Gene1', 'Gene5'),
    ('Gene5', 'Gene4'),
    ('Gene2', 'Gene5'),
    ('Gene1', 'Gene4'),
    ('Gene2', 'Gene1'),
    ('Gene2', 'Gene3')
]

# 샘플 ID 생성 (S1, S2, S3, S4 ...)
samples = [f"S{i+1}" for i in range(len(labels))]

# expr_df 생성 (샘플 x 유전자)
expr_df = pd.DataFrame({
    "GENE1": Gene1,
    "GENE2": Gene2,
    "GENE3": Gene3,
    "GENE4": Gene4
}, index=samples)

# labels_df 생성 (sample, label)
labels_df = pd.DataFrame({
    "sample": samples,
    "label": labels
})

# edges_df 생성
edges_df = pd.DataFrame(edges_tuples, columns=["source", "target"])

print("expr_df:")
print(expr_df, "\n")
print("labels_df:")
print(labels_df, "\n")
print("edges_df:")
print(edges_df)
