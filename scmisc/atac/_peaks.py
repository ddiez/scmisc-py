def as_sparse(df, row_column, col_column, val_column):
  import pandas as pd
  from scipy.sparse import csr_matrix

  if isinstance(row_column, int):
    row_column = df.iloc[:, row_column]
  else:
    row_column = df[row_column]
  
  if isinstance(col_column, int):
    col_column = df.iloc[:, col_column]
  else:
    col_column = df[col_column]
  
  if isinstance(val_column, int):
    val_column = df.iloc[:, val_column]
  else:
    val_column = df[val_column]

  # Factorize the row and col names to get numerical indices
  row_indices, row_labels = pd.factorize(row_column)
  col_indices, col_labels = pd.factorize(col_column)

  # Convert the DataFrame to a sparse matrix
  values = val_column.values
  
  X = csr_matrix((values, (row_indices, col_indices)))

  return {
    "X": X,
    "row_labels": row_labels,
    "col_labels": col_labels
  }

def range_to_list(x):
  return x.replace(":", "-").split("-")



# def range_to_list(ranges, index):
#   if isinstance(index, int):
#     chr, start, end = ranges.iloc[index]
  
#   if isinstance(index, str):
#     chr, start, end = ranges.loc[index]
    
#   return chr, start, end

def get_query(tb, range):
  #chr, start, end = range_to_list(range)
  #if len(range) == 3:
  #  chr, start, end = range
  #  out = tb.query(chr, int(start), int(end))

  if isinstance(range, str):
    out = tb.querys(range)
  
  if isinstance(range, int):
    out = tb.queryi(range)

  return out


def get_interval_barcodes(tb, range, whitelist=None):
  import pandas as pd
  import numpy as np
  from collections import Counter

  out = get_query(tb, range)

  barcode = [ i[3] for i in out ]

  if whitelist is not None:
    barcode = np.array(barcode)
    barcode = barcode[np.isin(barcode, whitelist)]

  barcode = Counter(barcode)
  barcode = pd.DataFrame({
    "barcode": list(barcode.keys()),
    "count": list(barcode.values())
  })
  
  barcode.insert(0, "peak", range)

  return barcode

def count_intervals(tb, ranges, whitelist=None):
  import pandas as pd

  res = []
  
  for range in ranges:
    tmp = get_interval_barcodes(tb, range, whitelist=whitelist)
    res.append(tmp)
  
  res = pd.concat(res)
  
  return res

def feature_matrix(tb, ranges, whitelist=None):
  res = count_intervals(tb, ranges, whitelist=whitelist)
  res = as_sparse(res, "barcode", "peak", "count")
  return res