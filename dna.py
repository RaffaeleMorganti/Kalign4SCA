import numpy as np
from scipy.stats import norm


def tracesToFile(traces, file = "out.afa"):
  """
  converts traces to DNA sequence according the method described
  parameters
    traces - a (n*l) numpy array containing n misaligned power traces of length l
    file - filename to store the output of the conversion
  """
   x,y = traces.shape
  letters = np.zeros((x,y-2),"S1")
  length = y-2
  t = traces
  letters[(t[:,:y-2] <= t[:,1:y-1])&(t[:,2:y] >= t[:,1:y-1])] = "U" # 123
  letters[(t[:,:y-2] <= t[:,2:y  ])&(t[:,2:y] <= t[:,1:y-1])] = "V" # 132
  letters[(t[:,:y-2] <= t[:,1:y-1])&(t[:,2:y] <= t[:,0:y-2])] = "W" # 231
  letters[(t[:,:y-2] >= t[:,1:y-1])&(t[:,2:y] <= t[:,1:y-1])] = "X" # 321
  letters[(t[:,:y-2] >= t[:,2:y  ])&(t[:,2:y] >= t[:,1:y-1])] = "Y" # 312
  letters[(t[:,:y-2] >= t[:,1:y-1])&(t[:,2:y] >= t[:,0:y-2])] = "Z" # 213
  # write fasta like file
  with open(file,"w") as f:
    for (index, string) in enumerate(letters.view(f"S{length}")):
      f.write(f"> {index}\n{string[0].decode()}\n")


def tracesFromFile(traces, file="aligned.afa", cutoff = 0.5):
  """
  align traces according to the dna-like alignment provided
  parameters
    traces - the same value used in tracesToFile()
    cutoff - discard time indices with higher ratio of missing
    file - filename of the aligned sequences
  returns - 2d numpy array of aligned traces
  """
  def str2array(string, trace):
    k,j = 0,0
    val = np.zeros(len(string) + 2,trace.dtype)
    for c in string:
      if c != '-':
        val[j:j+2] = trace[k:k+2]
        k += 1
      j += 1
    return val

  data = []
  # insert gaps
  with open(file,"r") as f:
    i, string = 0, ""
    for line in f:
      if line[0] == '>':
        data.append(str2array(string,traces[i]))
        i, string = int(line[1:]), ""
      else:
        string += line.strip()
    data.append(str2array(string,traces[i]))
  # remove columns below cutoff
  data = np.array(data[1:])
  data = data[:, (data==0).sum(0) < len(data) * cutoff]
  # fill missing values
  for i in range(2,data.shape[1]-1):
    x = np.all(data[:,i:i+2]==0, 1)
    data[x, i:i+2] = data[x, i-2:i]
  for i in range(1,data.shape[1]):
    data[data[:,i]==0,i] = data[data[:,i]==0,i-1]

  return data
