import RNA
import pandas as pd
from collections import Counter
#train col 18 et 16, 2400 ind et 100%homlogie
# test col 6 et 5, 3634ind et 49% homologie
data = pd.read_json('test.json', lines = True)
ss_list = data.loc[:,"structure"]
seq_list = data.loc[:,"sequence"]


id_list = []
dif_list = []
ss_id = []
ss_dif = []

with open("data/vienna_test", "r") as fillout:
	for i, seq in enumerate(seq_list):
		(ss, mfe) = RNA.fold(seq)
		if ss == ss_list[i]:
			id_list.append(i)
			ss_id.append(ss)
		else:
			x=0
			for j in range(len(seq)):
				if ss_list[i][j] == ss[j]:
					x += 1
			x = round(x/len(ss),2)
			dif_list.append(x)

			ss_dif.append(ss)
		fillout.write(str(ss) + "\n")
	select = [i for i in range(len(dif_list)) if dif_list[i] > 0.8]
	print(len(select)/len(dif_list))
	print(Counter(dif_list).most_common(50))


data2 = pd.read_json('train.json', lines = True)
ss_list2 = data2.loc[:,"structure"]
seq_list2 = data2.loc[:,"sequence"]
id_list2 = []
dif_list2 = []
ss_id2 = []
ss_dif2 = []

with open("data/vienna_train", "r") as fillout:
	for i, seq in enumerate(seq_list2):
		(ss, mfe) = RNA.fold(seq)
		if ss == ss_list2[i]:
			id_list2.append(i)
			ss_id2.append(ss)
		else:
			dif_list2.append(i)
			ss_dif2.append(ss)
		fillout.write(str(ss) + "\n")

print("train : {}\t test : {}".format(len(id_list2)/len(ss_list2), len(id_list)/len(ss_list)))