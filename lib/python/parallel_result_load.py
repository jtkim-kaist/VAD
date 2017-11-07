import pickle


with open('/home/sbie/storage2/result/result_proposed_soft.p', 'rb') as f:
    result = pickle.load(f)

print("result_load")