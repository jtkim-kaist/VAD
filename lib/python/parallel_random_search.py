import VAD_Proposed as VR
import numpy as np
import tensorflow as tf
import pickle
from multiprocessing import Process, Queue

'''random search script'''

distribution_num = 7
test_num = 1
max_epoch = 101  # 351
gpu_0_append = 4


def get_parameter(min_val, max_val, shape):

    x = np.random.rand(shape[0], shape[1])
    x = (max_val - min_val) * x + min_val

    return x


def vad_fnc(p_initLr, p_clip_threshold, p_device, p_max_epoch, output):

    result = []

    for i in range(len(p_initLr)):

        item = []
        print(i)
        tf.reset_default_graph()
        VR.config(c_initLr=p_initLr[i], c_clip_threshold=p_clip_threshold[i], c_device=p_device,
                  c_max_epoch=p_max_epoch)

        accuracy_list, variance_list = VR.main()
        item.append(p_initLr[i])
        item.append(p_clip_threshold[i])
        item.append(p_device)
        item.append(accuracy_list)
        item.append(variance_list)

        result.append(item)

    output.put(result)

# def fnc1(x, y, output):
#     z = x + y
#     output.put(z)


def fnctot():

    init_lr = get_parameter(1e-4, 1e-3, (distribution_num, test_num))
    init_lr = init_lr.tolist()
    clip_threshold = get_parameter(11, 12, (distribution_num, test_num))  # 0.01 ~ 1
    clip_threshold = clip_threshold.tolist()
    queue_list = []
    procs = []

    for i in range(distribution_num):
        queue_list.append(Queue())  # define queues for saving the outputs of functions
        if i < gpu_0_append:
            procs.append(Process(target=vad_fnc, args=(init_lr[i], clip_threshold[i], '/gpu:0', max_epoch, queue_list[i]))) # define process
        else:
            procs.append(Process(target=vad_fnc, args=(init_lr[i], clip_threshold[i], '/gpu:1', max_epoch, queue_list[i])))

    for p in procs:  # process start
        p.start()

    M_list = []

    for i in range(distribution_num):  # save results from queues and close queues
        M_list.append(queue_list[i].get())
        queue_list[i].close()

    for p in procs:  # close process
        p.join()

    return M_list


if __name__ == "__main__":
    result = fnctot()
    print("********************END********************")
    with open('/home/sbie/storage2/result/result_proposed_soft.p', 'wb') as f:
        pickle.dump(result, f)
