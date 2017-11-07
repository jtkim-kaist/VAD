import tensorflow as tf
import numpy as np
import utils_jskim as utils
import re
import data_reader_RNN as dr
import os
import matplotlib.pyplot as plt
from tensorflow.contrib import rnn
from sklearn import metrics
import time

FLAGS = tf.flags.FLAGS
SEED = 1
tf.set_random_seed(SEED)
tf.flags.DEFINE_string('mode', "train", "mode : train/ test [default : train]")

file_dir = "/home/sbie/storage2/VAD_Database/SE_TIMIT_MRCG_0328"
input_dir = file_dir
output_dir = file_dir + "/Labels"
valid_file_dir = "/home/sbie/storage2/VAD_Database/NX_TIMIT_MRCG_big"

# valid_file_dir = "/media/jskim/F440795840792312/database/jskim/NX_TIMIT_MRCG_big"

norm_dir = input_dir

logs_dir = "/home/sbie/github/VAD_Project_test/VAD_LSTM/logs_LSTM"
initial_logs_dir = "/home/sbie/github/VAD_Project_test/VAD_LSTM/logs_LSTM"
ckpt_name = "/LSTM"

reset = True  # remove all existed logs and initialize log directories
device = "/gpu:0"
mode = 'test'
if mode is 'test':
    reset = False

if reset:

    os.popen('rm -rf ' + logs_dir + '/*')
    os.popen('mkdir ' + logs_dir + '/train')
    os.popen('mkdir ' + logs_dir + '/valid')

summary_list = ["cost", "accuracy_SNR_-5", "accuracy_SNR_0", "accuracy_SNR_5", "accuracy_SNR_10",
                "accuracy_across_all_SNRs"]

learning_rate = 0.0001
eval_num_batches = 2e4
SMALL_NUM = 1e-4
max_epoch = int(1e5)
dropout_rate = 0.5

decay = 0.9  # batch normalization decay factor
w = 5  # w default = 19
u = 9  # u default = 9
eval_th = 0.6
th = 0.5
lstm_cell_size = 256
num_hidden_1 = 128
num_hidden_2 = 256
num_hidden_3 = 128

seq_size = 20
batch_num = 200
batch_size = batch_num*seq_size # batch_size = 32
valid_batch_size = batch_size

# assert (w-1) % u == 0, "w-1 must be divisible by u"

width = 768
num_features = 768  # MRCG feature
bdnn_winlen = (((w-1) / u) * 2) + 3

# bdnn_inputsize = int(bdnn_winlen * num_features)
bdnn_inputsize = num_features
bdnn_outputsize = 2#int(bdnn_winlen)
initLr = 1e-5
scope_name = 'RNN_scope'
eval_type = 2


def test_config(c_test_dir, c_norm_dir, c_initial_logs_dir, c_batch_size_eval, c_data_len):

    global test_file_dir
    global norm_dir
    global initial_logs_dir
    global ckpt_name
    global valid_batch_size
    global batch_size
    global data_len
    global batch_num
    test_file_dir = c_test_dir
    norm_dir = c_norm_dir
    initial_logs_dir = c_initial_logs_dir
    valid_batch_size = batch_size = c_batch_size_eval
    data_len = c_data_len
    batch_num = batch_size/seq_size

def affine_transform(x, output_dim, name=None):
    """
    affine transformation Wx+b
    assumes x.shape = (batch_size, num_features)
    """

    w = tf.get_variable(name + "_w", [x.get_shape()[1], output_dim],
                        initializer=tf.truncated_normal_initializer(stddev=0.02))
    b = tf.get_variable(name + "_b", [output_dim], initializer=tf.constant_initializer(0.0))

    return tf.matmul(x, w) + b


def rnn_in(inputs, batch_num, seq_size,delay):

    batch_num = int(batch_num)
    seq_size = int(seq_size)
    delay = int(delay)
    temp1 = tf.reshape(inputs,[-1,num_features])
    temp2 = tf.reshape(temp1[0:batch_size,:],[batch_num,seq_size,num_features])
    # temp3 = tf.zeros([batch_num,delay,num_features], tf.float32)
    temp3 = temp2[1:batch_num,0:delay,:]
    temp4 = tf.reshape(temp1[batch_size:batch_size+delay,:],[-1,delay,num_features])
    temp5 = tf.concat([temp3,temp4],0)
    return tf.concat([temp2,temp5],1)


def inference(inputs, keep_prob, is_training=True, reuse=None):

    # initialization
    # c1_out = affine_transform(inputs, num_hidden_1, name="hidden_1")
    # inputs_shape = inputs.get_shape().as_list()
    with tf.variable_scope(scope_name):
        print('inference time')
        # print(inputs.get_shape().as_list())
        in_rnn = rnn_in(inputs, batch_num, seq_size, w)
        # in_rnn = tf.reshape(inputs,[-1, seq_size+w, num_features])
        stacked_rnn = []
        for iiLyr in range(3):
            stacked_rnn.append(tf.nn.rnn_cell.LSTMCell(num_units=lstm_cell_size, state_is_tuple=True))
        MultiLyr_cell = tf.nn.rnn_cell.MultiRNNCell(cells=stacked_rnn, state_is_tuple=True)

        outputs, _state = tf.nn.dynamic_rnn(MultiLyr_cell, in_rnn, time_major=False, dtype=tf.float32)
        outputs = tf.reshape(outputs[:, 0:seq_size, :],[-1,lstm_cell_size])

        outputs = tf.nn.dropout(outputs, keep_prob=keep_prob)

        # # h1_out = affine_transform(inputs, num_hidden_1, name="hidden_1")
        # lh1_out = utils.batch_norm_affine_transform(outputs, num_hidden_1, name="lhidden_1", decay=decay,
        #                                             is_training=is_training)
        # lh1_out = tf.nn.relu(lh1_out)
        # lh1_out = tf.nn.dropout(lh1_out, keep_prob=keep_prob)

        logits = affine_transform(outputs, bdnn_outputsize, name="output1")
        # logits = tf.sigmoid(logits)
        logits = tf.reshape(logits, [-1, int(bdnn_outputsize)])
        print(logits.get_shape().as_list())
        print('asdf')
    return logits


def train(loss_val, var_list):

    lrDecayRate = .95
    lrDecayFreq = 200

    global_step = tf.Variable(0, trainable=False)
    lr = tf.train.exponential_decay(initLr, global_step, lrDecayFreq, lrDecayRate, staircase=True)

    optimizer = tf.train.AdamOptimizer(lr)
    grads = optimizer.compute_gradients(loss_val, var_list=var_list)

    return optimizer.apply_gradients(grads, global_step=global_step)


def bdnn_prediction(bdnn_batch_size, logits, threshold=th):

    result = np.zeros((bdnn_batch_size, 1))
    indx = np.arange(bdnn_batch_size) + 1
    indx = indx.reshape((bdnn_batch_size, 1))
    indx = utils.bdnn_transform(indx, w, u)
    indx = indx[w:(bdnn_batch_size-w), :]
    indx_list = np.arange(w, bdnn_batch_size - w)

    for i in indx_list:
        indx_temp = np.where((indx-1) == i)
        pred = logits[indx_temp]
        pred = np.sum(pred)/pred.shape[0]
        result[i] = pred

    result = np.trim_zeros(result)
    result = result >= threshold

    return result.astype(int)


def summary_generation(eval_file_dir):

    summary_dic = {}

    noise_list = os.listdir(eval_file_dir)
    noise_list = sorted(noise_list)
    summary_dic["summary_ph"] = summary_ph = tf.placeholder(dtype=tf.float32)

    for name in noise_list:

        with tf.variable_scope(name):
            for summary_name in summary_list:
                    summary_dic[name+"_"+summary_name] = tf.summary.scalar(summary_name, summary_ph)

    with tf.variable_scope("Averaged_Results"):

        summary_dic["cost_across_all_noise_types"] = tf.summary.scalar("cost_across_all_noise_types", summary_ph)
        summary_dic["accuracy_across_all_noise_types"]\
            = tf.summary.scalar("accuracy_across_all_noise_types", summary_ph)
        summary_dic["variance_across_all_noise_types"]\
            = tf.summary.scalar("variance_across_all_noise_types", summary_ph)
        summary_dic["AUC_across_all_noise_types"]\
            = tf.summary.scalar("AUC_across_all_noise_types", summary_ph)
    return summary_dic


def full_evaluation(m_eval, sess_eval, batch_size_eval, eval_file_dir, summary_writer, summary_dic, itr):

    mean_cost = []
    mean_accuracy = []
    mean_auc = []
    mean_time = []

    print("-------- Performance for each of noise types --------")

    noise_list = os.listdir(eval_file_dir)
    noise_list = sorted(noise_list)

    summary_ph = summary_dic["summary_ph"]

    for i in range(len(noise_list)):
        print("full time evaluation, now loading : %d",i)
        noise_name = '/' + noise_list[i]
        eval_input_dir = eval_file_dir + noise_name
        eval_output_dir = eval_file_dir + noise_name + '/Labels'
        ##########################################


        eval_calc_dir = eval_file_dir + noise_name + '/test_result' # for Final layer information saving



        ##########################################
        eval_data_set = dr.DataReader(eval_input_dir, eval_output_dir, norm_dir, w=w, u=u, name="eval")

        eval_cost, eval_accuracy, eval_list, eval_auc, auc_list, eval_time = evaluation(m_eval, eval_data_set, sess_eval,  batch_size_eval, noise_list[i], save_dir = eval_calc_dir)

        print("--noise type : " + noise_list[i])
        print("cost: %.3f, accuracy across all SNRs: %.3f, auc across all SNRs: %.3f " % (eval_cost, eval_accuracy, eval_auc))

        print('accuracy wrt SNR:')

        print('SNR_-5 : %.3f, SNR_0 : %.3f, SNR_5 : %.3f, SNR_10 : %.3f' % (eval_list[0], eval_list[1],
                                                                            eval_list[2], eval_list[3]))

        print('AUC wrt SNR:')

        print('SNR_-5 : %.3f, SNR_0 : %.3f, SNR_5 : %.3f, SNR_10 : %.3f' % (auc_list[0], auc_list[1],
                                                                            auc_list[2], auc_list[3]))
        eval_summary_list = [eval_cost] + eval_list + [eval_accuracy] + [eval_auc]

        for j, summary_name in enumerate(summary_list):
            summary_str = sess_eval.run(summary_dic[noise_list[i]+"_"+summary_name],
                                        feed_dict={summary_ph: eval_summary_list[j]})
            summary_writer.add_summary(summary_str, itr)

        mean_cost.append(eval_cost)
        mean_accuracy.append(eval_accuracy)
        mean_auc.append(eval_auc)
        mean_time.append(eval_time)

    mean_cost = np.mean(np.asarray(mean_cost))
    var_accuracy = np.var(np.asarray(mean_accuracy))
    mean_accuracy = np.mean(np.asarray(mean_accuracy))
    mean_auc = np.mean(np.asarray(mean_auc))
    mean_time = np.mean(np.asarray(mean_time))


    summary_writer.add_summary(sess_eval.run(summary_dic["cost_across_all_noise_types"],
                                             feed_dict={summary_ph: mean_cost}), itr)
    summary_writer.add_summary(sess_eval.run(summary_dic["accuracy_across_all_noise_types"],
                                             feed_dict={summary_ph: mean_accuracy}), itr)
    summary_writer.add_summary(sess_eval.run(summary_dic["variance_across_all_noise_types"],
                                             feed_dict={summary_ph: var_accuracy}), itr)
    summary_writer.add_summary(sess_eval.run(summary_dic["AUC_across_all_noise_types"],
                                             feed_dict={summary_ph: mean_auc}), itr)

    print("-------- Performance across all of noise types --------")
    print("cost : %.3f" % mean_cost)
    print("******* averaged accuracy across all noise_types : %.3f *******" % mean_accuracy)
    print("******* variance of accuracies across all noise_types : %6.6f *******" % var_accuracy)
    print("******* variance of AUC across all noise_types : %6.6f *******" % mean_auc)
    print("******* mean time : %6.6f *******" % mean_time)
    return mean_auc


def evaluation(m_valid, valid_data_set, sess, eval_batch_size, noise_name, num_batches=eval_num_batches,
               save_dir = None):
    # num_samples = valid_data_set.num_samples
    # num_batches = num_samples / batch_size
    avg_valid_cost = 0.
    avg_valid_accuracy = 0.
    avg_valid_time = 0.
    # AUC = 0.
    itr_sum = 0.
    file_num_before = -1
    accuracy_list = [0 for i in range(valid_data_set._file_len)]
    cost_list = [0 for i in range(valid_data_set._file_len)]
    auc_list = [0 for i in range(valid_data_set._file_len)]
    time_list = [0 for i in range(valid_data_set._file_len)]
    itr_file = 0
    channel_values = []
    save_calc_dir = ''
    valid_name_before = ''
    # plt.figure()
    while True:

        # valid_name_before = save_dir + 'test'+ [file_num_before] + '.txt'
        valid_inputs, valid_labels = valid_data_set.next_batch(eval_batch_size)

        if valid_data_set.file_change_checker():

            accuracy_list[itr_file] = avg_valid_accuracy / itr_sum
            cost_list[itr_file] = avg_valid_cost / itr_sum
            auc_list[itr_file] = utils.plot_ROC2(channel_values,file_num_before, noise_name)
            time_list[itr_file] = avg_valid_time / ((itr_sum*batch_size)/16000)
            avg_valid_accuracy = 0.
            avg_valid_cost = 0.
            avg_valid_time = 0.
            channel_values = []
            # avg_valid_auc = 0.
            itr_sum = 0
            itr_file += 1
            valid_data_set.file_change_initialize()

        if valid_data_set.eof_checker():
            #######
            # AUC = utils.plot_ROC2(valid_name_before, save_calc_dir,file_num_before, noise_name)
            # f_eval.close()
            #######
            valid_data_set.reader_initialize()
            print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
            break

        # if eval_batch_size * itr_file > 5000:
        #     f_eval.close()
        #     break

        one_hot_vlabels = valid_labels.reshape((-1, 1))
        one_hot_vlabels = dense_to_one_hot(one_hot_vlabels, num_classes=2)
        feed_dict = {m_valid.inputs: valid_inputs, m_valid.labels: one_hot_vlabels,
                     m_valid.keep_probability: 1}
        start_time = time.time()
        logits_val = sess.run(m_valid.logits, feed_dict=feed_dict)
        check_time = time.time()-start_time
        valid_cost, valid_accuracy = sess.run([m_valid.cost, m_valid.accuracy], feed_dict=feed_dict)
        # print(valid_labels.shape)
        # print(logits_val)

        save_val = np.concatenate([logits_val,valid_labels], axis = 1)
        ###############################################

        file_num = valid_data_set._num_file
        channel_values.append(save_val)

        avg_valid_cost += valid_cost
        avg_valid_accuracy += valid_accuracy
        avg_valid_time += check_time
        # avg_valid_auc += valid_auc
        itr_sum += 1

        ############################################
        file_num_before = valid_data_set._num_file
        ############################################

    total_avg_valid_cost = np.asscalar(np.mean(np.asarray(cost_list)))
    total_avg_valid_accuracy = np.asscalar(np.mean(np.asarray(accuracy_list)))
    total_avg_valid_auc = np.asscalar(np.mean(np.asarray(auc_list)))
    total_avg_valid_time = np.asscalar(np.mean(np.asarray(time_list)))
    return total_avg_valid_cost, total_avg_valid_accuracy, accuracy_list, total_avg_valid_auc, auc_list, total_avg_valid_time


def dense_to_one_hot(labels_dense, num_classes=2):

    """Convert class labels from scalars to one-hot vectors."""
    # copied from TensorFlow tutorial
    num_labels = labels_dense.shape[0]
    index_offset = np.arange(num_labels) * num_classes
    labels_one_hot = np.zeros((num_labels, num_classes))
    labels_one_hot.flat[(index_offset + labels_dense.ravel()).astype(int)] = 1
    return labels_one_hot.astype(np.float32)


class Model(object):
    def __init__(self, is_training=True):
        self.keep_probability = tf.placeholder(tf.float32, name="keep_probabilty")
        self.inputs = inputs = tf.placeholder(tf.float32, shape=[batch_size + w, bdnn_inputsize],
                                              name="inputs")
        self.labels = labels = tf.placeholder(tf.float32, shape=[batch_size, bdnn_outputsize], name="labels")

        # set inference graph
        self.logits = logits = inference(inputs, self.keep_probability, is_training=is_training)  # (batch_size, bdnn_outputsize)
        # set objective function
        # self.cost = cost = tf.reduce_mean(tf.nn.sigmoid_cross_entropy_with_logits(labels=labels, logits=logits))
        pred = tf.argmax(logits, axis=1, name="prediction")
        pred = tf.cast(pred, tf.int32)
        truth = tf.cast(labels[:, 1], tf.int32)
        log_one = logits[:, 1]

        self.accuracy = tf.reduce_mean(tf.cast(tf.equal(pred, truth), tf.float32))
        self.cost = cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(labels=labels, logits = logits))
        # self.cost = cost = tf.reduce_mean(tf.square(labels - logits))
        # fpr, tpr, thresholds = metrics.roc_curve(np.array(truth), np.array(pred), pos_label=2)
        # self.auc = metrics.auc(fpr, tpr)
        # cost = tf.nn.softmax_cross_entropy_with_logits(labels=labels, logits=logits)
        # cost = tf.reduce_sum(tf.square(labels - logits), axis=1)
        # self.cost = cost = tf.reduce_mean(cost)

        # self.sigm = tf.sigmoid(logits)
        # set training strategy

        trainable_var = tf.trainable_variables()
        self.train_op = train(cost, trainable_var)


def main(argv=None):
    mean_acc = 0
    #                               Graph Part                                 #

    print("Graph initialization...")
    with tf.device(device):
        with tf.variable_scope("model", reuse=None):
            m_train = Model(is_training=True)
        with tf.variable_scope("model", reuse=True):
            m_valid = Model(is_training=False)

    print("Done")

    #                               Summary Part                               #

    print("Setting up summary op...")
    summary_ph = tf.placeholder(dtype=tf.float32)

    with tf.variable_scope("Aggregated_Training_Parts"):

        cost_summary_op = tf.summary.scalar("cost", summary_ph)
        accuracy_summary_op = tf.summary.scalar("accuracy", summary_ph)

    train_summary_writer = tf.summary.FileWriter(logs_dir + '/train/', max_queue=4)
    valid_summary_writer = tf.summary.FileWriter(logs_dir + '/valid/', max_queue=4)
    summary_dic = summary_generation(valid_file_dir)

    print("Done")

    #                               Model Save Part                            #

    print("Setting up Saver...")
    saver = tf.train.Saver()
    ckpt = tf.train.get_checkpoint_state(logs_dir)
    print("Done")

    #                               Session Part                               #

    sess_config = tf.ConfigProto(allow_soft_placement=True, log_device_placement=False)
    sess_config.gpu_options.allow_growth = True
    sess = tf.Session(config=sess_config)

    if ckpt and ckpt.model_checkpoint_path:  # model restore
        print("Model restored...")
        saver.restore(sess, initial_logs_dir+ckpt_name)
        print("Done")
    else:
        sess.run(tf.global_variables_initializer())  # if the checkpoint doesn't exist, do initialization

    train_data_set = dr.DataReader(input_dir, output_dir, norm_dir, w=w, u=u, name="train")  # training data reader initialization

    if mode is 'train':

        for itr in range(max_epoch):

            train_inputs, train_labels = train_data_set.next_batch(batch_size)
            one_hot_labels = train_labels.reshape((-1, 1))
            one_hot_labels = dense_to_one_hot(one_hot_labels, num_classes=2)
            feed_dict = {m_train.inputs: train_inputs, m_train.labels: one_hot_labels,
                         m_train.keep_probability: dropout_rate}

            sess.run(m_train.train_op, feed_dict=feed_dict)

            if itr % 50 == 0 and itr >= 0:
                train_cost, train_accuracy = sess.run([m_train.cost, m_train.accuracy], feed_dict=feed_dict)

                # train_cost, logits = sess.run([m_train.cost, m_train.logits], feed_dict=feed_dict)
                #
                # result = bdnn_prediction(batch_size, logits, threshold=th)
                # raw_indx = int(np.floor(train_labels.shape[1] / 2))
                # raw_labels = train_labels[:, raw_indx]
                # raw_labels = raw_labels.reshape((-1, 1))
                # train_accuracy = np.equal(result, raw_labels)
                # train_accuracy = train_accuracy.astype(int)
                # train_accuracy = np.sum(train_accuracy) / batch_size  # change to mean...

                print("Step: %d, train_cost: %.3f, train_accuracy=%3.3f" % (itr, train_cost, train_accuracy))

                train_cost_summary_str = sess.run(cost_summary_op, feed_dict={summary_ph: train_cost})
                train_accuracy_summary_str = sess.run(accuracy_summary_op, feed_dict={summary_ph: train_accuracy})
                train_summary_writer.add_summary(train_cost_summary_str, itr)  # write the train phase summary to event files
                train_summary_writer.add_summary(train_accuracy_summary_str, itr)

            # if train_data_set.eof_checker():
            if itr % 50 == 0 and itr > 0:

                saver.save(sess, logs_dir + "/model.ckpt", itr)  # model save
                print('validation start!')
                mean_acc = full_evaluation(m_valid, sess, valid_batch_size, valid_file_dir,
                                           valid_summary_writer, summary_dic, itr)
            if mean_acc > 0.968:
                print('finish!!')
                break
                # train_data_set.reader_initialize()
                # print('Train data reader was initialized!')  # initialize eof flag & num_file & start index

    elif mode is 'test':

        final_softout, final_label = utils.vad_test3(m_valid, sess, valid_batch_size, test_file_dir, norm_dir, data_len,

                                                     eval_type)

    if data_len is None:

        return final_softout, final_label

    else:

        return final_softout[0:data_len, :], final_label[0:data_len, :]
if __name__ == "__main__":
    tf.app.run()
