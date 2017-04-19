import tensorflow as tf
import numpy as np
import utils as utils
import re
import data_reader_bDNN as dr
import os
from tensorflow.contrib import rnn

FLAGS = tf.flags.FLAGS

tf.flags.DEFINE_string('mode', "train", "mode : train/ test [default : train]")

file_dir = "/home/sbie/storage/VAD_Database/SE_TIMIT_MRCG_0328"
input_dir = file_dir
output_dir = file_dir + "/Labels"

valid_file_dir = "/home/sbie/storage/VAD_Database/NX_TIMIT_MRCG_small"
valid_input_dir = valid_file_dir + "/Babble"
valid_output_dir = valid_file_dir + "/Babble/Labels"

norm_dir = input_dir

logs_dir = "/home/sbie/github/VAD_bDNN_baseline/logs_proposed"

reset = True  # remove all existed logs and initialize log directories
device = '/gpu:0'

if FLAGS.mode is 'test':
    reset = False

if reset:

    os.popen('rm -rf ' + logs_dir + '/*')
    os.popen('mkdir ' + logs_dir + '/train')
    os.popen('mkdir ' + logs_dir + '/valid')

summary_list = ["cost", "accuracy_SNR_-5", "accuracy_SNR_0", "accuracy_SNR_5", "accuracy_SNR_10",
                "accuracy_across_all_SNRs"]

SMALL_NUM = 1e-5
eval_num_batches = 2e5
max_epoch = int(1e6)
dropout_rate = 0.5
decay = 0.9  # batch normalization decay factor
w = 19  # w default = 19
u = 9  # u default = 9
eval_th = 0.6
th = 0.6
action_hidden_1 = 512
action_hidden_2 = 512

effective_batch_size = int(4096 * 2)  # default = 4096 * 2
train_batch_size = effective_batch_size + 2*w
valid_batch_size = train_batch_size

assert (w-1) % u == 0, "w-1 must be divisible by u"

num_features = 768  # for MRCG feature = 768
bdnn_winlen = (((w-1) / u) * 2) + 3

bdnn_inputsize = int(bdnn_winlen * num_features)
bdnn_outputsize = int(bdnn_winlen)

glimpse_hidden = 256
bp_hidden = 256
glimpse_out = bp_out = 256
nGlimpses = 6
lstm_cell_size = 256

global_step = tf.Variable(0, trainable=False)

initLr = 5e-4  # default : 1e-5
lrDecayRate = .95
lrDecayFreq = 200
lr = tf.train.exponential_decay(initLr, global_step, lrDecayFreq, lrDecayRate, staircase=True)


def affine_transform(x, output_dim, name=None):
    """
    affine transformation Wx+b
    assumes x.shape = (batch_size, num_features)
    """

    w = tf.get_variable(name + "_w", [x.get_shape()[1], output_dim], initializer=tf.truncated_normal_initializer(stddev=0.02))
    b = tf.get_variable(name + "_b", [output_dim], initializer=tf.constant_initializer(0.0))

    return tf.matmul(x, w) + b


def sw_sensor(inputs, bp):

    bp = tf.expand_dims(bp, axis=2)
    bp = tf.tile(bp, (1, 1, num_features))
    bp = tf.reshape(bp, (inputs.get_shape()[0].value, -1, 1))
    bp = tf.squeeze(bp)

    sw = bp * inputs

    return sw


def get_glimpse(inputs, bp, reuse=None, is_training=True):

    with tf.variable_scope("glimpse_net", reuse=reuse):

        glimpse_input = sw_sensor(inputs, bp)

        act_glimpse_hidden = tf.nn.relu(utils.batch_norm_affine_transform(glimpse_input, glimpse_hidden, decay=decay,
                                                                          name='glimpse_hidden', is_training=is_training))
        act_bp_hidden = tf.nn.relu(utils.batch_norm_affine_transform(bp, bp_hidden, decay=decay, name='bp_hidden',
                                                                     is_training=is_training))

        glimpse_feature = tf.nn.relu(utils.batch_norm_affine_transform(act_glimpse_hidden, glimpse_out, decay=decay,
                                                                       name='glimpse_out', is_training=is_training) +
                                     utils.batch_norm_affine_transform(act_bp_hidden, bp_out, decay=decay,
                                                                       name='bp_out', is_training=is_training))

    return glimpse_feature


def bernoulli_pmf(mean, sample):
    """
    calculate the probability of bernoulli process
    :param mean: mean. shape = (batch_size, num_sbs)
    :param sample: sample. shape = (batch_size, num_sbs)
    :return: p_br: shape = (batch_size, num_sbs)
    """

    p_br = tf.pow(mean, sample) * tf.pow(1 - mean, 1 - sample)
    return p_br


def train(loss_val, var_list):

    lrDecayRate = .96
    lrDecayFreq = 200
    momentumValue = .9

    global_step = tf.Variable(0, trainable=False)
    lr = tf.train.exponential_decay(learning_rate, global_step, lrDecayFreq, lrDecayRate, staircase=True)

    # define the optimizer
    # optimizer = tf.train.MomentumOptimizer(lr, momentumValue)
    # optimizer = tf.train.AdagradOptimizer(learning_rate)
    #
    optimizer = tf.train.AdamOptimizer(lr)
    grads = optimizer.compute_gradients(loss_val, var_list=var_list)

    return optimizer.apply_gradients(grads, global_step=global_step)


def bdnn_prediction(bdnn_batch_size, logits, threshold=th):

    result = np.zeros((int(bdnn_batch_size), 1))
    indx = np.arange(int(bdnn_batch_size)) + 1
    indx = indx.reshape((int(bdnn_batch_size), 1))
    indx = utils.bdnn_transform(indx, w, u)
    indx = indx[w:(int(bdnn_batch_size)-w), :]
    indx_list = np.arange(w, int(bdnn_batch_size) - w)

    for i in indx_list:
        indx_temp = np.where((indx-1) == i)
        pred = logits[indx_temp]
        pred = np.sum(pred)/pred.shape[0]
        result[i] = pred

    result = np.trim_zeros(result)
    result = result >= threshold

    return result.astype(np.float32)


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

    return summary_dic


def full_evaluation(m_eval, sess_eval, batch_size_eval, eval_file_dir, summary_writer, summary_dic, itr):

    mean_cost = []
    mean_accuracy = []

    print("-------- Performance for each of noise types --------")

    noise_list = os.listdir(eval_file_dir)
    noise_list = sorted(noise_list)

    summary_ph = summary_dic["summary_ph"]

    for i in range(len(noise_list)):

        noise_name = '/' + noise_list[i]
        eval_input_dir = eval_file_dir + noise_name
        eval_output_dir = eval_file_dir + noise_name + '/Labels'
        eval_data_set = dr.DataReader(eval_input_dir, eval_output_dir, norm_dir, w=w, u=u, name="eval")

        eval_cost, eval_accuracy, eval_list = evaluation(m_eval, eval_data_set, sess_eval, batch_size_eval)

        print("noise type : " + noise_list[i])
        print("cost: %.3f, accuracy across all SNRs: %.3f" % (eval_cost, eval_accuracy))

        print('accuracy wrt SNR:')

        print('SNR_-5 : %.3f, SNR_0 : %.3f, SNR_5 : %.3f, SNR_10 : %.3f' % (eval_list[0], eval_list[1],
                                                                            eval_list[2], eval_list[3]))
        eval_summary_list = [eval_cost] + eval_list + [eval_accuracy]

        for j, summary_name in enumerate(summary_list):
            summary_str = sess_eval.run(summary_dic[noise_list[i]+"_"+summary_name], feed_dict={summary_ph: eval_summary_list[j]})
            summary_writer.add_summary(summary_str, itr)

        mean_cost.append(eval_cost)
        mean_accuracy.append(eval_accuracy)

    mean_cost = np.mean(np.asarray(mean_cost))
    mean_accuracy = np.mean(np.asarray(mean_accuracy))

    summary_writer.add_summary(sess_eval.run(summary_dic["cost_across_all_noise_types"],
                                             feed_dict={summary_ph: mean_cost}), itr)
    summary_writer.add_summary(sess_eval.run(summary_dic["accuracy_across_all_noise_types"],
                                             feed_dict={summary_ph: mean_accuracy}), itr)

    print("-------- Performance across all of noise types --------")
    print("cost : %.3f" % mean_cost)
    print("** accuracy across all noise_types : %.3f" % mean_accuracy)


def evaluation(m_valid, valid_data_set, sess, eval_batch_size):
    # num_samples = valid_data_set.num_samples
    # num_batches = num_samples / batch_size
    avg_valid_cost = 0.
    avg_valid_accuracy = 0.
    itr_sum = 0.

    accuracy_list = [0 for i in range(valid_data_set._file_len)]
    cost_list = [0 for i in range(valid_data_set._file_len)]
    itr_file = 0

    while True:

        valid_inputs, valid_labels = valid_data_set.next_batch(eval_batch_size)

        if valid_data_set.file_change_checker():

            accuracy_list[itr_file] = avg_valid_accuracy / itr_sum
            cost_list[itr_file] = avg_valid_cost / itr_sum
            avg_valid_accuracy = 0.
            avg_valid_cost = 0.
            itr_sum = 0
            itr_file += 1
            valid_data_set.file_change_initialize()

        if valid_data_set.eof_checker():
            valid_data_set.reader_initialize()
            print('Valid data reader was initialized!')  # initialize eof flag & num_file & start index
            break

        feed_dict = {m_valid.inputs: valid_inputs, m_valid.labels: valid_labels,
                     m_valid.keep_probability: 1}

        valid_cost, valid_accuracy = sess.run([m_valid.cost, m_valid.reward], feed_dict=feed_dict)

        avg_valid_cost += valid_cost
        avg_valid_accuracy += valid_accuracy
        itr_sum += 1

    total_avg_valid_cost = np.asscalar(np.mean(np.asarray(cost_list)))
    total_avg_valid_accuracy = np.asscalar(np.mean(np.asarray(accuracy_list)))

    return total_avg_valid_cost, total_avg_valid_accuracy, accuracy_list


def dense_to_one_hot(labels_dense, num_classes=2):

    """Convert class labels from scalars to one-hot vectors."""
    # copied from TensorFlow tutorial
    num_labels = labels_dense.shape[0]
    index_offset = np.arange(num_labels) * num_classes
    labels_one_hot = np.zeros((num_labels, num_classes))
    labels_one_hot.flat[(index_offset + labels_dense.ravel()).astype(int)] = 1
    return labels_one_hot.astype(np.float32)


class Model(object):

    def __init__(self, batch_size, reuse=None, is_training=True):

        self.batch_size = batch_size
        self.keep_probability = tf.placeholder(tf.float32, name="keep_probabilty")
        self.inputs = tf.placeholder(tf.float32, shape=[batch_size, bdnn_inputsize],
                                              name="inputs")
        self.labels = tf.placeholder(tf.float32, shape=[batch_size, bdnn_outputsize], name="labels")
        self.is_training = is_training
        self.mean_bps = []
        self.sampled_bps = []
        self.baselines = []

        # set inference graph
        cell_outputs = self.inference(reuse)  # (batch_size, bdnn_outputsize)
        # set objective function

        self.cost, self.reward, self.train_op, self.avg_b, self.rminusb, self.sampled_bps_tensor, self.p_bps, self.lr\
            = self.calc_reward(cell_outputs)

    def inference(self, reuse=None):

        # initialization
        raw_inputs = self.inputs
        batch_size = self.batch_size
        keep_prob = self.keep_probability
        is_training = self.is_training

        lstm_cell = rnn.LayerNormBasicLSTMCell(lstm_cell_size, dropout_keep_prob=keep_prob, reuse=reuse)
        initial_state = lstm_cell.zero_state(batch_size, tf.float32)

        tf.set_random_seed(1)  # initialize the random seed at graph level

        init_sw = tf.random_uniform([batch_size, int(bdnn_winlen)], minval=0, maxval=1)
        self.mean_bps.append(init_sw)

        init_sw = tf.cast(tf.greater(init_sw, 0.5), tf.float32)
        self.sampled_bps.append(init_sw)

        reuse_recurrent = None

        init_glimpse = self.get_glimpse(raw_inputs, init_sw, reuse=reuse_recurrent)  # (batch_size, glimpse_out)

        inputs = [0] * nGlimpses
        outputs = [0] * nGlimpses
        glimpse = init_glimpse

        for time_step in range(nGlimpses):

            if time_step == 0:
                with tf.variable_scope("core_network", reuse=reuse_recurrent):
                    (cell_output, cell_state) = lstm_cell(glimpse, initial_state)
            else:
                reuse_recurrent = True
                with tf.variable_scope("core_network", reuse=reuse_recurrent):
                    (cell_output, cell_state) = lstm_cell(glimpse, cell_state)

            inputs[time_step] = glimpse
            outputs[time_step] = cell_output

            if time_step != nGlimpses - 1:  # not final time_step

                glimpse = self.get_next_input(cell_output, reuse=reuse_recurrent)

            else:  # final time_step
                with tf.variable_scope("baseline", reuse=reuse_recurrent):
                    cell_output_no_grad = tf.stop_gradient(cell_output)
                    # baseline = tf.sigmoid(
                    #     utils.batch_norm_affine_transform(cell_output_no_grad, 1, decay=decay, name='baseline',
                    #                                       is_training=is_training))
                    # baseline = tf.sigmoid(affine_transform(cell_output_no_grad, 1, name='baseline'))
                    baseline = utils.clipped_relu(affine_transform(cell_output_no_grad, 1, name='baseline'))

                    self.baselines.append(baseline)

        return outputs

    def get_glimpse(self, inputs, bp, reuse=None):

        is_training = self.is_training

        with tf.variable_scope("glimpse_net", reuse=reuse):
            glimpse_input = sw_sensor(inputs, bp)

            act_glimpse_hidden = tf.nn.relu(
                utils.batch_norm_affine_transform(glimpse_input, glimpse_hidden, decay=decay,
                                                  name='glimpse_hidden', is_training=is_training))
            act_bp_hidden = tf.nn.relu(utils.batch_norm_affine_transform(bp, bp_hidden, decay=decay, name='bp_hidden',
                                                                         is_training=is_training))

            glimpse_feature = tf.nn.relu(utils.batch_norm_affine_transform(act_glimpse_hidden, glimpse_out, decay=decay,
                                                                           name='glimpse_out',
                                                                           is_training=is_training) +
                                         utils.batch_norm_affine_transform(act_bp_hidden, bp_out, decay=decay,
                                                                           name='bp_out', is_training=is_training))
        return glimpse_feature

    def get_next_input(self, cell_output, reuse=None):

        raw_inputs = self.inputs
        is_training = self.is_training
        cell_output_no_grad = tf.stop_gradient(cell_output)

        with tf.variable_scope("baseline", reuse=reuse):

            # baseline = tf.sigmoid(
            #     utils.batch_norm_affine_transform(cell_output_no_grad, 1, decay=decay, name='baseline',
            #                                       is_training=is_training))
            # baseline = tf.sigmoid(affine_transform(cell_output_no_grad, 1, name='baseline'))
            baseline = utils.clipped_relu(affine_transform(cell_output_no_grad, 1, name='baseline'))

            self.baselines.append(baseline)

        with tf.variable_scope("selection_network", reuse=reuse):
            # mean_bp = tf.sigmoid(
            #     utils.batch_norm_affine_transform(cell_output_no_grad, int(bdnn_winlen), decay=decay, name='selection',
            #                                       is_training=is_training))
            mean_bp = utils.clipped_relu(
                utils.batch_norm_affine_transform(cell_output_no_grad, int(bdnn_winlen), decay=decay, name='selection',
                                                  is_training=is_training))
            self.mean_bps.append(mean_bp)

            rand_seq = tf.random_uniform(mean_bp.get_shape().as_list(), minval=0, maxval=1)
            sampled_bp = tf.cast(tf.greater(mean_bp, rand_seq), tf.float32)
            sampled_bp = tf.stop_gradient(sampled_bp)
            self.sampled_bps.append(sampled_bp)

        return get_glimpse(raw_inputs, sampled_bp, is_training)

    def action_network(self, outputs):
        is_training = self.is_training
        with tf.variable_scope("action_network"):

            h1_out = tf.nn.relu(utils.batch_norm_affine_transform(outputs, action_hidden_1,
                                                                  decay=decay, name='action_hidden_1', is_training=is_training))
            h1_out = tf.nn.dropout(h1_out, keep_prob=self.keep_probability)
            h2_out = tf.nn.relu(utils.batch_norm_affine_transform(h1_out, action_hidden_2,
                                                                  decay=decay, name='action_hidden_2', is_training=is_training))
            h2_out = tf.nn.dropout(h2_out, keep_prob=self.keep_probability)

        return h2_out

    def bdnn_prediction(self, logits, threshold=0.5):

        # batch_size = self.batch_size
        # result = tf.zeros(shape=(batch_size, 1), dtype=tf.float32)
        # indx = np.arange(batch_size) + 1
        # indx = indx.reshape((batch_size, 1))
        # indx = utils.bdnn_transform(indx, w, u)
        # indx = indx[w:(batch_size - w), :]
        # indx_list = np.arange(w, batch_size - w)
        #
        # for i in indx_list:
        #     indx_temp = np.where((indx - 1) == i)
        #     pred = logits[int(indx_temp)]
        #     pred = tf.reduce_mean(pred)
        #     result[i] = pred
        #
        # result = tf.py_func(self.np_trim_zeros, [result], tf.float32)
        # result = tf.greater(result, threshold)
        # result = tf.cast(result, tf.float32)

        batch_size_tensor = tf.constant(self.batch_size+2*w, dtype=tf.float32)
        th_tenor = tf.constant(threshold, dtype=tf.float32)

        result = tf.py_func(bdnn_prediction, [batch_size_tensor, logits, th_tenor], tf.float32)

        return result

    @staticmethod
    def np_trim_zeros(x):
        return np.trim_zeros(x)

    def calc_reward(self, outputs):

        batch_size = self.batch_size

        # consider the action at the last time step

        outputs = outputs[-1]
        outputs = tf.reshape(outputs, (batch_size, lstm_cell_size))

        # get the baseline

        b = tf.stack(self.baselines)
        b = tf.tile(b, [1, 1, int(bdnn_winlen)])
        b = tf.reshape(tf.transpose(b, [1, 0, 2]), [batch_size, nGlimpses * int(bdnn_winlen)])
        no_grad_b = tf.stop_gradient(b)

        # get the action

        action_out = self.action_network(outputs)
        logits = tf.sigmoid(affine_transform(action_out, int(bdnn_outputsize), name="softmax"))
        result = self.bdnn_prediction(logits, threshold=0.5)

        # convert list of tensors to one big tensor

        mean_bps = tf.concat(axis=0, values=self.mean_bps)
        mean_bps = tf.reshape(mean_bps, (nGlimpses, self.batch_size, int(bdnn_winlen)))
        mean_bps = tf.transpose(mean_bps, [1, 0, 2])

        sampled_bps = tf.concat(axis=0, values=self.sampled_bps)
        sampled_bps = tf.reshape(sampled_bps, (nGlimpses, self.batch_size, int(bdnn_winlen)))
        sampled_bps = tf.transpose(sampled_bps, [1, 0, 2])

        # reward for all examples in the batch

        raw_indx = int(np.floor(bdnn_outputsize / 2))
        raw_labels = self.labels[:, raw_indx]
        raw_labels = tf.reshape(raw_labels, shape=(-1, 1))
        R = tf.cast(tf.equal(result, raw_labels), tf.float32)
        R = tf.stop_gradient(R)
        R = tf.reshape(R, (batch_size, 1))
        R = tf.tile(R, [1, nGlimpses * int(bdnn_winlen)])
        reward = tf.reduce_mean(R)

        # select the window

        p_bps = bernoulli_pmf(mean_bps, sampled_bps)
        p_bps = tf.reshape(p_bps, (self.batch_size, nGlimpses * int(bdnn_winlen)))

        # define the cost function

        sv_part = -tf.square(self.labels - logits)
        rf_part = tf.log(p_bps + SMALL_NUM) * (R - no_grad_b)

        J = tf.concat(axis=1, values=[sv_part, rf_part])
        J = tf.reduce_sum(J, 1)
        J = J - tf.reduce_mean(tf.square(R - b), 1)
        J = tf.reduce_mean(J, 0)
        cost = -J

        var_list = tf.trainable_variables()
        grads = tf.gradients(cost, var_list)
        grads, _ = tf.clip_by_global_norm(grads, 0.5)
        optimizer = tf.train.AdamOptimizer(lr)
        train_op = optimizer.apply_gradients(zip(grads, var_list), global_step=global_step)

        return cost, reward, train_op, tf.reduce_mean(b), tf.reduce_mean(R - b), sampled_bps, tf.reduce_mean(p_bps), lr


def main(argv=None):

    #                               Graph Part                                 #

    print("Graph initialization...")
    with tf.device(device):
        with tf.variable_scope("model", reuse=None):
            m_train = Model(batch_size=effective_batch_size, reuse=None, is_training=True)
        with tf.variable_scope("model", reuse=True):
            m_valid = Model(batch_size=effective_batch_size, reuse=True, is_training=False)

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
        saver.restore(sess, ckpt.model_checkpoint_path)
        print("Done")
    else:
        sess.run(tf.global_variables_initializer())  # if the checkpoint doesn't exist, do initialization

    train_data_set = dr.DataReader(input_dir, output_dir, norm_dir, w=w, u=u, name="train")  # training data reader initialization

    if FLAGS.mode is 'train':

        for itr in range(max_epoch):

            train_inputs, train_labels = train_data_set.next_batch(train_batch_size)

            feed_dict = {m_train.inputs: train_inputs, m_train.labels: train_labels,
                         m_train.keep_probability: dropout_rate}

            sess.run(m_train.train_op, feed_dict=feed_dict)

            if itr % 10 == 0 and itr > 0:

                train_cost, train_reward, train_avg_b, train_rminusb, train_p_bps, train_lr \
                    = sess.run([m_train.cost, m_train.reward, m_train.avg_b, m_train.rminusb, m_train.p_bps, m_train.lr]
                               , feed_dict=feed_dict)
                sampled_bps_tensor = sess.run(m_train.sampled_bps_tensor, feed_dict=feed_dict)
                sampled_bps_tensor = np.mean(sampled_bps_tensor[:, -1, :], axis=1)

                print("Step: %d, cost: %.3f, accuracy: %3.3f, b: %3.3f, R-b: %3.3f, p_bps: %6.6f, lr: %7.6f" %
                      (itr, train_cost, train_reward, train_avg_b, train_rminusb, train_p_bps, train_lr))

                train_cost_summary_str = sess.run(cost_summary_op, feed_dict={summary_ph: train_cost})
                train_accuracy_summary_str = sess.run(accuracy_summary_op, feed_dict={summary_ph: train_reward})
                train_summary_writer.add_summary(train_cost_summary_str, itr)  # write the train phase summary to event files
                train_summary_writer.add_summary(train_accuracy_summary_str, itr)

            # if train_data_set.eof_checker():

            if itr % 100 == 0 and itr > 0:

                saver.save(sess, logs_dir + "/model.ckpt", itr)  # model save
                print('validation start!')

                full_evaluation(m_valid, sess, valid_batch_size, valid_file_dir, valid_summary_writer, summary_dic, itr)

                # train_data_set.initialize()
    elif FLAGS.mode is 'test':
        pass

if __name__ == "__main__":
    tf.app.run()


